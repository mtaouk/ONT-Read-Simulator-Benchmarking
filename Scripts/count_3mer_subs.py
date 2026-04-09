#!/usr/bin/env python3
"""
This script takes in read to reference alignments (along with the read and reference files), and it
outputs all counts for 3-mer substitutions (1-bp substitutions between two matching bases).
"""

import argparse
import collections
import gzip
import re

CANONICAL_BASES = {'A', 'C', 'G', 'T'}


def get_arguments():
    parser = argparse.ArgumentParser(description='Count 3-mer substitutions')

    parser.add_argument('fastq', type=str,
                        help='Reads in FASTQ format')
    parser.add_argument('reference', type=str,
                        help='Reference in FASTA format')
    parser.add_argument('paf', type=str,
                        help='Alignments in PAF format')

    return parser.parse_args()


def main():
    args = get_arguments()
    reads = load_fastq_as_dict(args.fastq)
    ref = load_fasta_as_dict(args.reference)
    alignments = keep_only_best_alignment_per_read(load_alignments(args.paf))

    sub_counts = collections.defaultdict(int)
    total_counts = collections.defaultdict(int)

    for a in alignments:
        read_seq = reads[a.query_name][a.query_start:a.query_end]
        ref_seq = ref[a.target_name][a.target_start:a.target_end]
        if a.strand == '-':
            read_seq = reverse_complement(read_seq)

        for i in range(len(ref_seq) - 2):
            ref_3mer = ref_seq[i:i + 3]
            if set(ref_3mer) - CANONICAL_BASES:
                continue
            total_counts[canonicalise_3mer(ref_3mer)] += 1

        expanded_cigar = get_expanded_cigar(a.cigar)
        aligned_read_seq, aligned_ref_seq = align_sequences(read_seq, ref_seq, expanded_cigar)
        assert len(aligned_read_seq) == len(aligned_ref_seq) == len(expanded_cigar)

        for i in range(1, len(expanded_cigar) - 1):
            cigar_3mer = expanded_cigar[i - 1:i + 2]
            if cigar_3mer != '=X=':
                continue

            ref_3mer = aligned_ref_seq[i - 1:i + 2]
            read_3mer = aligned_read_seq[i - 1:i + 2]
            if set(ref_3mer + read_3mer) - CANONICAL_BASES:
                continue

            ref_3mer, read_3mer = canonicalise_3mer_substitution(ref_3mer, read_3mer)
            sub_counts[(ref_3mer, read_3mer)] += 1

    print('ref_seq\tread_seq\tcount\terror_rate')
    for (ref_3mer, read_3mer), count in sorted(sub_counts.items()):
        error_rate = count / total_counts[ref_3mer]
        print(f'{ref_3mer}\t{read_3mer}\t{count}\t{error_rate}')


def align_sequences(read_seq, ref_seq, expanded_cigar):
    aligned_read_seq = []
    aligned_ref_seq = []
    read_i, ref_i = 0, 0

    for op in expanded_cigar:
        if op in {'M', '=', 'X'}:
            aligned_read_seq.append(read_seq[read_i])
            aligned_ref_seq.append(ref_seq[ref_i])
            read_i += 1
            ref_i += 1
        elif op == 'I':
            aligned_read_seq.append(read_seq[read_i])
            aligned_ref_seq.append('-')
            read_i += 1
        else:  # D
            aligned_read_seq.append('-')
            aligned_ref_seq.append(ref_seq[ref_i])
            ref_i += 1

    assert read_i == len(read_seq)
    assert ref_i == len(ref_seq)
    return ''.join(aligned_read_seq), ''.join(aligned_ref_seq)


def canonicalise_3mer(ref_3mer):
    rev_ref_3mer = reverse_complement(ref_3mer)
    return min(ref_3mer, rev_ref_3mer)


def canonicalise_3mer_substitution(ref_3mer, read_3mer):
    canonical_ref_3mer = canonicalise_3mer(ref_3mer)
    if canonical_ref_3mer != ref_3mer:
        return canonical_ref_3mer, reverse_complement(read_3mer)
    return ref_3mer, read_3mer


def load_alignments(filename):
    alignments = []
    with get_open_func(filename)(filename, 'rt') as paf_file:
        for line in paf_file:
            alignments.append(Alignment(line))
    return alignments


def get_expanded_cigar(cigar):
    return ''.join(op * int(size) for size, op in re.findall(r'(\d+)([MID=X])', cigar))


def keep_only_best_alignment_per_read(alignments):
    alignments_per_read = collections.defaultdict(list)
    for a in alignments:
        alignments_per_read[a.query_name].append(a)

    return [max(read_alignments, key=lambda a: a.matching_bases)
            for read_alignments in alignments_per_read.values()]


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type in {'bz2', 'zip'}:
        raise ValueError('Use plain text or gzip-compressed files')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open

def iterate_fastq(filename):
    with get_open_func(filename)(filename, 'rt') as fastq:
        for line in fastq:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith('@'):
                continue
            name = line[1:].split()[0]
            header = line
            sequence = next(fastq).strip().upper()
            _ = next(fastq)
            qualities = next(fastq).strip()
            yield name, header, sequence, qualities


def iterate_fasta(filename):
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    yield name.split()[0], ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            yield name.split()[0], ''.join(sequence)


def load_fastq_as_dict(filename):
    return {name: seq for name, _, seq, _ in iterate_fastq(filename)}


def load_fasta_as_dict(filename):
    return {name: seq for name, seq in iterate_fasta(filename)}


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        assert len(line_parts) >= 11

        self.query_name = line_parts[0]
        self.query_length = int(line_parts[1])
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.target_name = line_parts[5]
        self.target_length = int(line_parts[6])
        self.target_start = int(line_parts[7])
        self.target_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases

        self.query_cov = 100.0 * (self.query_end - self.query_start) / self.query_length

        self.cigar, self.alignment_score = None, None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])

    def __repr__(self):
        return self.query_name + ':' + str(self.query_start) + '-' + str(self.query_end) + \
               '(' + self.strand + '), ' + \
               self.target_name + ':' + str(self.target_start) + '-' + str(self.target_end) + \
               ' (' + ('%.3f' % self.percent_identity) + '%)'


if __name__ == '__main__':
    main()
