#!/usr/bin/env python3
"""
This script counts the instances of each qscore, with counts divided into this groups:
- match: bases where the read matches the reference and is not adjacent to a deletion
- mismatch: bases where the read and reference differ
- insertion: bases in the read which are not in the reference
- deletion: bases where the read matches the reference and is adjacent to a deletion

Deletions are a bit tricky, since a deleted base isn't in the read and so doesn't have any qscore
at all. This script therefore considers positions adjacent to a deletion (on either side) to be
deletion positions.

Usage:
  per_base_qscores.py reads.fastq.gz reference.fasta alignments.paf > qscore_counts.tsv
"""

import collections
import gzip
import re
import sys

# Turning on verbose mode displays the aligned sequences to stdout. For debugging purposes.
VERBOSE = False


def main():
    read_filename = sys.argv[1]
    ref_filename = sys.argv[2]
    paf_filename = sys.argv[3]

    reads = load_reads(read_filename)
    ref = load_reference(ref_filename)
    alignments = load_alignments(paf_filename)
    alignments = keep_only_best_alignment_per_read(alignments)

    # Distributions max out at Q93, since that's the highest possible FASTQ qscore.
    match_qscore_distribution = [0] * 94
    mismatch_qscore_distribution = [0] * 94
    insertion_qscore_distribution = [0] * 94
    deletion_adjacent_qscore_distribution = [0] * 94

    for a in alignments:
        if VERBOSE:
            print()
            print(a)
            print(a.cigar)

        # Trim the sequences to the aligned regions.
        ref_seq = ref[a.ref_name][a.ref_start:a.ref_end]
        read_seq = reads[a.read_name][0][a.read_start:a.read_end]
        read_qual = reads[a.read_name][1][a.read_start:a.read_end]

        # Flip the strand, if necessary.
        if a.strand == '-':
            read_seq = reverse_complement(read_seq)
            read_qual = read_qual[::-1]

        # Add gaps ('-' in sequences, spaces in qscores) so the sequences are aligned.
        aligned_ref_seq, aligned_read_seq, aligned_read_qual = align(ref_seq, read_seq, read_qual, a.cigar)

        if VERBOSE:
            print(aligned_ref_seq)
            print(aligned_read_seq)
            print(aligned_read_qual)
            print()

        # Count the qscores for each read base.
        for i, read_base in enumerate(aligned_read_seq):
            if read_base == '-':
                continue
            read_qual = qscore_letter_to_num(aligned_read_qual[i])
            ref_base = aligned_ref_seq[i]
            if ref_base == '-':
                insertion_qscore_distribution[read_qual] += 1
            elif read_base != ref_base:
                mismatch_qscore_distribution[read_qual] += 1
            else:  # match
                if adjacent_to_deletion(i, aligned_read_seq):
                    deletion_adjacent_qscore_distribution[read_qual] += 1
                else:
                    match_qscore_distribution[read_qual] += 1

    # Print the counts to stdout.
    print('qscore\tmatch\tmismatch\tinsertion\tdeletion')
    for i in range(94):
        print(f'{i}\t{match_qscore_distribution[i]}'
              f'\t{mismatch_qscore_distribution[i]}'
              f'\t{insertion_qscore_distribution[i]}'
              f'\t{deletion_adjacent_qscore_distribution[i]}')


def align(ref_seq, read_seq, read_qual, cigar):
    aligned_ref = []
    aligned_read = []
    aligned_qual = []

    r_i, q_i = 0, 0

    if re.search(r'[^\dMID=X]', cigar):
        raise ValueError(f"Unexpected CIGAR op in: {cigar}")

    for length_str, op in re.findall(r'(\d+)([MID=X])', cigar):
        n = int(length_str)

        if op in ('M', '=', 'X'):
            aligned_ref.append(ref_seq[r_i:r_i + n])
            aligned_read.append(read_seq[q_i:q_i + n])
            aligned_qual.append(read_qual[q_i:q_i + n])
            r_i += n
            q_i += n

        elif op == 'I':
            aligned_ref.append('-' * n)
            aligned_read.append(read_seq[q_i:q_i + n])
            aligned_qual.append(read_qual[q_i:q_i + n])
            q_i += n

        elif op == 'D':
            aligned_ref.append(ref_seq[r_i:r_i + n])
            aligned_read.append('-' * n)
            aligned_qual.append(' ' * n)
            r_i += n

        else:
            assert False

    assert r_i == len(ref_seq) and q_i == len(read_seq)
    aligned_ref = ''.join(aligned_ref)
    aligned_read = ''.join(aligned_read)
    aligned_qual = ''.join(aligned_qual)
    assert len(aligned_ref) == len(aligned_read) == len(aligned_qual)
    return aligned_ref, aligned_read, aligned_qual


def load_reads(filename):
    """
    Loads a FASTQ file. Returns a {name: (seq, qual)} dictionary.
    """
    reads = {}
    with get_open_func(filename)(filename, 'rt') as fastq:
        for line in fastq:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith('@'):
                continue
            name = line[1:].split()[0]
            seq = next(fastq).strip().upper()
            assert_only_canonical_bases(seq)
            _ = next(fastq)
            qual = next(fastq).strip()
            assert len(seq) == len(qual)
            reads[name] = (seq, qual)
    return reads


def load_reference(filename):
    """
    Loads a FASTA file. Returns a {name: seq} dictionary.
    """
    ref = {}
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    name_parts = name.split(maxsplit=1)
                    contig_name = name_parts[0]
                    ref[contig_name] = ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            name_parts = name.split(maxsplit=1)
            contig_name = name_parts[0]
            ref[contig_name] = ''.join(sequence)
    for seq in ref.values():
        assert_only_canonical_bases(seq)
    return ref


def assert_only_canonical_bases(seq: str) -> None:
    bad = set(seq) - {"A", "C", "G", "T"}
    assert not bad, f"Non-ACGT base(s) found: {''.join(sorted(bad))}"


def load_alignments(filename):
    """
    Loads a PAF file, returns a list of Alignment objects.
    """
    alignments = []
    with get_open_func(filename)(filename, 'rt') as paf:
        for line in paf:
            alignments.append(Alignment(line))
    return alignments


def keep_only_best_alignment_per_read(alignments):
    """
    Returns a culled list of Alignment objects, with only the best alignment (as judged by
    matching bases, PAF column 10) for each read.
    """
    alignments_per_read = collections.defaultdict(list)
    for a in alignments:
        alignments_per_read[a.read_name].append(a)
    return [max(read_alignments, key=lambda a: a.matching_bases)
            for read_alignments in alignments_per_read.values()]


def qscore_letter_to_num(q):
    num = ord(q) - 33
    assert 0 <= num < 94  # the highest normal character is '~' = Q93
    return num


def adjacent_to_deletion(i, aligned_read_seq):
    if i > 0 and aligned_read_seq[i-1] == '-':
        return True
    if i < len(aligned_read_seq)-1 and aligned_read_seq[i+1] == '-':
        return True
    return False


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
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('\nError: alignment file does not seem to be in PAF format')

        self.read_name = line_parts[0]
        self.read_start = int(line_parts[2])
        self.read_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases

        self.cigar = None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]

    def __repr__(self):
        return self.read_name + ':' + str(self.read_start) + '-' + str(self.read_end) + \
               '(' + self.strand + '), ' + \
               self.ref_name + ':' + str(self.ref_start) + '-' + str(self.ref_end) + \
               ' (' + ('%.3f' % self.percent_identity) + '%)'


if __name__ == '__main__':
    main()
