#!/usr/bin/env python3

import csv
import gzip
import math
import os
import re
import sys

from Bio import SeqIO


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} FASTQ PAF")

    fastq_path, paf_path = sys.argv[1:]
    best_alignments = load_best_alignments(paf_path)
    read_stats = collect_read_stats(fastq_path, best_alignments)
    write_read_stats(read_stats, sys.stdout)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def get_tool_name(path):
    basename = os.path.basename(path)
    if basename.endswith(".gz"):
        basename = basename[:-3]
    return os.path.splitext(basename)[0]


def get_paf_tag(fields, prefix, default=None):
    for field in fields:
        if field.startswith(prefix):
            return field[len(prefix):]
    return default


def parse_paf_alignment(fields):
    return {
        "read_name": fields[0],
        "read_length": int(fields[1]),
        "query_start": int(fields[2]),
        "query_end": int(fields[3]),
        "matching_bases": int(fields[9]),
        "alignment_block_length": int(fields[10]),
        "tags": fields[12:],
    }


def load_best_alignments(path):
    """Keep the alignment with most matching bases per read; first wins ties."""
    best_alignments = {}
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            alignment = parse_paf_alignment(fields)
            read_name = alignment["read_name"]
            previous = best_alignments.get(read_name)
            if previous is None or alignment["matching_bases"] > previous["matching_bases"]:
                best_alignments[read_name] = alignment

    return best_alignments


def count_cigar_operations(cigar):
    counts = {"I": 0, "D": 0, "=": 0, "X": 0}
    for length, operation in re.findall(r"(\d+)([ID=X])", cigar):
        counts[operation] += int(length)
    return counts


def count_alignment_outcomes(alignment, cigar):
    counts = count_cigar_operations(cigar)
    if "=" in cigar or "X" in cigar:
        matches, substitutions = counts["="], counts["X"]
    else:
        edit_distance = int(get_paf_tag(alignment["tags"], "NM:i:", "0"))
        matches = alignment["matching_bases"]
        substitutions = max(edit_distance - counts["I"] - counts["D"], 0)
    return matches, substitutions, counts["I"], counts["D"]


def calculate_empirical_stats(alignment):
    stats = dict.fromkeys([
        "empirical_accuracy", "empirical_sub_rate", "empirical_ins_rate", "empirical_del_rate"
    ])
    if alignment is None:
        return stats

    cigar = get_paf_tag(alignment["tags"], "cg:Z:")
    if cigar:
        matches, substitutions, insertions, deletions = count_alignment_outcomes(alignment, cigar)
        total = matches + substitutions + insertions + deletions
        if total > 0:
            stats["empirical_accuracy"] = matches / total
            stats["empirical_sub_rate"] = substitutions / total
            stats["empirical_ins_rate"] = insertions / total
            stats["empirical_del_rate"] = deletions / total
    elif alignment["alignment_block_length"] > 0:
        stats["empirical_accuracy"] = (
            alignment["matching_bases"] / alignment["alignment_block_length"]
        )
    return stats


def reported_stats_from_phred(phred_scores):
    if not phred_scores:
        return None, None
    mean_error = sum(10 ** (-q / 10) for q in phred_scores) / len(phred_scores)
    return 1 - mean_error, -10 * math.log10(mean_error)


def empirical_qscore_from_accuracy(empirical_accuracy):
    if empirical_accuracy is None:
        return None
    if empirical_accuracy >= 1.0:
        return math.inf
    return -10 * math.log10(1.0 - empirical_accuracy)


def calculate_alignment_coverage(read_length, alignment):
    if alignment is None:
        return {
            "unaligned_start": None,
            "unaligned_end": None,
            "aligned_length": 0,
            "aligned_fraction": 0.0,
            "aligned": False,
        }
    unaligned_start = alignment["query_start"]
    unaligned_end = alignment["read_length"] - alignment["query_end"]
    aligned_length = read_length - unaligned_start - unaligned_end
    return {
        "unaligned_start": unaligned_start,
        "unaligned_end": unaligned_end,
        "aligned_length": aligned_length,
        "aligned_fraction": aligned_length / read_length if read_length else 0.0,
        "aligned": True,
    }


def calculate_reported_stats(phred_scores, coverage):
    whole_accuracy, whole_qscore = reported_stats_from_phred(phred_scores)
    aligned_scores = []
    if coverage["aligned"]:
        # PAF query coordinates refer to the original read on either strand.
        start = coverage["unaligned_start"]
        end = len(phred_scores) - coverage["unaligned_end"]
        aligned_scores = phred_scores[start:end]
    aligned_accuracy, aligned_qscore = reported_stats_from_phred(aligned_scores)
    return {
        "reported_accuracy_whole_read": whole_accuracy,
        "reported_qscore_whole_read": whole_qscore,
        "reported_accuracy_aligned_region": aligned_accuracy,
        "reported_qscore_aligned_region": aligned_qscore,
    }


def calculate_read_stats(record, alignment):
    sequence = str(record.seq).upper()
    coverage = calculate_alignment_coverage(len(sequence), alignment)
    empirical = calculate_empirical_stats(alignment)
    reported = calculate_reported_stats(record.letter_annotations["phred_quality"], coverage)
    return {
        "read_name": record.id,
        "read_length": len(sequence),
        **coverage,
        **empirical,
        "empirical_qscore": empirical_qscore_from_accuracy(empirical["empirical_accuracy"]),
        **reported,
        "gc_content": (sequence.count("G") + sequence.count("C")) / len(sequence)
        if sequence else None,
    }


def collect_read_stats(fastq_path, best_alignments):
    tool_name = get_tool_name(fastq_path)
    read_stats = {}
    with open_text(fastq_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            stats = calculate_read_stats(record, best_alignments.get(record.id))
            stats["tool"] = tool_name
            read_stats[record.id] = stats
    return read_stats.values()


def write_read_stats(read_stats, output):
    columns = [
        "tool",
        "read_name",
        "read_length",
        "aligned",
        "unaligned_start",
        "unaligned_end",
        "aligned_length",
        "aligned_fraction",
        "empirical_accuracy",
        "empirical_qscore",
        "empirical_sub_rate",
        "empirical_ins_rate",
        "empirical_del_rate",
        "reported_accuracy_whole_read",
        "reported_qscore_whole_read",
        "reported_accuracy_aligned_region",
        "reported_qscore_aligned_region",
        "gc_content",
    ]
    writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(read_stats)


if __name__ == '__main__':
    main()
