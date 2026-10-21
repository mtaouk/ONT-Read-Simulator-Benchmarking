#!/usr/bin/env python3
"""
This script reports the substitution:insertion:deletion ratio from a PAF file. The PAF must have
been generated with --eqx CIGAR strings. Only the alignment with the most matching bases
per read is used; the first alignment wins ties.

Usage:
  sub_ins_del.py alignments.paf
"""

import re
import sys


def main():
    paf_filename = sys.argv[1]
    substitutions, insertions, deletions = count_error_types(paf_filename)
    total = substitutions + insertions + deletions
    substitutions /= total
    insertions /= total
    deletions /= total
    print(f'{paf_filename}\t{substitutions*100:.1f}%\t{insertions*100:.1f}%\t{deletions*100:.1f}%')


def load_best_alignments(path):
    best_alignments = {}
    with open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            alignment = {
                "read_name": fields[0],
                "matching_bases": int(fields[9]),
                "tags": fields[12:],
            }
            read_name = alignment["read_name"]
            previous = best_alignments.get(read_name)
            if previous is None or alignment["matching_bases"] > previous["matching_bases"]:
                best_alignments[read_name] = alignment
    return best_alignments


def count_error_types(paf_filename):
    counts = {"X": 0, "I": 0, "D": 0}
    for alignment in load_best_alignments(paf_filename).values():
        for part in alignment["tags"]:
            if part.startswith("cg:Z:"):
                for p in re.findall(r"\d+[IDX=M]", part[5:]):
                    op = p[-1]
                    if op in counts:
                        counts[op] += int(p[:-1])
                break
    return counts["X"], counts["I"], counts["D"]


if __name__ == '__main__':
    main()
