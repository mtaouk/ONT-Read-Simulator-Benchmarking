#!/usr/bin/env python3
"""Print assembly statistics, whole/trimmed identities and completeness as one TSV row."""

import sys
from collections import defaultdict
from itertools import accumulate
from pathlib import Path


def fasta_lengths(filename):
    records = {}
    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                records[header] = 0
            elif line.strip():
                records[header] += len(line.strip())
    return records


def alignment_stats(paf, ref_lengths):
    if not Path(paf).exists():
        return 0, "NA", []
    with open(paf) as f:
        alignments = [line.split() for line in f]
    # Longest reference span first, then highest identity; ties retain PAF order.
    alignments.sort(key=lambda a: (int(a[8]) - int(a[7]), int(a[9]) / int(a[10])), reverse=True)
    kept_intervals, reference_intervals = defaultdict(list), defaultdict(list)
    kept_alignments = []
    matches = bases = 0
    for a in alignments:
        target = a[5]
        start, end, length = int(a[7]), int(a[8]), ref_lengths[target]
        span = end - start
        start %= length
        # Wrap each alignment into one or two non-overlapping reference intervals.
        if span >= length:
            intervals = [(0, length)]
        else:
            end = start + span
            intervals = [(start, min(end, length))]
            if end > length:
                intervals.append((0, end - length))
        # Sum overlap across the origin, comparing with each kept alignment separately.
        if any(sum(max(0, min(right, stop) - max(left, begin))
                   for left, right in intervals for begin, stop in kept) > 10
               for kept in kept_intervals[target]):
            continue
        kept_intervals[target].append(intervals)
        reference_intervals[target].extend(intervals)
        kept_alignments.append(a)
        matches += int(a[9])
        bases += int(a[10])

    covered = 0
    for intervals in reference_intervals.values():
        previous_end = 0
        for start, end in sorted(intervals):
            covered += max(0, end - max(start, previous_end))
            previous_end = max(previous_end, end)
    coverage = covered / sum(ref_lengths.values())
    identity = f"{matches / bases:.12f}" if bases else "NA"
    return coverage, identity, kept_alignments


def main():
    tool, depth, assembler, assembly, paf_whole, paf_trimmed, reference = sys.argv[1:]
    if not Path(assembly).exists():
        print(tool, depth, assembler, *(["NA"] * 6), "false", sep="\t")
        return

    contigs = fasta_lengths(assembly)
    lengths = sorted(contigs.values(), reverse=True)
    size = sum(lengths)
    n50 = next((n for n, total in zip(lengths, accumulate(lengths)) if total >= size / 2), 0)
    complete = (len(lengths) == 2
                and all(abs(n - expected) <= expected * 0.01
                        for n, expected in zip(lengths, (2883001, 89544))))
    ref_lengths = {header.split()[0]: n for header, n in fasta_lengths(reference).items()}
    coverage, identity_whole, _ = alignment_stats(paf_whole, ref_lengths)
    identity_trimmed = "NA"
    if lengths and lengths[0] >= 2000000:
        _, identity, kept = alignment_stats(paf_trimmed, ref_lengths)
        if len(kept) == 1 and int(kept[0][2]) == 0 and int(kept[0][3]) == int(kept[0][1]):
            identity_trimmed = identity
    print(tool, depth, assembler, size, len(lengths), n50, f"{coverage:.12f}",
          identity_whole, identity_trimmed,
          str(complete).lower(), sep="\t")


if __name__ == "__main__":
    main()
