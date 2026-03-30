#!/usr/bin/env python3

import gzip
import csv
import math
import os
import re
import sys

import pandas as pd
from Bio import SeqIO

csv.field_size_limit(10**7)  # 10 million chars

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} FASTQ PAF", file=sys.stderr)
    sys.exit(1)

fastq_path = sys.argv[1]
paf_path = sys.argv[2]


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def get_tool_name(path):
    basename = os.path.basename(path)
    if basename.endswith(".gz"):
        basename = basename[:-3]
    for ext in (".fastq", ".fq"):
        if basename.lower().endswith(ext):
            return basename[: -len(ext)]
    return os.path.splitext(basename)[0]


def parse_cigar(cigar):
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    counts = {"M": 0, "I": 0, "D": 0, "=": 0, "X": 0}
    for length, op in ops:
        if op in counts:
            counts[op] += int(length)
    return counts


def get_nm(fields):
    for field in fields:
        if isinstance(field, str) and field.startswith("NM:i:"):
            return int(field.split(":")[-1])
    return 0


def reported_stats_from_phred(phred_scores):
    mean_error = sum(10 ** (-q / 10) for q in phred_scores) / len(phred_scores)
    reported_accuracy = 1 - mean_error
    reported_qscore = -10 * math.log10(mean_error)
    return reported_accuracy, reported_qscore


def empirical_qscore_from_accuracy(empirical_accuracy):
    if pd.isna(empirical_accuracy):
        return None
    if empirical_accuracy >= 1.0:
        return math.inf
    return -10 * math.log10(1.0 - empirical_accuracy)


def load_paf(path):
    best_rows = {}
    skipped_short = 0
    skipped_invalid = 0

    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 12:
                skipped_short += 1
                continue

            try:
                read_length = int(fields[1])
                unaligned_start = int(fields[2])
                read_end = int(fields[3])
                matches = int(fields[9])
                alignment_block_length = int(fields[10])
            except ValueError:
                skipped_invalid += 1
                continue

            cigar = None
            for field in fields[12:]:
                if field.startswith("cg:Z:"):
                    cigar = field[5:]
                    break

            empirical_accuracy = None
            empirical_sub_rate = None
            empirical_ins_rate = None
            empirical_del_rate = None
            if cigar:
                cigar_counts = parse_cigar(cigar)
                has_eqx = "=" in cigar or "X" in cigar
                if has_eqx:
                    total_len = (
                        cigar_counts["="]
                        + cigar_counts["X"]
                        + cigar_counts["I"]
                        + cigar_counts["D"]
                    )
                    if total_len > 0:
                        empirical_accuracy = cigar_counts["="] / total_len
                        empirical_sub_rate = cigar_counts["X"] / total_len
                        empirical_ins_rate = cigar_counts["I"] / total_len
                        empirical_del_rate = cigar_counts["D"] / total_len
                else:
                    nm = get_nm(fields[12:])
                    subs = max(nm - (cigar_counts["I"] + cigar_counts["D"]), 0)
                    total_len = matches + subs + cigar_counts["I"] + cigar_counts["D"]
                    if total_len > 0:
                        empirical_accuracy = matches / total_len
                        empirical_sub_rate = subs / total_len
                        empirical_ins_rate = cigar_counts["I"] / total_len
                        empirical_del_rate = cigar_counts["D"] / total_len
            elif alignment_block_length > 0:
                empirical_accuracy = matches / alignment_block_length

            read_name = fields[0]
            best_row = best_rows.get(read_name)
            if best_row is None or matches > best_row["matching_bases"]:
                best_rows[read_name] = {
                    "read_name": read_name,
                    "unaligned_start": unaligned_start,
                    "unaligned_end": read_length - read_end,
                    "empirical_accuracy": empirical_accuracy,
                    "empirical_sub_rate": empirical_sub_rate,
                    "empirical_ins_rate": empirical_ins_rate,
                    "empirical_del_rate": empirical_del_rate,
                    "aligned": True,
                    "matching_bases": matches,
                }

    if skipped_short or skipped_invalid:
        print(
            f"Skipped {skipped_short + skipped_invalid} malformed PAF lines "
            f"({skipped_short} short, {skipped_invalid} with invalid numeric fields).",
            file=sys.stderr,
        )

    return pd.DataFrame(
        best_rows.values(),
        columns=[
            "read_name",
            "unaligned_start",
            "unaligned_end",
            "empirical_accuracy",
            "empirical_sub_rate",
            "empirical_ins_rate",
            "empirical_del_rate",
            "aligned",
        ],
    )


# --- Load FASTQ and calculate read-level stats ---
fastq_dict = {}
with open_text(fastq_path) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        seq = str(record.seq).upper()
        reported_accuracy, reported_qscore = reported_stats_from_phred(
            record.letter_annotations["phred_quality"]
        )
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        fastq_dict[record.id] = {
            "read_length": len(seq),
            "reported_accuracy": reported_accuracy,
            "reported_qscore": reported_qscore,
            "gc_content": gc,
        }

fastq_df = (
    pd.DataFrame.from_dict(fastq_dict, orient="index")
    .reset_index()
    .rename(columns={"index": "read_name"})
)
fastq_df["tool"] = get_tool_name(fastq_path)

# --- Load PAF ---
df_paf = load_paf(paf_path)

# --- Merge FASTQ and PAF ---
merged_df = pd.merge(fastq_df, df_paf, on="read_name", how="left")
merged_df["aligned"] = merged_df["aligned"].astype("boolean").fillna(False)
merged_df["unaligned_start"] = merged_df["unaligned_start"].astype("Int64")
merged_df["unaligned_end"] = merged_df["unaligned_end"].astype("Int64")
merged_df["empirical_qscore"] = merged_df["empirical_accuracy"].apply(
    empirical_qscore_from_accuracy
)

# --- Columns to write ---
output_cols = [
    "tool",
    "read_name",
    "read_length",
    "unaligned_start",
    "unaligned_end",
    "empirical_accuracy",
    "empirical_qscore",
    "empirical_sub_rate",
    "empirical_ins_rate",
    "empirical_del_rate",
    "aligned",
    "reported_accuracy",
    "reported_qscore",
    "gc_content",
]

merged_df[output_cols].to_csv(sys.stdout, sep="\t", index=False)
