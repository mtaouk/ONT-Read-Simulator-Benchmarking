#!/usr/bin/env python3

import gzip
import csv
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


def parse_cigar(cigar):
    ops = re.findall(r'(\d+)([MID])', cigar)
    counts = {"M": 0, "I": 0, "D": 0}
    for length, op in ops:
        counts[op] += int(length)
    return counts


def get_nm(fields):
    for field in fields:
        if isinstance(field, str) and field.startswith("NM:i:"):
            return int(field.split(":")[-1])
    return 0


# --- Load FASTQ and calculate read-level stats ---
fastq_dict = {}
handle = gzip.open(fastq_path, "rt") if fastq_path.endswith(".gz") else open(fastq_path, "r")
for record in SeqIO.parse(handle, "fastq"):
    seq = str(record.seq).upper()
    mean_q = sum(record.letter_annotations["phred_quality"]) / len(record)
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    fastq_dict[record.id] = {"read_length": len(seq), "mean_qscore": mean_q, "gc_content": gc}
handle.close()

fastq_df = (
    pd.DataFrame.from_dict(fastq_dict, orient="index")
    .reset_index()
    .rename(columns={"index": "read_name"})
)

# --- Load PAF ---
df_paf = pd.read_csv(paf_path, sep="\t", header=None, engine="python", on_bad_lines="warn", dtype=str)

# convert numeric columns
for col in [1, 2, 3, 6, 7, 8, 9, 10, 11]:
    df_paf[col] = pd.to_numeric(df_paf[col], errors="coerce")

df_paf["read_name"] = df_paf[0].astype(str)
df_paf["unaligned_start"] = df_paf[2]
df_paf["unaligned_end"] = df_paf[1] - df_paf[3]
df_paf["actual_identity"] = (df_paf[9] / df_paf[10]).astype(float)

sub_rate_list, ins_rate_list, del_rate_list = [], [], []

for _, row in df_paf.iterrows():
    cigar_field = [x for x in row[12:] if isinstance(x, str) and x.startswith("cg:Z:")]
    if cigar_field:
        cigar_counts = parse_cigar(cigar_field[0].split(":")[-1])
        aligned_len = cigar_counts["M"] + cigar_counts["I"]
        nm = get_nm(row[12:])
        subs = max(nm - (cigar_counts["I"] + cigar_counts["D"]), 0)
        sub_rate_list.append(subs / aligned_len)
        ins_rate_list.append(cigar_counts["I"] / aligned_len)
        del_rate_list.append(cigar_counts["D"] / aligned_len)
    else:
        sub_rate_list.append(None)
        ins_rate_list.append(None)
        del_rate_list.append(None)

df_paf["sub_rate"] = sub_rate_list
df_paf["ins_rate"] = ins_rate_list
df_paf["del_rate"] = del_rate_list

df_paf = df_paf[["read_name", "unaligned_start", "unaligned_end", "actual_identity", "sub_rate", "ins_rate", "del_rate"]]

# --- Merge FASTQ and PAF ---
merged_df = pd.merge(fastq_df, df_paf, on="read_name", how="left")
merged_df["aligned"] = merged_df["actual_identity"].notna()

# --- Columns to write ---
output_cols = ["read_name", "read_length", "unaligned_start", "unaligned_end", "actual_identity", "sub_rate", "ins_rate", "del_rate", "aligned", "mean_qscore", "gc_content"]

merged_df[output_cols].to_csv(sys.stdout, sep="\t", index=False)
