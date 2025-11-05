#!/usr/bin/env python3

import pandas as pd
import re
from Bio import SeqIO
import numpy as np
import gzip
from pathlib import Path

paf_path = "/home/taouk/ONT_read_sim/mapped_sim_reads/Badread_primary.paf"

df = pd.read_csv(paf_path, sep='\t', header=None)

#clean data formats
df[[1, 2, 3, 6, 7, 8, 9, 10, 11]] = df[[1, 2, 3, 6, 7, 8, 9, 10, 11]].astype(int)


#append data we care about to end of the orginal dataframe
df["tool"] = "Badread"
df["read_name"] = df[0]
df["read_length"] = df[1]
df["unaligned_bases_at_start"] = df[2]
df["unaligned_bases_at_end"] = (df[1] - df[3])
df["actual_identity"] = (df[9] / df[10]).astype(float)


# helper functions
def parse_cigar(cigar):
    ops = re.findall(r'(\d+)([MID])', cigar)
    counts = {'M':0,'I':0,'D':0}
    for length, op in ops:
        counts[op] += int(length)
    return counts

def get_nm(fields):
    for f in fields:
        if f.startswith("NM:i:"):
            return int(f.split(":")[-1])
    return 0

# append substitution, insertion, deletion rates
sub_rate_list = []
ins_rate_list = []
del_rate_list = []

for idx, row in df.iterrows():
    cigar_field = [str(x) for x in row[12:] if str(x).startswith("cg:Z:")]
    if cigar_field:
        cigar_counts = parse_cigar(cigar_field[0].split(":")[-1])
        aligned_len = cigar_counts['M'] + cigar_counts['I']
        nm = get_nm(row[12:])
        subs = max(nm - (cigar_counts['I'] + cigar_counts['D']), 0)
        sub_rate_list.append(subs / aligned_len)
        ins_rate_list.append(cigar_counts['I'] / aligned_len)
        del_rate_list.append(cigar_counts['D'] / aligned_len)
    else:
        sub_rate_list.append(None)
        ins_rate_list.append(None)
        del_rate_list.append(None)

df["sub_rate"] = sub_rate_list
df["ins_rate"] = ins_rate_list
df["del_rate"] = del_rate_list


fastq_file = "/home/taouk/ONT_read_sim/Badread/simulated_reads.fastq.gz"

fastq_metrics = {}
with gzip.open(fastq_file, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        seq = str(record.seq).upper()
        phred_scores = record.letter_annotations["phred_quality"]
        mean_q = np.mean(phred_scores)
        gc = (seq.count("G") + seq.count("C")) / len(seq) if len(seq) > 0 else 0
        fastq_metrics[record.id] = {"mean_qscore": mean_q, "gc_content": gc}

# Add FASTQ metrics to dataframe
df["mean_qscore"] = df["read_name"].map(lambda x: fastq_metrics.get(x, {}).get("mean_qscore", np.nan))
df["gc_content"] = df["read_name"].map(lambda x: fastq_metrics.get(x, {}).get("gc_content", np.nan))

# Mark aligned reads
df["aligned"] = True  # all reads present in PAF are aligned

# include unaligned reads from FASTQ
all_reads = set(fastq_metrics.keys())
aligned_reads = set(df["read_name"])
unaligned_reads = all_reads - aligned_reads

if unaligned_reads:
    unaligned_df = pd.DataFrame({
        "read_name": list(unaligned_reads),
        "aligned": False,
        "mean_qscore": [fastq_metrics[r]["mean_qscore"] for r in unaligned_reads],
        "gc_content": [fastq_metrics[r]["gc_content"] for r in unaligned_reads]
    })
    # Keep same columns as df; fill others with NaN
    for col in df.columns:
        if col not in unaligned_df.columns:
            unaligned_df[col] = np.nan
    df = pd.concat([df, unaligned_df], ignore_index=True, sort=False)


# output
output_cols = ["tool", "read_name","read_length", "unaligned_bases_at_start", "unaligned_bases_at_end", "actual_identity", "sub_rate", "ins_rate", "del_rate", "aligned", "mean_qscore", "gc_content", "sub_rate"]
df[output_cols].to_csv("Badread_stats.tsv", sep='\t', index=False)