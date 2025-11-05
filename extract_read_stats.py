#!/usr/bin/env python3

import pandas as pd
import gzip
from Bio import SeqIO
import re
from pathlib import Path
import csv

# tool name -> (paf, fastq)
files = {
    "Badread": ("/home/taouk/ONT_read_sim/mapped_sim_reads/Badread_primary.paf", "/home/taouk/ONT_read_sim/Badread/simulated_reads.fastq.gz"),
    "LongISLND": ("/home/taouk/ONT_read_sim/mapped_sim_reads/LongISLND_primary.paf", "/home/taouk/ONT_read_sim/LongISLND/longislnd_reads.fastq.gz"),
    "lrsim": ("/home/taouk/ONT_read_sim/mapped_sim_reads/lrsim_primary.paf", "/home/taouk/ONT_read_sim/lrsim/lrsim_reads.fq.gz"),
    "NanoSim": ("/home/taouk/ONT_read_sim/mapped_sim_reads/NanoSim_primary.paf", "/home/taouk/ONT_read_sim/NanoSim/simulation/simulated_aligned_reads.fastq.gz"),
    "PBSIM3": ("/home/taouk/ONT_read_sim/mapped_sim_reads/PBSIM3_primary.paf", "/home/taouk/ONT_read_sim/PBSIM3/pbsm_reads.fastq.gz"),
    "simlord": ("/home/taouk/ONT_read_sim/mapped_sim_reads/simlord_primary.paf", "/home/taouk/ONT_read_sim/simlord/simlord_ont_like.fastq.gz"),
    "real_reads": ("/home/taouk/ONT_read_sim/mapped_sim_reads/real_reads_primary.paf", "/home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_shortnames.fastq.gz")
}

# empty df for stacking
all_dfs = []

csv.field_size_limit(10**7)  # 10 million chars, adjust if needed

for tool, (paf_path, fastq_path) in files.items():
    print(f"Processing {tool}...")

    df = pd.read_csv(
    paf_path,
    sep='\t',
    header=None,
    engine="python",
    on_bad_lines="warn"
)
    df[[1, 2, 3, 6, 7, 8, 9, 10, 11]] = df[[1, 2, 3, 6, 7, 8, 9, 10, 11]].astype(int)

    # add your existing derived fields
    df["read_name"] = df[0]
    df["read_length"] = df[1]
    df["unaligned_start"] = df[2]
    df["unaligned_end"] = df[1] - df[3]
    df["actual_identity"] = (df[9] / df[10]).astype(float)

    # --- parse cigar and rates (reuse your helper functions) ---
    def parse_cigar(cigar):
        ops = re.findall(r'(\d+)([MID])', cigar)
        counts = {'M':0,'I':0,'D':0}
        for length, op in ops:
            counts[op] += int(length)
        return counts

    def get_nm(fields):
        for f in fields:
            if isinstance(f, str) and f.startswith("NM:i:"):
                return int(f.split(":")[-1])
        return 0

    sub_rate_list, ins_rate_list, del_rate_list = [], [], []

    for _, row in df.iterrows():
        cigar_field = [x for x in row[12:] if isinstance(x, str) and x.startswith("cg:Z:")]
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

    # --- if there's a fastq, extract mean qscore + GC ---
    if fastq_path:
        qscore_dict = {}
        gc_dict = {}
        handle = gzip.open(fastq_path, "rt") if fastq_path.endswith(".gz") else open(fastq_path, "r")

        for record in SeqIO.parse(handle, "fastq"):
            mean_q = sum(record.letter_annotations["phred_quality"]) / len(record)
            seq = str(record.seq).upper()
            gc = (seq.count("G") + seq.count("C")) / len(seq)
            qscore_dict[record.id] = mean_q
            gc_dict[record.id] = gc
        handle.close()

        df["mean_qscore"] = df["read_name"].map(qscore_dict)
        df["gc_content"] = df["read_name"].map(gc_dict)
    else:
        df["mean_qscore"] = None
        df["gc_content"] = None

    df["tool"] = tool

    df["aligned"] = True

    # collect each df
    all_dfs.append(df)

# --- merge like rbind ---
combined = pd.concat(all_dfs, ignore_index=True)


# define the exact columns in the final output
output_cols = ["tool", "read_name", "read_length", "unaligned_start", "unaligned_end", "actual_identity", "sub_rate", "ins_rate", "del_rate", "aligned", "mean_qscore", "gc_content"]

# subset and write
combined[output_cols].to_csv("all_tools_stats.tsv", sep='\t', index=False)
print("✅ Combined stats written to all_tools_stats.tsv")
