#!/usr/bin/env python3

import pandas as pd
import gzip
from Bio import SeqIO
import re
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

all_dfs = []
csv.field_size_limit(10**7)  # 10 million chars

for tool, (paf_path, fastq_path) in files.items():
    print(f"Processing {tool}...")

    # --- Load FASTQ and calculate read-level stats ---
    fastq_dict = {}
    handle = gzip.open(fastq_path, "rt") if fastq_path.endswith(".gz") else open(fastq_path, "r")
    for record in SeqIO.parse(handle, "fastq"):
        seq = str(record.seq).upper()
        mean_q = sum(record.letter_annotations["phred_quality"]) / len(record)
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        fastq_dict[record.id] = {
            "read_length": len(seq),
            "mean_qscore": mean_q,
            "gc_content": gc
        }
    handle.close()
    
    fastq_df = pd.DataFrame.from_dict(fastq_dict, orient="index").reset_index().rename(columns={"index": "read_name"})
    fastq_df["tool"] = tool

    # --- Load PAF if available ---
    try:
        df_paf = pd.read_csv(
            paf_path, sep='\t', header=None, engine="python",
            on_bad_lines="warn", dtype=str  # force all as str
        )

        # convert numeric columns
        for c in [1, 2, 3, 6, 7, 8, 9, 10, 11]:
            df_paf[c] = pd.to_numeric(df_paf[c], errors='coerce')

        df_paf["read_name"] = df_paf[0].astype(str)
        df_paf["read_length_paf"] = df_paf[1]
        df_paf["unaligned_start"] = df_paf[2]
        df_paf["unaligned_end"] = df_paf[1] - df_paf[3]
        df_paf["actual_identity"] = (df_paf[9] / df_paf[10]).astype(float)

        # parse CIGAR
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

        for _, row in df_paf.iterrows():
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

        df_paf["sub_rate"] = sub_rate_list
        df_paf["ins_rate"] = ins_rate_list
        df_paf["del_rate"] = del_rate_list

        df_paf = df_paf[["read_name","unaligned_start","unaligned_end","actual_identity","sub_rate","ins_rate","del_rate"]]

        # --- Merge FASTQ and PAF ---
        merged_df = pd.merge(fastq_df, df_paf, on="read_name", how="left")
        merged_df["aligned"] = merged_df["actual_identity"].notna()

    except FileNotFoundError:
        merged_df = fastq_df.copy()
        merged_df["unaligned_start"] = None
        merged_df["unaligned_end"] = None
        merged_df["actual_identity"] = None
        merged_df["sub_rate"] = None
        merged_df["ins_rate"] = None
        merged_df["del_rate"] = None
        merged_df["aligned"] = False

    all_dfs.append(merged_df)

# --- Combine all tools ---
combined = pd.concat(all_dfs, ignore_index=True)

# --- Columns to write ---
output_cols = ["tool","read_name","read_length","unaligned_start","unaligned_end","actual_identity","sub_rate","ins_rate","del_rate","aligned","mean_qscore","gc_content"]

combined[output_cols].to_csv("all_tools_stats.tsv", sep="\t", index=False)
print("✅ Combined stats written to all_tools_stats.tsv")
