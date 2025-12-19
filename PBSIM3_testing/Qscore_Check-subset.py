from Bio import SeqIO
import gzip
import hashlib


def load_q_strings(fastq, max_reads=1000):
    qs = []
    with gzip.open(fastq, "rt") as handle:
        for i, rec in enumerate(SeqIO.parse(handle, "fastq")):
            if i >= max_reads:
                break
            q = ''.join(chr(q + 33) for q in rec.letter_annotations["phred_quality"])
            qs.append(q)
    return qs

real_qs = load_q_strings("/home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz")
sim_qs  = load_q_strings("/home/taouk/ONT_read_sim/PBSIM3/18122025/pbsm_reads.fastq.gz")

hits = 0
for sq in sim_qs:
    for rq in real_qs:
        if sq in rq:
            hits += 1
            break

print(f"Sim Q strings that are substrings of real Qs: {hits}/{len(sim_qs)}")
