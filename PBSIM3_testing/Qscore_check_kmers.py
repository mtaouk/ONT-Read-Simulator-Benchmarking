import gzip
import random
from Bio import SeqIO

# >>>>>>>>>>>> INPUT FILES GO HERE <<<<<<<<<<<<
REAL_FASTQ = "/home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz"
SIM_FASTQ  = "/home/taouk/ONT_read_sim/PBSIM3/18122025/pbsm_reads.fastq.gz"
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def load_qs(fastq):
    qs = []
    with gzip.open(fastq, "rt") as handle:
        for rec in SeqIO.parse(handle, "fastq"):
            q = ''.join(chr(q + 33) for q in rec.letter_annotations["phred_quality"])
            qs.append(q)
    return qs

def q_kmers(q, k=20):
    return {q[i:i+k] for i in range(len(q) - k + 1)}

# Load inputs
real_qs = load_qs(REAL_FASTQ)
sim_qs  = load_qs(SIM_FASTQ)

# Build k-mer universe from real reads
real_kmers = set()
for q in real_qs:
    real_kmers |= q_kmers(q)

# Test random simulated reads
total = 0
found = 0

for q in random.sample(sim_qs, 200):
    for kmer in q_kmers(q):
        total += 1
        if kmer in real_kmers:
            found += 1

print(f"{found/total:.3f} fraction of simulated Q k-mers found in real reads")
