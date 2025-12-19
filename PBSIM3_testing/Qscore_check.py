from Bio import SeqIO
import gzip
import hashlib

def get_q_hashes(fastq):
    hashes = {}
    with gzip.open(fastq, "rt") as handle:
        for rec in SeqIO.parse(handle, "fastq"):
            q = ''.join(chr(q + 33) for q in rec.letter_annotations["phred_quality"])
            h = hashlib.md5(q.encode()).hexdigest()
            hashes.setdefault(h, []).append(len(q))
    return hashes

real_q = get_q_hashes("/home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz")
sim_q  = get_q_hashes("/home/taouk/ONT_read_sim/mapped_sim_reads/Badread.fastq.gz")

shared = set(real_q.keys()) & set(sim_q.keys())
print(f"Exact Q-string matches: {len(shared)}")
