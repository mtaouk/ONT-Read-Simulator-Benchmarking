#!/bin/bash

OUTDIR="/home/taouk/ONT_read_sim/unicycler_assemblies"
mkdir -p "$OUTDIR"

declare -A reads=(
  [Badread]="/home/taouk/ONT_read_sim/mapped_sim_reads/Badread.fastq.gz"
  [LongISLND]="/home/taouk/ONT_read_sim/mapped_sim_reads/LongISLND.fastq.gz"
  [lrsim]="/home/taouk/ONT_read_sim/mapped_sim_reads/lrsim.fastq.gz"
  [NanoSim]="/home/taouk/ONT_read_sim/mapped_sim_reads/NanoSim.fastq.gz"
  [PBSIM3]="/home/taouk/ONT_read_sim/mapped_sim_reads/PBSIM3.fastq.gz"
  [simlord]="/home/taouk/ONT_read_sim/mapped_sim_reads/simlord.fastq.gz"
  [real]="/home/taouk/ONT_read_sim/mapped_sim_reads/real.fastq.gz"
)

for name in "${!reads[@]}"; do
  echo "Assembling $name..."
  unicycler \
    -l "${reads[$name]}" \
    -o "$OUTDIR/$name" \
    --mode bold \
    --threads 30 &
done

wait