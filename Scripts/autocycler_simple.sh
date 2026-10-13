#!/usr/bin/env bash

# This script is a wrapper for running a simple Autocycler assembly with no read subsampling step.

# Usage:
#   autocycler_simple.sh <read_fastq> <threads> <jobs>

# Copyright 2026 Ryan Wick (rrwick@gmail.com)
# Licensed under the GNU General Public License v3.
# See https://www.gnu.org/licenses/gpl-3.0.html.

# Ensure script exits on error.
set -e

# Get arguments.
reads=$1                 # input reads FASTQ
threads=$2               # threads per job
jobs=$3                  # number of simultaneous jobs

# Input assembly jobs that exceed this time limit will be killed
max_time="4h"

# Validate input parameters
if [[ -z "$reads" || -z "$threads" || -z "$jobs" ]]; then
    echo "Usage: $0 <read_fastq> <threads> <jobs> [read_type]" 1>&2
    exit 1
fi
if [[ ! -f "$reads" ]]; then
    echo "Error: Input file '$reads' does not exist." 1>&2
    exit 1
fi
if (( threads > 100 )); then threads=100; fi  # Autocycler won't work with more than 100 threads

mkdir -p assemblies
rm -f assemblies/jobs.txt
for assembler in raven myloasm miniasm flye metamdbg; do
    echo "autocycler helper $assembler --reads $reads --out_prefix assemblies/$assembler --threads $threads --min_depth_rel 0.1" >> assemblies/jobs.txt
done
set +e
nice -n 19 parallel --jobs "$jobs" --joblog assemblies/joblog.tsv --results assemblies/logs --timeout "$max_time" < assemblies/jobs.txt
set -e

autocycler compress -i assemblies -a autocycler_out -t "$threads" --max_contigs 1000 2>> autocycler.stderr

autocycler cluster -a autocycler_out --max_contigs 1000 2>> autocycler.stderr

shopt -s nullglob
clusters=(autocycler_out/clustering/qc_pass/cluster_*)
if (( ${#clusters[@]} == 0 )); then
    echo "No clusters passed QC" >> autocycler.stderr
    exit 1
fi

for c in "${clusters[@]}"; do
    autocycler trim -c "$c" -t "$threads" 2>> autocycler.stderr
    autocycler resolve -c "$c" 2>> autocycler.stderr
done

autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa -r "$reads" -t "$threads" 2>> autocycler.stderr
