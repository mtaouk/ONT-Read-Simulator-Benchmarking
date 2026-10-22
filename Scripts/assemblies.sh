#!/usr/bin/env bash

# This script runs Flye, metaMDBG, miniasm, Myloasm and Raven assemblies, using Autocycler helper.

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
