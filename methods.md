# ONT-Read-Simulator-Benchmarking

Moving the *Listeria innocua* genome and reads from Ryan's directory
into mine:

```{bash}
mv /home/wickr/ONT_read_sim/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz

mv /home/wickr/ONT_read_sim/Listeria_innocua_SAMN46906078_reference.fasta Listeria_innocua_SAMN46906078_reference.fasta
```

First, I extracted some basic statistics for the reads:

```{bash}
# Extract lengths and compute N50
zcat Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz | awk 'NR%4==2 {print length($0)}' | sort -nr > read_lengths.txt

# Compute N50
awk '{
    sum+=$1; 
    lengths[NR]=$1
} END {
    target=sum/2; 
    running=0; 
    for(i=1;i<=NR;i++){ 
        running+=lengths[i]; 
        if(running>=target){print "N50:", lengths[i]; break} 
    }
}' read_lengths.txt
#N50: 8923

#get length stats quick
zcat Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz \
  | awk 'NR%4==2 {
      L = length($0)
      n++
      sum += L
      if (L < min || min == "") min = L
      if (L > max) max = L
    }
    END {print "reads:",n, "mean:",sum/n, "min:",min, "max:",max}'
#reads: 680482 mean: 5588.16 min: 5 max: 826029
```

## Read Simulation

The following is the code I used for each read simulator tool to
generate the reads, including installation and any preparation steps.
YAML files for each conda environment used in the following code is
available in Environments

### Nanosim

```{bash}
conda activate nanosim

# was getting a weird error because header names were too long
zcat /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz | awk 'NR%4==1{sub(/^@.*/,"@"NR/4)}1' | gzip > /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_shortnames.fastq.gz

# run genome characterisation stage
/home/taouk/bin/NanoSim/src/read_analysis.py genome -i /home/taouk/ONT_read_sim/longislnd/real_reads.fastq -rg /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta -o training2/training -t 60 --fastq

#run read simulation stage
/home/taouk/bin/NanoSim/src/simulator.py genome -rg /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta -o simulation2/simulated -x 100 --seed 3 -c training2/training -t 32 --fastq -hp 5

#compress final read files
gzip simulation/simulated_aligned_reads.fastq > simulation/simulated_aligned_reads.fastq.gz
gzip simulation/simulated_unaligned_reads.fastq > simulation/simulated_unaligned_reads.fastq.gz
```

### LongISLND

```{bash}
# Create the alignments and indices
conda activate mapping

minimap2 -a -x map-ont -t 32 /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz | samtools sort > real_reads.fastq.bam
samtools index real_reads.fastq.bam
samtools faidx /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta
zcat /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_shortnames.fastq.gz \
> real_reads.fastq

conda activate longislnd

# Train a LongISLND model
/home/taouk/bin/longislnd/sample.py --input_suffix fastq.bam --read_type fastq --model_dir longislnd_model --reference /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta

# Generate LongISLND reads
/home/taouk/bin/longislnd/simulate.py --movie_id ONT --read_type fastq --model_dir longislnd_model --fasta /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta --coverage 100

# merge
cat out/*.fq | gzip > longislnd_reads.fastq.gz

# Clean up
rm -r out real_reads.fastq.bam* reference.fasta.fai
```

### PBSIM3

```{bash}
conda activate pbsim3

#run
pbsim --strategy wgs --method sample --sample /home/taouk/ONT_read_sim/longislnd/real_reads.fastq --genome /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta --depth 100  --prefix ont_sample --difference-ratio 39:24:36

# merge
cat *.fq.gz > pbsim_reads.fastq.gz
```

### Badread

```{bash}
# align reads to reference
conda activate mapping

minimap2 -x map-ont  -c /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz -t 32 > real_reads.paf

conda activate badread

# train for errors
badread error_model --reads /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz --reference /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta --alignment real_reads.paf > ont_model.err

# train for qscore
badread qscore_model --reads /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz --reference /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta --alignment real_reads.paf > ont_model.qsc

# simulate without perameters
badread simulate --reference /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta --error_model ont_model.err --qscore_model ont_model.qsc --quantity 100x --length 8500,8500 | gzip > badread_reads.fastq.gz 

```

### LRSim

```{bash}
### lrsim
conda activate lrsim_env

#simulation
lrsim -t 32 -m /home/taouk/bin/lrsim/models/HG002_ONT_UL.lrsm -d 100 --error=0.01 --rvarianceratio=8 --bfvarianceratio=1 --blocksize=200 --eratio 50:30:20 --fixedreadlength=12000 --lengthfloatratio=0.2  /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta | gzip > lrsim_reads.fq.gz

```

### Simlord

```{bash}
conda activate simlord

# make reads
simlord --read-reference /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta --coverage 100 --sample-readlength-from-fastq /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz --max-passes 12 --prob-ins 0.05 --prob-del 0.06 --prob-sub 0.02 --probability-threshold 0.18 --sqrt-params 0.5 0.25 --norm-params 0 0.3 --gzip simlord_reads

```

## Data Gathering

The following code is how I extracted statistics from the simulated
reads for analysis.

```{bash}
mkdir -p mapped_sim_reads
cd mapped_sim_reads

conda activate seqtk

# subset the real reads to make them about 100x coverage
seqtk sample -s100 /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz 50000 > real_reads_subset.fastq
gzip real_reads_subset.fastq

# copying over whichever read set I want to use for the data extraction, because sometimes I edit the upstream code and change the names but I want the "final" ones to be in one place for clarity and easy coding

cp /home/taouk/ONT_read_sim/badread/badread_reads.fastq.gz Badread.fastq.gz
cp /home/taouk/ONT_read_sim/longislnd/longislnd_reads.fastq.gz LongISLND.fastq.gz
cp /home/taouk/ONT_read_sim/lrsim/lrsim_reads.fq.gz lrsim.fastq.gz
cp /home/taouk/ONT_read_sim/nanosim/simulation/simulated_aligned_reads.fastq.gz NanoSim.fastq.gz
cp /home/taouk/ONT_read_sim/pbsim/pbsim_reads.fastq.gz PBSIM3.fastq.gz
cp /home/taouk/ONT_read_sim/simlord/simlord_reads.fastq.gz simlord.fastq.gz
mv real_reads_subset.fastq.gz real.fastq.gz

conda activate mapping

### mapping reads to reference:
ref="/home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta"

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
  echo "Aligning $name..."
  minimap2 -t 8 -c -eqx "$ref" "${reads[$name]}" > "${name}.paf" &
done

wait

# clean .pafs so that only the primary alignments for each read are kept
grep -v "tp:A:S" Badread.paf > Badread_primary.paf
grep -v "tp:A:S" LongISLND.paf > LongISLND_primary.paf
grep -v "tp:A:S" lrsim.paf > lrsim_primary.paf
grep -v "tp:A:S" NanoSim.paf > NanoSim_primary.paf
grep -v "tp:A:S" PBSIM3.paf > PBSIM3_primary.paf
grep -v "tp:A:S" simlord.paf > simlord_primary.paf
grep -v "tp:A:S" real.paf > real_primary.paf

# run code that gets the stats
conda activate python_tools

python extract_read_stats.py

```

extract_read_stats.py is available in Scripts. The results:
all_tools_stats.tsv is available in Results.

### This code is for making assemblies and then extracting stats

```{bash}

##### Assemblies 

conda activate flye

cd flye_assemblies

bash assemblies.sh

# get stats for assemblies
quast /home/taouk/ONT_read_sim/flye_assemblies/*/assembly.fasta -r /home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta --large -o quast_flye_ref -t 8
```

assemblies.sh is available in Scripts. The QUAST results html is
available in Results.

### This code is for base level statistics

```{bash}
ref="/home/taouk/ONT_read_sim/reads/Listeria_innocua_SAMN46906078_reference.fasta"

declare -A reads=(
  [Badread]="Badread.fastq.gz"
  [LongISLND]="LongISLND.fastq.gz"
  [lrsim]="lrsim.fastq.gz" 
  [NanoSim]="NanoSim.fastq.gz" 
  [PBSIM3]="PBSIM3.fastq.gz"
  [simlord]="simlord.fastq.gz" 
  [real]="real.fastq.gz" 
)

for name in "${!reads[@]}"; do
  echo "Calculating $name..."
  python per_base_qscores.py "${reads[$name]}" "$ref" "${name}.paf" > $name.qscore_counts.tsv
done

wait
```

per_base_qscores.py is available in Scripts. The outputs are available
in Results.

### Extracting some basic read level statistics

```{bash}
### read stats summary:
cd /home/taouk/ONT_read_sim/mapped_sim_reads
seqkit stats -a \
  real.fastq.gz \
  Badread.fastq.gz \
  LongISLND.fastq.gz \
  lrsim.fastq.gz \
  NanoSim.fastq.gz \
  PBSIM3.fastq.gz \
  simlord.fastq.gz \
  > length_stats.tsv
```

length_stats.tsv is available in Results.

Everything was run again for an **Enterobacter** genome. All scripts and
results for that analysis can be found in the Enterobacter folder.
