#!/bin/bash

cd ONT_read_sim
mkdir reads
cd reads

# Extract lengths and compute N50
zcat Enterobacter_hormaechei_SAMN31246718_shortnames.fastq.gz | awk 'NR%4==2 {print length($0)}' | sort -nr > read_lengths.txt

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


#get length stats quick
zcat /home/taouk/ONT_read_sim/lrsim/lrsim_reads.fq.gz \
  | awk 'NR%4==2 {
      L = length($0)
      n++
      sum += L
      if (L < min || min == "") min = L
      if (L > max) max = L
    }
    END {print "reads:",n, "mean:",sum/n, "min:",min, "max:",max}'
#







### NanoSim
# installation
git clone https://github.com/bcgsc/NanoSim.git
conda create --name nanosim python=3.7
conda activate nanosim
mamba install scikit-learn=0.22.1 six samtools pysam pybedtools minimap2 joblib htseq genometools-genometools piecewise-regression last numpy=1.21.5 scipy=1.7.3 sam2pairwise
conda install pip
conda install -c conda-forge regex

# setup
mkdir ONT_read_sim/NanoSim
cd ONT_read_sim/NanoSim

# was getting a weird error because header names were too long
zcat ../reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz \
| awk 'NR%4==1{sub(/^@.*/,"@"NR/4)}1' \
| gzip > ../reads/Enterobacter_hormaechei_SAMN31246718_shortnames.fastq.gz

# run genome characterisation stage
/home/taouk/bin/NanoSim/src/read_analysis.py genome -i ../reads/Enterobacter_hormaechei_SAMN31246718_shortnames.fastq.gz -rg ../reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -o training/training -t 30 --fastq

#run read simulation stage
/home/taouk/bin/NanoSim/src/simulator.py genome -rg ../reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -o simulation/simulated -x 100 --seed 3 -c training/training -t 32 --fastq -med 12000 -sd 0.40

gzip simulation/simulated_aligned_reads.fastq > simulation/simulated_aligned_reads.fastq.gz
gzip simulation/simulated_unaligned_reads.fastq > simulation/simulated_unaligned_reads.fastq.gz

# Compute N50
zcat simulation/simulated_aligned_reads.fastq.gz | awk 'NR%4==2 {print length($0)}' | sort -nr | awk '{
    sum+=$1; lengths[NR]=$1
} END {
    target=sum/2; running=0;
    for(i=1;i<=NR;i++){
        running+=lengths[i];
        if(running>=target){print "N50:", lengths[i]; break}
    }
}'

#run read simulation stage without specifying (hope it goes off of reads)
/home/taouk/bin/NanoSim/src/simulator.py genome -rg ../reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -o simulation_15122025/simulated -x 100 --seed 3 -c training/training -t 32 --fastq
gzip simulated_aligned_reads.fastq > nanosim_reads.fastq.gz












### LongISLND
# Create a conda environment
mamba create -n longislnd python=2.7 openjdk maven
conda activate longislnd

# Clone LongISLND
git clone https://github.com/bioinform/longislnd

# Fix the URL in the build script:
cd longislnd
sed -i 's|http://www.hdfgroup.org/ftp/HDF5/prev-releases/HDF-JAVA/hdf-java-2.11/bin|https://support.hdfgroup.org/ftp/HDF5/prev-releases/HDF-JAVA/hdf-java-2.11/bin|' linux_build.sh

# Build LongISLND
./linux_build.sh

mkdir ONT_read_sim/LongISLND
cd ONT_read_sim/LongISLND

#LongISLND needs read files to be unzipped and in the same directory as the tool is running
zcat /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_shortnames.fastq.gz \
> real_reads.fastq

# Create the alignments and indices
conda create -n mapping python=3.9 minimap2 samtools -c conda-forge -c bioconda 
conda activate mapping

minimap2 -a -x map-ont -t 32 /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta real_reads.fastq | samtools sort > real_reads.fastq.bam
samtools index real_reads.fastq.bam
samtools faidx /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta

# Train a LongISLND model
/home/taouk/bin/longislnd/sample.py --input_suffix fastq.bam --read_type fastq --model_dir longislnd_model --reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta

# Generate LongISLND reads
/home/taouk/bin/longislnd/simulate.py --movie_id ONT --read_type fastq --model_dir longislnd_model --fasta /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --coverage 100
# merge
cat out/*.fq | gzip > longislnd_reads.fastq.gz

# Clean up
rm -r out real_reads.fastq.bam* reference.fasta.fai

# Compute N50
zcat longislnd_reads.fastq.gz | awk 'NR%4==2 {print length($0)}' | sort -nr | awk '{
    sum+=$1; lengths[NR]=$1
} END {
    target=sum/2; running=0;
    for(i=1;i<=NR;i++){
        running+=lengths[i];
        if(running>=target){print "N50:", lengths[i]; break}
    }
}'








### PBSIM3
conda activate pbsim3

mkdir ONT_read_sim/PBSIM3
cd ONT_read_sim/PBSIM3

# simulate reads
pbsim --strategy wgs --method sample --sample /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz --depth 100 --genome /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta

# merge
cat *.fq.gz > pbsm_reads.fastq.gz

# if I just use default settings, its giving an N50 of 250 which isn't right. Even though the reads are 12k. Its reading the reads wrong so I have to specify.
pbsim --strategy wgs --method errhmm --errhmm /home/taouk/ONT_read_sim/PBSIM3/ERRHMM-ONT.model  --genome /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --depth 100 --length-mean 12000 --length-sd 5000 --accuracy-mean 0.99 --difference-ratio 39:24:36 --hp-del-bias 1 --prefix sd

# trying to optimise even more
pbsim --strategy wgs --method errhmm --errhmm /home/taouk/ONT_read_sim/PBSIM3/models/ERRHMM-ONT-HQ.model --qshmm  /home/taouk/ONT_read_sim/PBSIM3/models/QSHMM-ONT-HQ.model --genome /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --depth 100 --length-mean 12000 --length-sd 5000 --length-min 100 --length-max 200000 --accuracy-mean 0.99 --difference-ratio 39:24:36 --hp-del-bias 1 --prefix sd

# merge
cat *.fq.gz > pbsm_reads.fastq.gz

### more optimisation
pbsim --strategy wgs --method errhmm --errhmm /home/taouk/ONT_read_sim/PBSIM3/models/ERRHMM-ONT-HQ.model --qshmm  /home/taouk/ONT_read_sim/PBSIM3/models/QSHMM-ONT-HQ.model --genome /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --depth 100 --length-mean 12000 --length-sd 5000 --length-min 100 --length-max 200000 --accuracy-mean 0.99 --difference-ratio 39:24:36 --hp-del-bias 1 --prefix sd


### try AGAIN to get Q scores
pbsim --strategy wgs --method qshmm --errhmm /home/taouk/ONT_read_sim/PBSIM3/models/ERRHMM-ONT-HQ.model --qshmm  /home/taouk/ONT_read_sim/PBSIM3/models/QSHMM-ONT-HQ.model --genome /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --depth 100 --accuracy-mean 0.99 --difference-ratio 39:24:36 --hp-del-bias 1 --pass-num 1 --fastq --prefix sd

### try sampling from reads again
pbsim --strategy wgs --method sample --sample /home/taouk/ONT_read_sim/LongISLND/real_reads.fastq --genome/home/taouk/NGtransmission/NCCP11945.fa --depth 100  --prefix ont_sample --difference-ratio 39:24:36

pbsim --strategy wgs --method sample --sample /home/taouk/ONT_read_sim/LongISLND/real_reads.fastq --genome/home/taouk/NGtransmission/NCCP11945.fa --depth 100  --prefix ont_sample --difference-ratio 39:24:36


# Compute N50
zcat pbsm_reads.fastq.gz | awk 'NR%4==2 {print length($0)}' | sort -nr | awk '{
    sum+=$1; lengths[NR]=$1
} END {
    target=sum/2; running=0;
    for(i=1;i<=NR;i++){
        running+=lengths[i];
        if(running>=target){print "N50:", lengths[i]; break}
    }
}'








### simON-reads
conda activate simON-reads

mkdir ONT_read_sim/simON-reads
cd ONT_read_sim/simON-reads

# simulate reads
# this is the simplest tool. It will give ~90–95% accuracy, cannot be changed. I can only specify the number of reads total, not coverage or length.
simON_reads.py -i /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -n 10000 --enable_homopolymer_error --enable_sequencing_error > simON_reads.fastq

#zip
gzip simON_reads.fastq > simON_reads.fastq.gz








### Badread
conda activate badread

mkdir ONT_read_sim/Badread
cd ONT_read_sim/Badread

# align reads to reference
minimap2 -x map-ont  -c /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz -t 32 > real_reads.paf

# train for errors
badread error_model --reads /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz --reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --alignment real_reads.paf > ont_model.err

# train for qscore
badread qscore_model --reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --reads /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz --alignment real_reads.paf > ont_model.qsc

# simulate
badread simulate --reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --error_model ont_model.err --qscore_model ont_model.qsc --quantity 50x --length 15000,10000 --identity 94,2.5 | gzip > simulated_reads.fastq.gz

# simulate without parameters
badread simulate --reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --error_model ont_model.err --qscore_model ont_model.qsc --quantity 50x | gzip > simulated_reads_15122025.fastq.gz









### lrsim
conda activate lrsim_env

bash taouk/bin/lrsim/unzipmodels.sh

mkdir ONT_read_sim/lrsim
cd ONT_read_sim/lrsim

#simulation
lrsim -t 32 -m /home/taouk/bin/lrsim/models/HG002_ONT_UL.lrsm -d 100 --error=0.01 --rvarianceratio=8 --bfvarianceratio=1 --blocksize=200 --eratio 50:30:20 --fixedreadlength=12000 --lengthfloatratio=0.2  /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta | gzip > lrsim_reads.fq.gz









### simlord
# installation
conda create -n simlord python=3.10 pip numpy scipy cython 
conda activate simlord
pip install pysam 
conda install -c bioconda dinopy 
pip install simlord

mkdir ONT_read_sim/simlord
cd ONT_read_sim/simlord

#simulation
simlord --read-reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -n 50000 -fl 15000 -pi 0.08 -pd 0.08 -ps 0.02 --no-sam simlord_reads

#zip
gzip simlord_reads.fastq > simlord_reads.fastq.gz

# ONT linke
simlord --read-reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --coverage 100 --lognorm-readlength 0.8 0 12000 --prob-ins 0.04 --prob-del 0.08 --prob-sub 0.01 --no-sam simlord_ont_like

#zip
gzip simlord_ont_like.fastq > simlord_ont_reads.fastq.gz

# More ONT accurate 15122025
simlord --read-reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --coverage 100 --sample-readlength-from-fastq /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz --prob-ins 0.03 --prob-del 0.06 --prob-sub 0.004 --probability-threshold 0.18 --no-sam --gzip simlord_15122025_like




### Squigilator
conda activate squigulator_env

mkdir ONT_read_sim/squigulator
cd ONT_read_sim/squigulator

squigulator /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -x dna-r10-prom -o squigulator_reads.blow5 -f 30












### Data gathering
mkdir -p mapped_sim_reads
cd mapped_sim_reads

conda activate mapping

# subset the real reads to make them about 100x coverage
seqtk sample -s100 Enterobacter_hormaechei_SAMN31246718_shortnames.fastq.gz 50000 > real_reads_subset.fastq.gz


# copying over whichever read set I want to use for the data extraction, because sometimes I edit the upstream code and change the names but I want the "final" ones to be in one place for clarity and easy coding

cp /home/taouk/ONT_read_sim/Badread/simulated_reads_15122025.fastq.gz Badread.fastq.gz
cp /home/taouk/ONT_read_sim/LongISLND/longislnd_reads.fastq.gz LongISLND.fastq.gz
cp /home/taouk/ONT_read_sim/lrsim/lrsim_reads.fq.gz lrsim.fastq.gz
cp /home/taouk/ONT_read_sim/NanoSim/simulation_15122025/simulated_aligned_reads.fastq.gz NanoSim.fastq.gz
cp /home/taouk/ONT_read_sim/PBSIM3/18122025/pbsm_reads.fastq.gz PBSIM3.fastq.gz
cp /home/taouk/ONT_read_sim/simlord/simlord_15122025_like.fastq.gz simlord.fastq.gz
cp /home/taouk/ONT_read_sim/reads/real_reads_subset.fastq.gz real.fastq.gz


### mapping reads to reference:
```
ref="/home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta"

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
```

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






###### PBSIM3 test using gono reads

minimap2 -t 8 -c -eqx /home/taouk/NGtransmission/NCCP11945.fa  /home/taouk/ONT_read_sim/PBSIM3/test18122025/gono_reads.fastq.gz > gono_reads.paf

grep -v "tp:A:S" gono_reads.paf > gono_primary.paf

python TEST_extract_read_stats.py