#!/bin/bash

cd ONT_read_sim
mkdir reads
cd reads


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

# run read simulation stage
/home/taouk/bin/NanoSim/src/simulator.py genome -rg ../reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -o simulation_fastq/simulated -x 100 --seed 3 -c training/training -t 32 --fastq

gzip simulation_fastq/simulated_aligned_reads.fastq > simulation_fastq/simulated_aligned_reads.fastq.gz
gzip simulation_fastq/simulated_unaligned_reads.fastq > simulation_fastq/simulated_unaligned_reads.fastq.gz

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


### PBSIM3
conda activate pbsim3

mkdir ONT_read_sim/PBSIM3
cd ONT_read_sim/PBSIM3

# simulate reads
pbsim --strategy wgs --method sample --sample /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz --depth 100 --genome /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta
# merge
cat *.fq.gz > pbsm_reads.fastq.gz


### simON-reads
conda activate simON-reads

mkdir ONT_read_sim/simON-reads
cd ONT_read_sim/simON-reads

# simulate reads
simON_reads.py -i /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta > simON_reads.fastq
gzip simON_reads.fastq > simON_reads.fastq.gz


### Badread
conda activate badread

mkdir ONT_read_sim/Badread
cd ONT_read_sim/Badread

badread simulate --reference /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta --quantity 50x --glitches 10000,10,10  --junk_reads 0.1 --random_reads 0.1 --chimeras 0.1 --identity 20,3 | gzip > badread_reads.fastq.gz


### lrsim
conda activate lrsim_env

bash taouk/bin/lrsim/unzipmodels.sh

mkdir ONT_read_sim/lrsim
cd ONT_read_sim/lrsim

#simulation
lrsim -t 32 -m /home/taouk/bin/lrsim/models/HG002_ONT_UL.lrsm -d 100 /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta > lrsim_reads.fq
gzip lrsim_reads.fq > lrsim_reads.fq.gz


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
gzip simlord_reads.fastq > simlord_reads.fastq.gz


### Squigilator
conda activate squigulator_env

mkdir ONT_read_sim/squigulator
cd ONT_read_sim/squigulator

squigulator /home/taouk/ONT_read_sim/reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -x dna-r10-prom -o squigulator_reads.blow5 -f 30