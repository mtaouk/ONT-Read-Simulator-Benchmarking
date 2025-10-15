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
mkdir NanoSim
cd NanoSim

# was getting a weird error because header names were too long
zcat ../reads/Enterobacter_hormaechei_SAMN31246718_sup5.2.0.fastq.gz \
| awk 'NR%4==1{sub(/^@.*/,"@"NR/4)}1' \
| gzip > ../reads/Enterobacter_hormaechei_SAMN31246718_shortnames.fastq.gz

# run genome characterisation stage
/home/taouk/bin/NanoSim/src/read_analysis.py genome -i ../reads/Enterobacter_hormaechei_SAMN31246718_shortnames.fastq.gz -rg ../reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -o training/training -t 30

# run read simulation stage
/home/taouk/bin/NanoSim/src/simulator.py genome -rg ../reads/Enterobacter_hormaechei_SAMN31246718_reference.fasta -o simulation/simulated -x 100 --seed 3 -c training/training -t 32


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

cd /home/taouk/ONT_read_sim
mkdir LongISLND
cd LongISLND

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
/home/taouk/bin/longislnd/simulate.py --movie_id ONT --read_type fastq --model_dir longislnd_model --fasta reference.fasta --coverage 100

cat out/*.fq | gzip > longislnd_reads.fastq.gz

# Clean up
rm -r out real_reads.fastq.bam* reference.fasta.fai