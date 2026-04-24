# ONT-Read-Simulator-Benchmarking


## Basecalling

This is the sample: [SAMN46906078](https://www.ncbi.nlm.nih.gov/biosample/SAMN46906078). And while there are reads on SRA, I'd like to use freshly basecalled reads instead to have data as up-to-date as possible. At the time of writing, Dorado v1.1.1 and v5.2.0 basecalling models are the most current.

I did the basecalling on [Spartan](https://dashboard.hpc.unimelb.edu.au) because it has big GPUs. We may just use sup reads, but I'll create fast and hac as well, just in case.
```bash
cd /data/scratch/projects/punim1894/O2024-029
sbatch --job-name=dorado_basecaller --time=10:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 -p gpu-h100 --gres=gpu:1 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado basecaller --kit-name SQK-RBK114-96 fast pod5 > dorado_v1.1.1_fast5.2.0.bam"
sbatch --job-name=dorado_basecaller --time=20:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 -p gpu-h100 --gres=gpu:1 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado basecaller --kit-name SQK-RBK114-96 hac pod5 > dorado_v1.1.1_hac5.2.0.bam"
sbatch --job-name=dorado_basecaller --time=40:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 -p gpu-h100 --gres=gpu:1 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado basecaller --kit-name SQK-RBK114-96 sup pod5 > dorado_v1.1.1_sup5.2.0.bam"
```

Summary table:
```bash
cd /data/scratch/projects/punim1894/O2024-029

sbatch --job-name=dorado_summary --time=1:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado summary dorado_v1.1.1_fast5.2.0.bam > dorado_v1.1.1_fast5.2.0.tsv"
sbatch --job-name=dorado_summary --time=1:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado summary dorado_v1.1.1_hac5.2.0.bam > dorado_v1.1.1_hac5.2.0.tsv"
sbatch --job-name=dorado_summary --time=1:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado summary dorado_v1.1.1_sup5.2.0.bam > dorado_v1.1.1_sup5.2.0.tsv"
```

Demultiplex:
```bash
cd /data/scratch/projects/punim1894/O2024-029

sbatch --job-name=dorado_demux --time=4:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado demux --output-dir dorado_v1.1.1_fast5.2.0 --no-classify -t 8 dorado_v1.1.1_fast5.2.0.bam"
sbatch --job-name=dorado_demux --time=4:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado demux --output-dir dorado_v1.1.1_hac5.2.0 --no-classify -t 8 dorado_v1.1.1_hac5.2.0.bam"
sbatch --job-name=dorado_demux --time=4:00:00 --ntasks=1 --mem=64000 --cpus-per-task=8 --wrap "~/programs/dorado-1.1.1-linux-x64/bin/dorado demux --output-dir dorado_v1.1.1_sup5.2.0 --no-classify -t 8 dorado_v1.1.1_sup5.2.0.bam"
```

Grab the _Listeria innocua_ reads from barcode 05 - will be used as training data and as the real reads:
```bash
cd /data/scratch/projects/punim1894/O2024-029

samtools cat -o Listeria_innocua_SAMN46906078_fast5.2.0.bam dorado_v1.1.1_fast5.2.0/*_barcode05.bam
samtools cat -o Listeria_innocua_SAMN46906078_hac5.2.0.bam dorado_v1.1.1_hac5.2.0/*_barcode05.bam
samtools cat -o Listeria_innocua_SAMN46906078_sup5.2.0.bam dorado_v1.1.1_sup5.2.0/*_barcode05.bam
```

Convert to FASTQ:
```bash
cd /data/scratch/projects/punim1894/O2024-029

samtools fastq -T '*' Listeria_innocua_SAMN46906078_fast5.2.0.bam | tr '\t' ' ' | paste - - - - | shuf | tr '\t' '\n' | pigz -p16 > Listeria_innocua_SAMN46906078_fast5.2.0.fastq.gz
samtools fastq -T '*' Listeria_innocua_SAMN46906078_hac5.2.0.bam | tr '\t' ' ' | paste - - - - | shuf | tr '\t' '\n' | pigz -p16 > Listeria_innocua_SAMN46906078_hac5.2.0.fastq.gz
samtools fastq -T '*' Listeria_innocua_SAMN46906078_sup5.2.0.bam | tr '\t' ' ' | paste - - - - | shuf | tr '\t' '\n' | pigz -p16 > Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz
```

By using `shuf`, I'm ensuring there isn't any structure in the FASTQ file. This means that taking the first n reads from the file will be an unbiased subset.


## Prep training and test data

Get the reference genome and add circularity indicators to the reference FASTA headers (Badread will use this, other tools may not):
```bash
cd ~/2025-08_ONT_read_simulator_benchmark
cp ~/2025-04_Autocycler_paper/Listeria_innocua/reference.fasta reference.fasta
sed -i "s/chromosome/chromosome circular=true/" reference.fasta
sed -i "s/plasmid/plasmid circular=true/" reference.fasta
samtools faidx reference.fasta
```

Copy the reads from Spartan to the MDU servers:
```bash
cd ~/2025-08_ONT_read_simulator_benchmark
mkdir real_reads
cd real_reads
scp spartan:/data/scratch/projects/punim1894/O2024-029/Listeria_innocua_SAMN46906078_sup5.2.0.fastq.gz real_reads.fastq.gz
```

Some basic read stats:
```bash
cd ~/2025-08_ONT_read_simulator_benchmark/real_reads
seqtk size real_reads.fastq.gz
seqkit stats real_reads.fastq.gz
fast_count real_reads.fastq.gz  # for N50
```

Results:
* Read count: 680,482 bp
* Total length: 3,802,641,279 bp
* Mean length: 5,588 bp
* N50 length: 8,923 bp
* Max length: 826,029 bp (almost certainly a junk read)

Split the real reads into training and test sets. Since they have been shuffled, I can just take the first 53093 reads for the test set. This value was chosen (via some trial and error) to be close to 297,254,500 bp which is 100x. The remaining reads go in the training set.
```bash
zcat real_reads.fastq.gz | sed 's/[[:space:]].*$//' | paste - - - - | head -n 53093 | tr '\t' '\n' > test_reads.fastq
zcat real_reads.fastq.gz | sed 's/[[:space:]].*$//' | paste - - - - | tail -n +53094 | tr '\t' '\n' > training_reads.fastq
```

Produce a short-names version of the reads (avoid a problem later with NanoSim):
```bash
cd ~/2025-08_ONT_read_simulator_benchmark/real_reads
cat training_reads.fastq | awk 'NR%4==1{sub(/^@.*/,"@"int(NR/4+1))}1' > training_reads_shortnames.fastq
```




## Read characterisation

In this step, I gather a bunch of stats about the training reads that I can use in the read simulation commands.

Some basic read stats:
```bash
cd ~/2025-08_ONT_read_simulator_benchmark/real_reads
seqtk size training_reads.fastq
seqkit stats training_reads.fastq
fast_count training_reads.fastq
seqkit fx2tab -n -l training_reads.fastq | csvtk -H -t summary -f 2:mean,2:stdev
```

Align to the reference:
```bash
cd ~/2025-08_ONT_read_simulator_benchmark/real_reads

conda activate mapping
minimap2 -c -x map-ont --eqx -t 32 ../reference.fasta training_reads.fastq | grep "tp:A:P" > training_reads.paf
minimap2 -a -x map-ont -t 32 ../reference.fasta ../real_reads/training_reads.fastq | samtools view -u -F 0x904 | samtools sort > training_reads.bam
samtools index training_reads.bam
```

Mean and stdev for identity and Q-score:
```python
paf_filename = "training_reads.paf"

import math, statistics

identities, qscores = [], []
with open(paf_filename) as f:
    for line in f:
        parts = line.split("\t")
        identity = int(parts[9]) / int(parts[10])
        identities.append(identity)
        if identity < 1.0:
            qscores.append(-10 * math.log10(1.0 - identity))

print(statistics.mean(identities))
print(statistics.stdev(identities))
print(statistics.mean(qscores))
print(statistics.stdev(qscores))
```

substitution:insertion:deletion ratio:
```python
paf_filename = "training_reads.paf"

import re
counts = {"X": 0, "I": 0, "D": 0}

with open(paf_filename, "rt") as f:
    for line in f:
        for part in line.rstrip().split("\t"):
            if part.startswith("cg:Z:"):
                for p in re.findall(r"\d+[IDX=M]", part[5:]):
                    op = p[-1]
                    if op in counts:
                        counts[op] += int(p[:-1])
                break

total = counts["X"] + counts["I"] + counts["D"]
print(f"{int(round(counts['X'] / total * 1000))}:{int(round(counts['I'] / total * 1000))}:{int(round(counts['D'] / total * 1000))}")
```


```bash
samtools stats training_reads.bam | awk -F '\t' '$1=="SN" && $2=="error rate:" {print $3}'
samtools stats -r ../reference.fasta training_reads.bam | awk -F '\t' '$1=="MPC" {for (i=4; i<=NF; i++) s+=$i} END {print s+0}'
samtools stats training_reads.bam | awk -F '\t' '$1=="ID" && $2==1 {print $3; exit}'


S=$(samtools stats -r ../reference.fasta training_reads.bam | awk -F '\t' '$1=="MPC" {for (i=4; i<=NF; i++) s+=$i} END {print s+0}')
I1=$(samtools stats training_reads.bam | awk -F '\t' '$1=="ID" && $2==1 {print $3; exit}')
SUB=$(awk -v s="$S" -v i="$I1" 'BEGIN {print (i>0 ? int(100*s/i + 0.5) : 0)}')
```

Results:
* Read count: 627,389
* Total length: 3,505,375,301 bp
* Read length:
  * mean: 5,587 bp
  * stdev: 6,402 bp
  * N50: 8,919 bp
  * max: 789,802 bp (almost certainly a junk read)
* Read accuracy:
  * error rate mean: 0.03313
  * identity mean: 0.9669
  * identity stdev: 0.0546 (also the stdev for error rate)
  * qscore mean: 19.22
  * qscore stdev: 6.23
* error types:
  * sub:ins:del ratio: 400:189:410
  * multiplying those by the mean error rate of 0.03313 gives:
    * sub rate: 0.0133
    * ins rate: 0.00626
    * del rate: 0.0136





## Install tools

For all of the `conda` commands below, I had conda-forge and bioconda in my conda channels.

Some tools were easy to install just with conda:
```bash
conda create -n mapping minimap2=2.30 samtools=1.22
conda create -n nanosim nanosim=3.2.3
conda create -n pbsim3 pbsim3=3.0.5
conda create -n badread badread=0.4.1
conda create -n simlord simlord=1.0.4
```

LongISLND didn't have a bioconda recipe, so it took a bit more work:
```bash
conda create -n longislnd openjdk=11.0 python=2.7 maven
conda activate longislnd
cd ~/programs
git clone https://github.com/bioinform/longislnd.git
cd longislnd
sed -i 's#http://www.hdfgroup.org/ftp/HDF5/prev-releases/HDF-JAVA/hdf-java-2.11/bin#https://support.hdfgroup.org/ftp/HDF5/releases/HDF-JAVA/hdf-java-2.11/bin#' linux_build.sh
bash linux_build.sh
cp sample.py simulate.py LongISLND.jar ~/miniconda3/envs/longislnd/bin/
mkdir -p ~/miniconda3/envs/longislnd/bin/build
cp -r build/lib ~/miniconda3/envs/longislnd/bin/build/
```

The lrsim on bioconda is a different tool, so installing this one manually:
```bash
conda create -n lrsim bzip2 libcurl libdeflate xz zlib openssl
conda activate lrsim
cd ~/programs
git clone https://github.com/CoREse/lrsim
cd lrsim
git submodule update --init --recursive
sed -i 's|autoreconf -i \&\& ./configure|autoreconf -i \&\& ./configure --disable-bz2|' Makefile
sed -i 's|^HTSLIB_LIBS = .*|HTSLIB_LIBS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib -lz -lm -llzma -lcurl -lpthread -lcrypto -ldeflate|' Makefile
make
cp lrsim extractModel.py "$CONDA_PREFIX/bin/"
```

For each tool, I installed the latest release, with the following exception:
* For LongISLND, the last tagged version is 0.9.5 from 9 years ago, so it seems likely that there will be no more releases. So I cloned from GitHub to get the very latest version, which has a couple little fixes after 0.9.5.





## Read Simulation

The following is the code I used for each read simulator tool to generate the reads, including installation and any preparation steps.

Whenever I had to combine multiple FASTQs together, I shuffle them like this: `tr '\t' ' ' | paste - - - - | shuf | tr '\t' '\n'`. Ensuring the reads are in a random order makes it easier to take samples of the reads later, e.g. I can take the first 1000 reads and they shouldn't be biased.


### Run NanoSim

https://github.com/BirolLab/NanoSim

* `-hp` enables homopolymer length simulation.
* `-k 5` sets the minimum homopolyer length where expansion/contraction will be applied. Chosen because the help text says 'a typical k is 5'
* `-x 100` to get 100x read depth.
* `-t 32` is the thread count.

```bash
cd ~/2025-08_ONT_read_simulator_benchmark
mkdir nanosim
cd nanosim
conda activate nanosim

read_analysis.py genome -i ../real_reads/training_reads_shortnames.fastq -rg ../reference.fasta -o training/training -t 32 --fastq -hp
simulator.py genome -rg ../reference.fasta -c training/training -o simulation/simulation --fastq -hp -k 5 -x 100 -t 32
cat simulation/simulation_*.fastq | tr '\t' ' ' | paste - - - - | shuf | tr '\t' '\n' > nanosim.fastq
```


### LongISLND

https://bioinform.github.io/longislnd

Parameters:
* `--coverage 100` to get 100x read depth.

```bash
cd ~/2025-08_ONT_read_simulator_benchmark
mkdir longislnd
cd longislnd

conda activate mapping
cp ../real_reads/training_reads.fastq .
minimap2 -a -x map-ont -t 32 ../reference.fasta training_reads.fastq | samtools sort > training_reads.fastq.bam
samtools index training_reads.fastq.bam

conda activate longislnd
sample.py --input_suffix fastq.bam --read_type fastq --model_dir longislnd_model --reference ../reference.fasta
simulate.py --movie_id ONT --read_type fastq --model_dir longislnd_model --fasta ../reference.fasta --coverage 100
cat out/*.fq | tr '\t' ' ' | paste - - - - | shuf | tr '\t' '\n' > longislnd.fastq
rm -r out training_reads*  # clean up
```


### PBSIM3

https://github.com/yukiteruono/pbsim3

PBSIM3 has three mode options: `errhmm`, `qshmm` and `sample`
* `errhmm` models errors, but it doesn't generate Q-strings (just `!`).
  * Not ideal for realism.
* `qshmm` models Q-strings and then adds errors based on the Q-scores.
  * Not bad, but the next option is better.
* `sample` copies Q-string from the training reads and then adds errors based on the Q-scores.
  * Better than `qshmm` because the Q-strings are now perfectly realistic (they are real Q-strings).
  * So this one is best, and that's what I'll use.

Parameters:
* `--depth 100` to get 100x read depth.
* `--difference-ratio 400:189:410` is the sub:ins:del ratio. PBSIM3 doesn't care about the scaling here, just the relative values.

```bash
cd ~/2025-08_ONT_read_simulator_benchmark
mkdir pbsim3
cd pbsim3
conda activate pbsim3

pbsim --strategy wgs --method sample --sample ../real_reads/training_reads.fastq --genome ../reference.fasta --depth 100 --prefix ont_sample --difference-ratio 400:189:410
zcat *.fq.gz | tr '\t' ' ' | paste - - - - | shuf | tr '\t' '\n' > pbsim3.fastq
```


### Badread

https://github.com/rrwick/Badread

Parameters:
* `--quantity 100x` to get 100x read depth.
* `--length 5587,6402` is the mean and stdev length
* `--identity 19.22,6.23` is the mean and stdev identity (expressed as Q-scores)

```bash
cd ~/2025-08_ONT_read_simulator_benchmark
mkdir badread
cd badread

conda activate mapping
minimap2 -x map-ont -c ../reference.fasta ../real_reads/training_reads.fastq -t 32 > training_reads.paf

conda activate badread
badread error_model --reads ../real_reads/training_reads.fastq --reference ../reference.fasta --alignment training_reads.paf > error_model
badread qscore_model --reads ../real_reads/training_reads.fastq --reference ../reference.fasta --alignment training_reads.paf > qscore_model
badread simulate --reference ../reference.fasta --error_model error_model --qscore_model qscore_model --quantity 100x --length 5587,6402 --identity 19.22,6.23 > badread.fastq
```


### lrsim

https://github.com/CoREse/lrsim

Parameters:
* `-d 100` to get 100x read depth.
* `--error=0.03313` is from the mean error rate of the training reads.
* `--rvarianceratio 1.648` is the error rate stdev over the mean error rate: 0.0546 / 0.03313. This option is confusingly named, as lrsim treats it as stdev/mean, not variance/mean.
* `--eratio 100:217:212` is ins:del:sub, so I took the `400:189:410` sub:ins:del ratio that I got, rearranged the order and normalised to ins=100, which is what lrsim expects.
* `--nolengthfloat` avoids adding extra variation on top of the trained length distribution.

```bash
cd ~/2025-08_ONT_read_simulator_benchmark
mkdir lrsim
cd lrsim

conda activate lrsim

minimap2 -a -x map-ont -t 32 ../reference.fasta ../real_reads/training_reads.fastq | samtools view -u -F 0x904 | samtools sort > training_reads.bam
samtools index training_reads.bam
samtools stats training_reads.bam | extractModel.py > training.lrsm

lrsim -t 32 -m training.lrsm -d 100 --error=0.03313 --rvarianceratio 1.648 --eratio 100:217:212 --nolengthfloat ../reference.fasta > lrsim.fastq
```


### Simlord

https://bitbucket.org/genomeinformatics/simlord

Parameters:
* `--sample-readlength-from-fastq` makes it copy the read length distribution from the training reads.
* `--prob-sub 0.0133` is the substitution rate in the training reads.
* `--prob-ins 0.00626` is the insertion rate in the training reads.
* `--prob-del 0.0136` is the deletion rate in the training reads.
* `--no-sam` because I'm just interested in the simulated reads, not their alignment.
* `--norm-params -0.1 0.8` increases the accuracy variance so the reads aren't too clustered around Q15.
* `--max-passes 1` is to turn off PacBio-style CCS.

```bash
cd ~/2025-08_ONT_read_simulator_benchmark
mkdir simlord
cd simlord

conda activate simlord

simlord --read-reference ../reference.fasta --coverage 100 --sample-readlength-from-fastq ../real_reads/training_reads.fastq --prob-sub 0.0133 --prob-ins 0.00626 --prob-del 0.0136 --max-passes 1 --norm-params -0.1 0.8 --no-sam simlord
```





## Analysis

Gather all read sets in one directory:
```bash
cd ~/2025-08_ONT_read_simulator_benchmark
mkdir analysis
cp real_reads/test_reads.fastq analysis/real.fastq
cp nanosim/nanosim.fastq analysis/
cp longislnd/longislnd.fastq analysis/
cp pbsim3/pbsim3.fastq analysis/
cp badread/badread.fastq analysis/
cp lrsim/lrsim.fastq analysis/
cp simlord/simlord.fastq analysis/
```

Stats on the full read sets:
```bash
cd ~/2025-08_ONT_read_simulator_benchmark/analysis
seqkit stats -a *.fastq
```

Per-read analysis:
```bash
conda activate mapping

cd ~/2025-08_ONT_read_simulator_benchmark/analysis
for r in badread longislnd lrsim nanosim pbsim3 real simlord; do
    minimap2 -t 32 -c --eqx ../reference.fasta "$r".fastq | grep "tp:A:P" > "$r".paf
    python3 ../scripts/extract_read_stats.py "$r".fastq "$r".paf > "$r".tsv
done

head -n1 real.tsv > read_analysis_results.tsv
tail -n+2 badread.tsv longislnd.tsv lrsim.tsv nanosim.tsv pbsim3.tsv real.tsv simlord.tsv >> read_analysis_results.tsv
```

Per-base analysis:
```bash
cd ~/2025-08_ONT_read_simulator_benchmark/analysis
for r in badread longislnd lrsim nanosim pbsim3 real simlord; do
    python3 ../scripts/per_base_qscores.py "$r".fastq ../reference.fasta "$r".paf > "$r".qscore_counts.tsv
done
```





# Error analysis

Install Pomoxis:
```bash
conda create -n pomoxis pomoxis pandas=2.2.3
```

Align reads (BAM format):
```bash
conda activate mapping

cd ~/2025-08_ONT_read_simulator_benchmark/analysis
for r in badread longislnd lrsim nanosim pbsim3 real simlord; do
    minimap2 -t 32 -a -x map-ont --eqx --MD ../reference.fasta "$r".fastq | samtools view -u -F 0x904 | samtools sort > "$r".bam
    samtools index "$r".bam
done
```

Homopolymers:
```bash
conda activate pomoxis

cd ~/2025-08_ONT_read_simulator_benchmark/analysis
for r in badread longislnd lrsim nanosim pbsim3 real simlord; do
    assess_homopolymers count -o "$r".homopolymers -t 32 "$r".bam
done

{
    printf "tool\t"
    cut -f1,14,15 real.homopolymers/hp_correct_vs_len.txt | head -n1
    for r in badread longislnd lrsim nanosim pbsim3 real simlord; do
        awk -F'\t' -v OFS='\t' -v r="$r" 'NR > 1 {print r, $1, $14, $15}' \
            "$r".homopolymers/hp_correct_vs_len.txt
    done
} > homopolymers.tsv
```

3-mer substitutions:
```bash
cd ~/2025-08_ONT_read_simulator_benchmark/analysis
for r in real badread longislnd lrsim nanosim pbsim3 simlord; do
    python3 ../scripts/count_3mer_subs.py "$r".fastq ../reference.fasta "$r".paf > "$r".3mer_subs
done
```
