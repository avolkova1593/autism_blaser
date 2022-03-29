#!/bin/bash
#SBATCH --job-name=mouse_alignment # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=angelina.volkova@nyulangone.org # Where to send mail
#SBATCH --ntasks=4 # Run on a single CPU
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --ntasks=4
#SBATCH -p cpu_short

module load fastqc/0.11.7
module load trimgalore/0.5.0
module load python/cpu/2.7.15-ES ### CutAdapt is hidden in here module load bbmap/38.25
module load samtools/1.9
module load subread/1.6.3
module load bbmap/38.25

trim_galore \
--q 30 \
--phred33  \
-o /gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed  \
--fastqc /gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/fastq/${1}

bbmap.sh \
-Xmx26G \
ref=/gpfs/data/ruggleslab/databases/mm10/genome.fa \
in=/gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed/fastq/${1} \
outm=/gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed/fastq/${1}.sam \
minid=0.90 \
ambiguous=best \
nodisk

samtools view -S -b /gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed/fastq/${1}.sam > /gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed/fastq/${1}.bam 
samtools sort /gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed/fastq/${1}.bam -o /gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed/fastq/${1}_sorted.bam

featureCounts -s 2 -p -B -C -P \
--ignoreDup \
--primary \
-a /gpfs/data/ruggleslab/databases/mm10/mm10.RefSeq.whole.gene.gtf \
-g gene_id \
-o /gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed/${1}_feature_counts /gpfs/data/ruggleslab/microbiome/projects/autism/2016-07-06/trimmed/sorted/${1}

