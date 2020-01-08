#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=0-11


hostname
date

# load software
module load bwa/0.7.17
module load samtools/1.9
module load samblaster/0.1.24

# input, output files and directories
INDIR=../../results/trimmed_WGS

OUTDIR=../../results/aligned_WGS
mkdir -p $OUTDIR

# indexed reference genome
REFERENCE=/labs/Wegrzyn/EAB_github/Genome_Stats/remove_2kb/fraxinus_pennsylvanica_21Jun2018_formatted.fa-filtered

# fastq array
FASTQS=($(find $INDIR -name "*1P.fastq.gz"))

# fastq and bam files for this task
FQ1=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]} | sed 's/.*\///')
FQ2=$(echo $FQ1 | sed 's/1P/2P/')
OUTFILE=$(echo $FQ1 | sed 's/_1P.fastq.gz/.bam/')

# sample ID, read group
SAM=$(echo $OUTFILE | sed 's/\..*//')
RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)

# alignment pipe
bwa mem -t 4 -R $RG $REFERENCE $INDIR/$FQ1 $INDIR/$FQ2 | \
samblaster | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$SAM - >$OUTDIR/$OUTFILE

samtools index $OUTDIR/$OUTFILE

date
