#!/bin/bash
#SBATCH --job-name=bwa_R1
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-90]%20


######
# ALIGN ONLY DDRAD R1s TO GET
# PRELIMINARY DATA FOR MEETING
#####

hostname
date

# load software
module load bwa/0.7.17
module load samtools/1.9
module load samblaster/0.1.24

# input, output files and directories
INDIR=../../results/demultiplexed_fastqs

OUTDIR=../../results/align_R1_ddRAD
mkdir -p $OUTDIR

# indexed reference genome
REFERENCE=/labs/Wegrzyn/EAB_github/Genome_Stats/remove_2kb/fraxinus_pennsylvanica_21Jun2018_formatted.fa-filtered

# fastq array
FASTQS=($(find $INDIR -name "*1.fq.gz" | grep -v "\.rem\."))

# fastq and bam files for this task
FQ1=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]})
OUTFILE=$(echo $FQ1 | sed 's/.1.fq.gz/.bam/' | sed 's/.*\///')

# sample ID, read group
SAM=$(echo $OUTFILE | sed 's/\..*//')
RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)

# alignment pipe
bwa mem -t 4 -R $RG $REFERENCE $FQ1 | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$SAM - >$OUTDIR/$OUTFILE

samtools index $OUTDIR/$OUTFILE

date
