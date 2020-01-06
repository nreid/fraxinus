#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=0-10

hostname
date

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# load software
module load fastqc


#input/output directories, supplementary files

INDIR=/labs/Wegrzyn/EAB_github/Data
FASTQS=($(ls -1 $INDIR/Sample*R1*fastq.gz))

OUTDIR=../../results/fastqc
mkdir -p $OUTDIR

FQ1=${FASTQS[$SLURM_ARRAY_TASK_ID]}
FQ2=$(echo $FQ1 | sed 's/_R1_/_R2_/')

# run fastqc. "*fq" tells it to run on all fastq files in directory "../rawdata/"
fastqc -t 6 -o $OUTDIR $FQ1 $FQ2