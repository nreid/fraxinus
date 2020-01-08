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
#SBATCH --array=[0-90]%20

hostname
date

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# load software
module load fastqc


#input/output directories, supplementary files

INDIR=../../results/demultiplexed_fastqs
FASTQS=($(find $INDIR -name "*1.fq.gz" | grep -v "rem.1"))

OUTDIR=../../results/fastqc
mkdir -p $OUTDIR

FQ1=${FASTQS[$SLURM_ARRAY_TASK_ID]}
FQ2=$(echo $FQ1 | sed 's/1.fq.gz/2.fq.gz/')

# run fastqc. 
fastqc -t 6 -o $OUTDIR $FQ1 $FQ2