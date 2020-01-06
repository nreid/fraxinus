#!/bin/bash
#SBATCH --job-name=process_radtags
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-5]

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

module load stacks/2.41

#input/output directories, supplementary files
INDIR=/labs/Wegrzyn/EAB_github/Data

POOLS=($(ls -1 $INDIR/ddRAD*R1*fastq.gz))
BARCODES=($(ls -1 ../../metadata/*barcode*))

# make demultiplexed directory if it doesn't exist
OUTDIR=../../results/demultiplexed_fastqs/pool_$SLURM_ARRAY_TASK_ID
mkdir -p $OUTDIR

FASTQ1=$(echo ${POOLS[$SLURM_ARRAY_TASK_ID]})
FASTQ2=$(echo $FASTQ1 | sed 's/R1/R2/')
BC=$(echo ${BARCODES[$SLURM_ARRAY_TASK_ID]})

echo demultiplexing file pair $FASTQ1 using barcode set $BC

process_radtags \
-1 $FASTQ1 \
-2 $FASTQ2 \
-b $BC \
-o $OUTDIR \
-i gzfastq \
-y gzfastq \
-e pstI \
-c \
-q \
-s 20 \



