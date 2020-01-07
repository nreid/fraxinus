#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=0-10

module load java
module load Trimmomatic/0.36

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

#input/output directories, supplementary files

ADAPT=../../metadata/wgs_adapters.fa

INDIR=/labs/Wegrzyn/EAB_github/Data
FASTQS=($(find $INDIR/ -name "Sample*R1*fastq.gz"))
INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]} | sed 's/.*\///')
OUTFILE=$(echo $INFILE | sed 's/_S._R1_.*/.fastq.gz/')

# make trimmed directory if it doesn't exist
OUTDIR=../../results/trimmed_WGS
mkdir -p $OUTDIR

java -jar $Trimmomatic PE \
-threads 10 \
-basein $INDIR/$INFILE \
-baseout $OUTDIR/$OUTFILE \
ILLUMINACLIP:$ADAPT:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:20 \
MINLEN:45

