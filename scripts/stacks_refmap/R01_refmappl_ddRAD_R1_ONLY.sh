#!/bin/bash 
#SBATCH --job-name=refmap.pl
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=30G
#SBATCH --qos=general
#SBATCH --partition=general


hostname
date

# load software
module load stacks/2.41

############

# this script generates preliminary results from R1 ddRAD data ONLY

############

# input, output files, directories

INDIR=../../results/align_R1_ddRAD
OUTDIR=../../results/stacks/refmap
mkdir -p $OUTDIR

# popmap file

POPMAP=../../metadata/R1_popmap.txt
ls ${INDIR}/*bam | grep -oP "FP...." | sed 's/$/\tnopop/' >$POPMAP

# refmap.pl -s option is broken. 
ref_map.pl \
--samples $INDIR \
--popmap $POPMAP \
-o $OUTDIR \
-T 10 \
-X "populations:--vcf" \
-X "populations:--fasta-samples" \
-X "populations:--fasta-loci"
