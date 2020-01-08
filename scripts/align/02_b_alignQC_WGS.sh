#!/bin/bash
#SBATCH --job-name=alignQC_WGS
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general


module load samtools/1.7
module load R/3.6.0

# This script does some QC on the alignments. 
# should add in the generation of coverage/targets files

INDIR=../../results/aligned_WGS

# make a list of bam files
find $INDIR -name "*bam" >$INDIR/bams.list

# make an output directory
OUTDIR=../../results/aligned_WGS_stats
mkdir -p $OUTDIR

# samtools bam statistics
for file in $(find $INDIR -name "*bam"); 
do samtools stats $file >${file}.stats
echo $file;
done

mv $INDIR/*stats $OUTDIR

# put the basic stats all in one file. 

FILES=($(find $OUTDIR -name "*bam.stats" | sort))

grep "^SN" ${FILES[0]} | cut -f 2 > $OUTDIR/SN.txt
for file in ${FILES[@]}
do paste $OUTDIR/SN.txt <(grep ^SN $file | cut -f 3) > $OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt
done

# add a header
find $OUTDIR -name "*bam.stats" | sort | sed 's/.bam.*//' | sed 's/.*\///' | tr "\n" "\t" | sed 's/\t$/\n/'>$OUTDIR/SN2.txt
# cat "\n" >>$OUTDIR/SN2.txt
cat $OUTDIR/SN.txt >>$OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt




