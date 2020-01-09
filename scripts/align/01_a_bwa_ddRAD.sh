#!/bin/bash
#SBATCH --job-name=bwa
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


# load software
module load bwa/0.7.17
module load samtools/1.9
module load umi_tools/1.0.0

SLURM_ARRAY_TASK_ID=0

#####
# Define a bunch of variables for this array task
#####

# input, output directories
INDIR=../../results/demultiplexed_fastqs
OUTDIR=../../results/aligned_tempddRAD
mkdir -p $OUTDIR
DEDUPDIR=../../results/aligned_ddRAD
mkdir -p $DEDUPDIR

# fastq array
FASTQS=($(find $INDIR -name "*1.fq.gz" | grep -v "\.rem\."))

# input, output files for this array task
FQ1=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]})
FQ2=$(echo $FQ1 | sed 's/1.fq.gz/2.fq.gz/')
BAM=$(echo $FQ1 | sed 's/.1.fq.gz/.bam/' | grep -oP "FP.*")

# named pipes for this task to redirect umi-tagged fastqs
P1=$OUTDIR/$(echo $FQ1 | grep -oP "FP......")
P2=$OUTDIR/$(echo $FQ2 | grep -oP "FP......")

mkfifo $P1
mkfifo $P2

# reference genome
REFERENCE=/labs/Wegrzyn/EAB_github/Genome_Stats/remove_2kb/fraxinus_pennsylvanica_21Jun2018_formatted.fa-filtered
# sample ID, read group
SAM=$(echo $BAM | sed 's/.bam//')
RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)

######
# Execute umi_tools and bwa
# uses named pipes to avoid writing intermediate fastq files
######

# run umi tools
# two problems with this tool
	# bc-pattern - there is no umi on R1, but you have to specify it anyway. X means don't cut it off. 
	# /1 and /2 tags added by stacks to read names cause an error. removing them with process substitution. 

umi_tools extract \
--bc-pattern=XXXX \
--bc-pattern2=NNNNNNNNNN \
-I <(zcat $FQ1 | sed 's/..$//;n;n;n') \
-S >(cat >$P1) \
--read2-in=<(zcat $FQ2 | sed 's/..$//;n;n;n') \
--read2-out=>(cat >$P2) &

# run bwa, sort, index

bwa mem -t 4 -R $RG $REFERENCE $P1 $P2 | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$SAM - >$OUTDIR/$BAM

samtools index $OUTDIR/$BAM
 
rm $P1
rm $P2

umi_tools dedup -I $OUTDIR/$BAM --paired -S $DEDUPDIR/$BAM
samtools index $DEDUPDIR/$BAM


