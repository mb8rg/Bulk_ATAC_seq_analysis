#!/bin/bash
#SBATCH -A rivanna_fallahi_lab		# account
#SBATCH -p standard	# partition/queue
#SBATCH --nodes=1		# number of compute nodes
#SBATCH --ntasks=1		# number of program instances
#SBATCH --time=04:00:00		# max time before job cancels
#SBATCH --mem-per-cpu=16GB               # memory
#SBATCH --cpus-per-task=8

# Create directory and set variables
PROJECT_PATH=/scratch/mb8rg/20240814_AP1_perturbations_Bulk_ATACseq
DATA_PATH=$PROJECT_PATH/Trimmed_reads
HG38_GENOME_INDEX_PATH=$PROJECT_PATH/Hg38_gencode_v46_genomeindex
OUTPUT_PATH=$PROJECT_PATH/Alignment/sam
mkdir -p $OUTPUT_PATH/bowtie2_summary

# Get list of all data files
FILES_FOR=($(ls -1 $DATA_PATH/*_R1_001_val_1.fq.gz))
FILES_REV=($(ls -1 $DATA_PATH/*_R2_001_val_2.fq.gz))

# Use Slurm Array number to select file for this job
FILE_FOR=${FILES_FOR[$SLURM_ARRAY_TASK_ID]}
FILE_REV=${FILES_REV[$SLURM_ARRAY_TASK_ID]}
SAMPLE_ID=($(basename ${FILE_FOR%%_L001_R*}))

echo "Processing file $FILE_FOR and $FILE_REV"
echo "Sample ID: $SAMPLE_ID"
start=$SECONDS
module purge
module load bowtie2/2.5.1

# Align reads to human hg38 reference genome index
bowtie2 \
 -p 8 --end-to-end --very-sensitive \
 --no-mixed --phred33 -X 2000 \
 -x $HG38_GENOME_INDEX_PATH/hg38_gencode_v46_genomeindex \
 -1 $FILE_FOR -2 $FILE_REV \
 -S $OUTPUT_PATH/${SAMPLE_ID}_bowtie2.sam \
 &> $OUTPUT_PATH/bowtie2_summary/${SAMPLE_ID}_bowtie2.txt

end=$SECONDS
echo "duration:$((end-start)) seconds."
