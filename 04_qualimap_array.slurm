#!/bin/bash
#SBATCH -A rivanna_fallahi_lab		# account
#SBATCH -p standard	# partition/queue
#SBATCH --nodes=1		# number of compute nodes
#SBATCH --ntasks=1		# number of program instances
#SBATCH --time=01:00:00		# max time before job cancels
#SBATCH --mem-per-cpu=9GB               # memory
#SBATCH -c 4

# Create directory and set variables
PROJECT_PATH=/scratch/mb8rg/20240814_AP1_perturbations_Bulk_ATACseq
DATA_PATH=$PROJECT_PATH/Alignment/sam
BAM_FILE_DATA_PATH=$PROJECT_PATH/Alignment/bam
mkdir -p $BAM_FILE_DATA_PATH

# Get list of all data files
SAM_FILES=($(ls -1 $DATA_PATH/*_bowtie2.sam))
GTF_FILE=$PROJECT_PATH/Hg38_gencode_v46_gtf/gencode.v46.primary_assembly.annotation.gtf

# Use Slurm Array number to select file for this job
SAM_FILE=${SAM_FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE_ID=($(basename ${SAM_FILE%%_bowtie2.sam}))

# Create output directory for each sample
OUTPUT_PATH=$PROJECT_PATH/Qualimap_results_filtered_reads_only/$SAMPLE_ID
mkdir -p $OUTPUT_PATH

echo "Processing file $SAM_FILE"
echo "Sample ID: $SAMPLE_ID"
start=$SECONDS
module purge
unset DISPLAY
module load gcc/13.3.0 samtools/1.17 qualimap/2.2.1

# Generate a BAM file sorted by location
samtools sort -@ 4 $SAM_FILE -o $BAM_FILE_DATA_PATH/$SAMPLE_ID.sorted.bam

# Run qualimap on each file in the array
qualimap bamqc -nt 4 \
 -bam $BAM_FILE_DATA_PATH/$SAMPLE_ID.sorted.bam \
 --feature-file $GTF_FILE \
 --paint-chromosome-limits \
 --genome-gc-distr HUMAN \
 -outdir $OUTPUT_PATH \
 --java-mem-size=32G

end=$SECONDS
echo "duration:$((end-start)) seconds."
