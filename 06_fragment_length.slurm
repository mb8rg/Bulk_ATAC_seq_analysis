#!/bin/bash
#SBATCH -A rivanna_fallahi_lab		# account
#SBATCH -p standard	# partition/queue
#SBATCH --nodes=1		# number of compute nodes
#SBATCH --ntasks=1		# number of program instances
#SBATCH --time=00:30:00		# max time before job cancels
#SBATCH --mem=64GB               # memory
#SBATCH --cpus-per-task=4

# Create directory and set variables
PROJECT_PATH=/scratch/mb8rg/20240814_AP1_perturbations_Bulk_ATACseq
DATA_PATH=$PROJECT_PATH/Alignment_maxins1000/sam
OUTPUT_PATH=$DATA_PATH/Fragment_length
mkdir -p $OUTPUT_PATH

# Get list of all sam data files
FILES=($(ls -1 $DATA_PATH/*_bowtie2.sam))

# Use Slurm Array number to select file for this job
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE_ID=($(basename ${FILE%%_bowtie2.sam}))

echo "Processing file $FILE"
echo "Sample ID: $SAMPLE_ID"
start=$SECONDS

module purge
module load gcc
module load samtools

# Extract the 9th column from the alignment sam file which is the fragment length
samtools view -@ 4 -F 0x04 $FILE | \
awk -F'\t' 'function abs(x){return ((x<0.0) ? -x : x)} {print abs($9)}' | \
sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' \
> $OUTPUT_PATH/${SAMPLE_ID}_fragmentLen.txt

end=$SECONDS
echo "duration:$((end-start)) seconds."
