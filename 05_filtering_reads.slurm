#!/bin/bash
#SBATCH -A rivanna_fallahi_lab		# account
#SBATCH -p standard	# partition/queue
#SBATCH --nodes=1		# number of compute nodes
#SBATCH --ntasks=1		# number of program instances
#SBATCH --time=02:00:00		# max time before job cancels
#SBATCH --mem=64GB               # memory
#SBATCH --cpus-per-task=8

# Create directory and set variables
PROJECT_PATH=/scratch/mb8rg/20240814_AP1_perturbations_Bulk_ATACseq
DATA_PATH=$PROJECT_PATH/Alignment/
TEMP_PATH=$DATA_PATH/Temp_path
mkdir -p $DATA_PATH/Picard_metrics
mkdir -p $DATA_PATH/bam
mkdir -p $TEMP_PATH

# Get list of all sam data files
FILES=($(ls -1 $DATA_PATH/sam/*_bowtie2.sam))

# Use Slurm Array number to select file for this job
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE_ID=($(basename ${FILE%%_bowtie2.sam}))

echo "Processing file $FILE"
echo "Sample ID: $SAMPLE_ID"
start=$SECONDS

module purge
module load gcc
module load samtools
module load java
module load picard

# Remove mitochondrial reads and sort by coordinate
samtools view -@ 8 -h $FILE | grep -v chrM | samtools sort -@ 8 -o $DATA_PATH/bam/$SAMPLE_ID.rmChrM.bam
echo "Finished removing mitochondrial reads and sorting"

# Mark duplicates with a flag
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
 --REMOVE_DUPLICATES false \
 --INPUT $DATA_PATH/bam/${SAMPLE_ID}.rmChrM.bam \
 --OUTPUT $DATA_PATH/bam/${SAMPLE_ID}.dup_marked.bam \
 --CREATE_INDEX true \
 --VALIDATION_STRINGENCY LENIENT \
 --TMP_DIR $TEMP_PATH \
 --METRICS_FILE $DATA_PATH/Picard_metrics/${SAMPLE_ID}_picard.rmDup.txt
echo "Finished marking duplicates"

# Remove duplicates, low quality reads (q 30), reads unmapped,
# mate unmapped, not primary alignment, reads failing platform and
# duplicates
samtools view -h -b -F 1804 -f 2 $DATA_PATH/bam/${SAMPLE_ID}.dup_marked.bam > $DATA_PATH/bam/${SAMPLE_ID}.filtered.bam
samtools index -@ 8 $DATA_PATH/bam/${SAMPLE_ID}.filtered.bam
echo "Finished filtering low quality reads and duplicates"

end=$SECONDS
echo "duration:$((end-start)) seconds."
