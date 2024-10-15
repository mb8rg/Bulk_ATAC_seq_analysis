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
PEAKS_PATH=$PROJECT_PATH/MACS2_peak_calling

# Get list of all sam data files
FILES=($(ls -1 $DATA_PATH/sam/*_bowtie2.sam))

# Use Slurm Array number to select file for this job
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE_ID=($(basename ${FILE%%_bowtie2.sam}))

echo "Sample_ID $SAMPLE_ID"

module purge
module load samtools bedtools

# Total number of reads
total_reads=`samtools view -@ 8 $FILE -O sam | wc -l`
echo "Total_reads $total_reads"

# Number of reads that map to mitochondrial DNA
rmChrM=`samtools view -@ 8 $DATA_PATH/bam/${SAMPLE_ID}.rmChrM.bam -O sam | wc -l`
echo "Mitochondrial_reads $((total_reads-rmChrM))"

# Duplicates
dups=`samtools view -@ 8 -f 1024 $DATA_PATH/bam/${SAMPLE_ID}.dup_marked.bam -O sam | wc -l`
echo "Total_duplicates $dups"

# Total number of reads after filtering for: duplicates, reads with q less than 30,
# reads unmapped, mate unmapped, not primary alignment, read fails platform/vendor quality check
filtered=`samtools view -@ 8 $DATA_PATH/bam/${SAMPLE_ID}.filtered.bam -O sam | wc -l`
echo "Reads_after_all_filtering $filtered"

# Total number of peaks called
peaks=`wc -l $PEAKS_PATH/$SAMPLE_ID/${SAMPLE_ID}_peaks.narrowPeak`
echo "Total_peaks ${peaks% *}"

# Calculate reads in peaks
reads_in_peaks=`bedtools sort -i $PEAKS_PATH/$SAMPLE_ID/${SAMPLE_ID}_peaks.narrowPeak \
  | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
  -a $DATA_PATH/bam/${SAMPLE_ID}.filtered.bam -b stdin -ubam | samtools view -O sam | wc -l`
echo "Reads_in_peaks $reads_in_peaks"

# FRiP score
FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${filtered}"}")
echo "FRiP $FRiP"
