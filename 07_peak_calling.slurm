#!/bin/bash
#SBATCH -A rivanna_fallahi_lab		# account
#SBATCH -p standard	# partition/queue
#SBATCH --nodes=1		# number of compute nodes
#SBATCH --ntasks=1		# number of program instances
#SBATCH --time=02:00:00		# max time before job cancels
#SBATCH --mem=64GB               # memory

# Create directory and set variables
PROJECT_PATH=/scratch/mb8rg/20240814_AP1_perturbations_Bulk_ATACseq
DATA_PATH=$PROJECT_PATH/Alignment/
OUTPUT_PATH=$PROJECT_PATH/MACS2_peak_calling
mkdir -p $OUTPUT_PATH

# Get list of all sam data files
FILES=($(ls -1 $DATA_PATH/bam/*.filtered.bam))

# Use Slurm Array number to select file for this job
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE_ID=($(basename ${FILE%%.filtered.bam}))

echo "Processing file $FILE"
echo "Sample ID: $SAMPLE_ID"
start=$SECONDS

module purge
module load gcc python macs2

# Peak calling with MACS2
# --format BAMPE: MACS2 will pileup the real fragment length instead of an estimate
macs2 callpeak -t $FILE \
 --format BAMPE \
 -n ${SAMPLE_ID} \
 --outdir $OUTPUT_PATH/${SAMPLE_ID} \
-g hs -B --SPMR --call-summits --keep-dup all


end=$SECONDS
echo "duration:$((end-start)) seconds."
