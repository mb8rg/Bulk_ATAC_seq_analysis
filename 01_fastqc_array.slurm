#!/bin/bash
#SBATCH -A rivanna_fallahi_lab		# account
#SBATCH -p standard	# partition/queue
#SBATCH --nodes=1		# number of compute nodes
#SBATCH --ntasks=1		# number of program instances
#SBATCH --time=00:10:00		# max time before job cancels
#SBATCH --mem=6GB               # memory

# Create directory and set variables
PROJECT_PATH=/scratch/mb8rg/20240814_AP1_perturbations_Bulk_ATACseq/
DATA_PATH=$PROJECT_PATH/Original_FASTQ_files_from_miSeq
OUTDIR=$PROJECT_PATH/FASTQC_pretrim_output
mkdir -p $OUTDIR

# Get list of all data files
FILES=($(ls -1 $DATA_PATH/*.fastq.gz))

# Use Slurm Array number to select file for this job
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Processing file $FILE"
module purge
module load fastqc/0.11.5

fastqc $FILE --outdir $OUTDIR

