#!/bin/bash
#SBATCH -A rivanna_fallahi_lab		# account
#SBATCH -p standard	# partition/queue
#SBATCH --nodes=1		# number of compute nodes
#SBATCH --ntasks=1		# number of program instances
#SBATCH --time=00:40:00		# max time before job cancels
#SBATCH --mem-per-cpu=9GB               # memory
#SBATCH -c 4

# Create directory and set variables
PROJECT_PATH=/scratch/mb8rg/20240814_AP1_perturbations_Bulk_ATACseq
DATA_PATH=$PROJECT_PATH/Original_FASTQ_files_from_NextSeq_Run_1
OUTDIR=$PROJECT_PATH/Trimmed_reads
FASTQC_OUTDIR=$PROJECT_PATH/FASTQC_aftertrim_output
mkdir -p $OUTDIR
mkdir -p $FASTQC_OUTDIR

# Get list of all data files
FILES_FOR=($(ls -1 $DATA_PATH/*_R1_001.fastq.gz))
FILES_REV=($(ls -1 $DATA_PATH/*_R2_001.fastq.gz))

# Use Slurm Array number to select file for this job
FILE_FOR=${FILES_FOR[$SLURM_ARRAY_TASK_ID]}
FILE_REV=${FILES_REV[$SLURM_ARRAY_TASK_ID]}

echo "Processing file $FILE_FOR and $FILE_REV"
start=$SECONDS
module purge
module load trimgalore/0.6.4

trim_galore --cores 4 --nextseq 20 --nextera \
 --length 20 --paired --fastqc_args "--outdir $FASTQC_OUTDIR" \
 --output_dir $OUTDIR/ $FILE_FOR $FILE_REV

end=$SECONDS
echo "duration:$((end-start)) seconds."
