#!/bin/bash
#SBATCH -A rivanna_fallahi_lab		# account
#SBATCH -p standard	# partition/queue
#SBATCH --nodes=1		# number of compute nodes
#SBATCH --ntasks=1		# number of program instances
#SBATCH --time=05:00:00		# max time before job cancels
#SBATCH --mem=120GB               # memory
#SBATCH --cpus-per-task=8

# Create directory and set variables
PROJECT_PATH=/scratch/mb8rg/20240814_AP1_perturbations_Bulk_ATACseq
HG38_GENOME_INDEX_PATH=$PROJECT_PATH/Hg38_gencode_v46_genomeindex
mkdir -p $HG38_GENOME_INDEX_PATH

start=$SECONDS
module purge
module load gcc
module load bowtie2/2.5.1

# Build the bowtie2 reference human hg38 genome index
bowtie2-build --threads 8 \
$PROJECT_PATH/Hg38_gencode_v46_fasta/GRCh38.primary_assembly.genome.fa \
$HG38_GENOME_INDEX_PATH/hg38_gencode_v46_genomeindex

end=$SECONDS
echo "duration:$((end-start)) seconds."
