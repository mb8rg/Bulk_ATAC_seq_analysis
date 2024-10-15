# Bulk ATAC-seq library preparation 

## Cell line list (3 biological replicates each)
1. COLO858 WT
2. COLO858 Cas9 Control
3. COLO858 Cas9 FOS KO-B8S
4. COLO858 Cas9  FOSL1 KO-E2S
5. COLO858 Cas9 FOSL2 KO-A2S
6. COLO858 Cas9 / Positive ORF Control
7. COLO858 Cas9 FOSL1 KO-E2S / FOSL1 OE
8. COLO858 Positive ORF Control
9. COLO858 FOS OE
10. COLO858 FOSL2 OE
11. COLO858 JUN OE
12. COLO858 NTC shRNA #2 Control
13. COLO858 JUN shRNA #3
14. COLO858 JUND shRNA #1
15. COLO858 JUNB shRNA #1
16. COLO858 Cas9 / NTC shRNA #2 Control
17. COLO858 Cas9 FOSL2 KO-A2S / JUN shRNA #3
18. COLO858 Cas9 FOS KO-B8S / JUND shRNA #1
19. COLO858 Cas9 FOS KO-B8S / JUNB shRNA #1
20. COLO858 Cas9 JUN KO-A4 / JUND shRNA #1

## Library preparation protocol






# Bulk ATAC-seq analysis
## Required Software
- FastQC
  - verson 0.11.5
  - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Trim Galore 
  - version 0.6.4
  - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Cutadapt
  - version 3.4
- gcc
  - version 13.3.0
- Bowtie2
  - version 2.5.1
- QualiMap
   - version 2.2.1
- Samtools
  - version 1.17
  - https://www.htslib.org/
- python
  - version 3.11.4
- java
  - version 21
- PICARD
  - version 2.27.5
- MACS2
  - version 2.2.7.1
- Bedtools
  - version 2.30.0
    
## 1. Sequencing quality check
FASTQC (version 0.11.5) was used to determine the quality of sequencing and to check for presence of adapters. I run the script using the following command: `sbatch --array=0-119 01_fastqc_array.slurm`. The script takes each fastq file as an input and outputs html documents with information about the quality of reads.

The following is part of the "01_fastqc_array.slurm" file:
```bash
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
```
The quality scores across bases are good and there is some adapter content matching Nextera Transposase Sequence towards the ends of the reads.

## 2. Trimming of adapters and low quality reads
Adapter trimming was performed using Trim Galore (version 0.6.4).  I run "02_adapter_trimming_array.slurm" with the following command: `sbatch --array=0-59 02_adapter_trimming_array.slurm`. Note: the array range is from 0 to 59 because each array number corresponds to each PAIR of forward and reverse FASTQ files. 

Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) is a wrapper script that utilizes Cutadapt to trim reads and FastQC to check the quality of reads after trimming. The following is part of the "02_adapter_trimming_array.slurm" with explainations for some of the arugements used. 

The following is part of the "02_adapter_trimming_array.slurm" file:
```bash
echo "Processing file $FILE_FOR and $FILE_REV"
start=$SECONDS
module purge
module load trimgalore/0.6.4

trim_galore --cores 4 \ # uses 4 cores
 --nextseq 20 \ # 2-color high quality G-trimming enabled, with quality cutoff of 20. This argument is used here because the samples were sequenced using Illumina NextSeq, which uses 2-color chemistry. It interprets lack of signal as "G" base. 
 --nextera \ # The Tn5 transposase was loaded with Nextra adapters. The trimming will be done specifically for these adapters.
 --length 20 \ # Remove any reads that are less than 20 nt long after trimming
 --paired \ # Specifies that paireds-end sequencing was performed and the trimming is done for each pair at a time
 --fastqc_args "--outdir $FASTQC_OUTDIR" \
 --output_dir $OUTDIR/ $FILE_FOR $FILE_REV

end=$SECONDS
echo "duration:$((end-start)) seconds."
```
The FASTQC files indicated that the trimming successfully removed adapters.

## 3. Alignment to human genome
Bowtie2 (version 2.5.1) was used for alignement. The reference genome index was build using GRCh38 primary assembly genome fasta and annotation gtf files from gencode release version 46. The alignment requires two steps: (1) generation of genome indexes and (2) read alignment.

### 3a. Generating Genome Indexes
I used Human CRCh38 primary assembly genome fasta and annotation from gencode release version 46. I made folders named "hg38_gencode_v46_fasta" and "hg38_gencode_v46_gtf" and downloaded the FASTA and GTF files using the following commands:  
`wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz`  
`wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz`  
Both files were decompressed using `gunzip`.  

To generate the genome index for alignmnet, I run "03a_bowtie_build_ref_genome_index.slurm" with the following command: `sbatch 03a_bowtie_build_ref_genome_index.slurm`. The indexing of the genome is only done once and the same index will be used for all of the samples.

The following is part of the "03a_bowtie_build_ref_genome_index.slurm" file:
```bash
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
```
## 3b. Read alignment 
Each pair of trimmed reads is aligned to the hg38 reference genome index using Bowtie2 aligner. I run the "03b_bowtie_alignment_array.slurm" using the following command: `sbatch --array=0-59 03b_bowtie_alignment_array.slurm`. The script uses very-sensitive mode and end-to-end read alignment, which searches for alignment using the full reads. The `--no-mixed` option means that the aligner will only look for paired-end alignments for each pair of reads. `-X 2000` argument tells the aligner to looks for reads up to 2000 bp long. I set the upper limit very high because I want to obtain the fragment length distribution in step 6 in order to determine is the samples have both nucleosome-free and mono-nucleosome peaks. 

The following is part of the "03b_bowtie_alignment_array.slurm" file:
```bash
echo "Processing file $FILE_FOR and $FILE_REV"
echo "Sample ID: $SAMPLE_ID"
start=$SECONDS
module purge
module load bowtie2/2.5.1

# Align reads to human hg38 reference genome index
bowtie2 \
 -p 8 \
 --end-to-end \ # end-to-end read alignment
 --very-sensitive \
 --no-mixed --phred33 -X 2000 \ 
 -x $HG38_GENOME_INDEX_PATH/hg38_gencode_v46_genomeindex \
 -1 $FILE_FOR -2 $FILE_REV \
 -S $OUTPUT_PATH/${SAMPLE_ID}_bowtie2.sam \
 &> $OUTPUT_PATH/bowtie2_summary/${SAMPLE_ID}_bowtie2.txt

end=$SECONDS
echo "duration:$((end-start)) seconds."
```

## 4. Checking the quality of the read alignmnet
I checked the quality of read alignment using QualiMap using bamqc argument. I run the "04_qualimap_array.slurm" using the following command: `sbatch --array=0-59 04_qualimap_array.slurm`.

The following is part of the "04_qualimap_array.slurm" file:
```bash
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
```

## 5. Filtering reads
Before peak calling, I performed filtering using "05_filtering_reads.slurm" script, which performs the following:  
(1) Remove mitochondrial reads and sorting by coordinate using samtools. The mitochondrial reads are excluded using `grep -v chrM`.  
(2) The PCR and optical duplicates are flaged using PICARD.   
(3) Samtools is used to filter out reads that have flags matching the following: duplicates, low quality reads (q 30), reads unmapped, mate unmapped, not primary alignment, reads failing platform.  
(4) The final filtered bam file is indexed using samtools. The indexed file ending in `.bai` is necessary for viewing the bam file in IGV.   

The following is part of the "05_filtering_reads.slurm" file:
```bash
module purge
module load gcc samtools java picard

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
```

## 6. Peak calling
Peak calling was performed on the filtered bam files from step 5 using MACS2 callpeak function with the `--format BAMPE` argument, which tells MACS2 to pileup the real fragment length instead of an estimate. The `-g hs` indicates that human genome should be used for the effective genome size. `-B` argument indicates to save the  save  extended  fragment pileup, and local lambda tracks at every bp into a bedGraph file and `-SPMR` tells MACS2 to normalize these pileip files to signal per million  reads. 


The following is part of the "06_peak_calling.slurm" file:
``` bash
# Create directory and set variables
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
```

## 7. Fragment length 
I used the "07_fragment_length.slurm" file to extract fragment lengths using samtools view function. To plot the fragment length distribution, I imported the txt files to R and made fragment length vs read count line plots using ggplot2. This step is optional. Technically the QualiMap files have a fragment length distribution plot at the end of each html file. I used this script to compile all of the fragment length distribution plots into one figure. 

The following is part of the "07_fragment_length.slurm" file:
```bash
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
```
![Fragment_length](https://github.com/user-attachments/assets/9de13f84-378e-4099-ad11-360d79e1e927)

## 8. Extracting metrics
I used the "08_extract_metrics.slurm" script to extract the following metrics:  
(1) Total number of reads  
(2) Number of reads that map to mitochondrial DNA   
(3) Number of duplicates   
(4) Total number of reads after filtering for: duplicates, reads with q less than 30, reads unmapped, mate unmapped, not primary alignment, read fails platform/vendor quality check  
(5) Total number of peaks called   
(6) FRiP score   

The following is part of the "08_extract_metrics.slurm" file:
```bash
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
```






