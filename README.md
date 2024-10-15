[alignment_filtering_and_peak_calling_metrics.csv](https://github.com/user-attachments/files/17370545/alignment_filtering_and_peak_calling_metrics.csv)# Bulk ATAC-seq library preparation 

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




[Uploading alignment_filtering_and_Sample_ID,Total_reads,Mitochondrial_reads,Fraction of mitochondrial reads,Total duplicates after removing mitochondrial reads,Fraction of duplicates,Reads after all filtering,Fraction of reads after all filtering,Total_peaks,Reads_in_peaks,FRiP
S01_WT_rep1,"63,351,254","2,720,412",0.04,"27,496,014",0.45,"32,628,340",0.52,"78,435","12,503,472",0.38
S02_WT_rep2,"51,924,022","2,335,422",0.04,"23,961,560",0.48,"25,220,722",0.49,"70,918","8,495,751",0.34
S03_WT_rep3,"54,538,720","1,879,398",0.03,"25,849,558",0.49,"26,386,504",0.48,"76,565","10,391,221",0.39
S04_Cas9_Control_rep1,"87,322,432","4,087,800",0.05,"36,669,672",0.44,"45,841,556",0.52,"92,460","17,897,591",0.39
S05_Cas9_Control_rep2,"57,844,100","3,712,074",0.06,"25,756,934",0.48,"27,863,978",0.48,"75,156","9,686,504",0.35
S06_Cas9_Control_rep3,"62,201,682","3,727,686",0.06,"28,073,510",0.48,"29,895,706",0.48,"82,044","11,245,911",0.38
S07_FOS_KO_rep1,"76,989,020","3,014,414",0.04,"35,119,092",0.47,"38,173,192",0.5,"89,717","14,421,313",0.38
S08_FOS_KO_rep2,"42,543,836","1,820,156",0.04,"19,485,542",0.48,"20,892,096",0.49,"71,211","7,308,321",0.35
S09_FOS_KO_rep3,"44,088,420","2,203,482",0.05,"20,230,608",0.48,"21,292,422",0.48,"73,820","8,100,709",0.38
S10_FOSL1_KO_rep1,"43,548,952","559,106",0.01,"21,230,316",0.49,"21,391,618",0.49,"70,896","6,897,133",0.32
S11_FOSL1_KO_rep2,"45,486,460","652,900",0.01,"22,051,960",0.49,"22,370,218",0.49,"57,640","5,479,587",0.24
S12_FOSL1_KO_rep3,"47,036,462","781,120",0.02,"22,710,458",0.49,"23,152,740",0.49,"77,728","7,654,053",0.33
S13_FOSL2_KO_rep1,"49,929,148","1,731,484",0.03,"23,626,212",0.49,"24,114,256",0.48,"64,095","6,874,144",0.29
S14_FOSL2_KO_rep2,"49,546,774","2,002,592",0.04,"23,316,132",0.49,"23,775,910",0.48,"63,249","6,843,179",0.29
S15_FOSL2_KO_rep3,"51,496,330","1,574,846",0.03,"24,810,240",0.5,"24,647,622",0.48,"67,399","7,882,187",0.32
S16_Cas9_and_ORF_Control_rep1,"47,817,602","1,739,950",0.04,"22,156,688",0.48,"23,529,350",0.49,"68,902","8,824,212",0.38
S17_Cas9_and_ORF_Control_rep2,"49,562,684","2,134,654",0.04,"23,291,980",0.49,"23,694,228",0.48,"66,485","8,219,845",0.35
S18_Cas9_and_ORF_Control_rep3,"47,038,818","2,610,110",0.06,"21,890,138",0.49,"22,093,372",0.47,"60,572","6,889,255",0.31
S19_FOSL1_KO_OE_rep1,"45,249,486","1,467,868",0.03,"21,481,646",0.49,"21,934,744",0.48,"66,906","7,949,340",0.36
S20_FOSL1_KO_OE_rep2,"43,594,376","709,554",0.02,"20,704,626",0.48,"21,839,562",0.5,"67,257","8,005,549",0.37
S21_FOSL1_KO_OE_rep3,"46,656,442","1,062,230",0.02,"22,396,308",0.49,"22,844,970",0.49,"71,286","9,037,305",0.4
S22_ORF_Control_rep1,"44,934,126","1,226,734",0.03,"21,673,348",0.5,"21,535,106",0.48,"59,979","5,518,608",0.26
S23_ORF_Control_rep2,"43,516,870","1,277,762",0.03,"21,331,320",0.51,"20,519,496",0.47,"58,954","5,286,829",0.26
S24_ORF_Control_rep3,"44,899,932","1,599,790",0.04,"21,352,990",0.49,"21,549,806",0.48,"63,861","6,005,121",0.28
S25_FOS_OE_rep1,"34,964,682","687,488",0.02,"15,772,238",0.46,"18,134,650",0.52,"48,036","3,583,364",0.2
S26_FOS_OE_rep2,"47,404,598","1,344,548",0.03,"22,434,690",0.49,"23,159,582",0.49,"56,583","5,642,362",0.24
S27_FOS_OE_rep3,"49,472,754","1,440,354",0.03,"23,787,426",0.5,"23,785,916",0.48,"60,609","6,241,088",0.26
S28_FOSL2_OE_rep1,"43,102,758","1,566,272",0.04,"19,107,654",0.46,"22,023,582",0.51,"59,621","5,468,997",0.25
S29_FOSL2_OE_rep2,"48,851,482","1,882,834",0.04,"22,550,656",0.48,"23,941,232",0.49,"57,047","5,746,812",0.24
S30_FOSL2_OE_rep3,"47,850,038","1,564,138",0.03,"21,943,150",0.47,"23,876,228",0.5,"54,385","5,307,088",0.22
S31_JUN_OE_rep1,"47,652,516","1,322,812",0.03,"22,526,206",0.49,"23,321,630",0.49,"58,882","5,397,419",0.23
S32_JUN_OE_rep2,"42,184,768","1,009,610",0.02,"19,965,988",0.48,"20,796,702",0.49,"60,556","5,009,642",0.24
S33_JUN_OE_rep3,"41,296,444","903,382",0.02,"19,315,402",0.48,"20,635,604",0.5,"46,349","3,010,044",0.15
S34_NTC_KD_Control_rep1,"41,424,878","1,252,454",0.03,"18,756,276",0.47,"20,960,346",0.51,"70,025","7,103,129",0.34
S35_NTC_KD_Control_rep2,"44,466,646","1,815,688",0.04,"20,848,080",0.49,"21,309,126",0.48,"73,154","7,879,433",0.37
S36_NTC_KD_Control_rep3,"49,040,760","1,602,650",0.03,"21,419,812",0.45,"25,476,480",0.52,"72,138","8,634,791",0.34
S37_JUN_KD_rep1,"48,255,490","1,674,190",0.03,"21,822,102",0.47,"24,250,058",0.5,"64,201","7,649,299",0.32
S38_JUN_KD_rep2,"44,367,202","2,295,052",0.05,"20,373,224",0.48,"21,132,152",0.48,"65,554","6,963,995",0.33
S39_JUN_KD_rep3,"45,935,528","2,166,462",0.05,"19,729,096",0.45,"21,899,176",0.48,"62,106","7,512,712",0.34
S40_JUNB_KD_rep1,"44,133,918","966,838",0.02,"20,575,860",0.48,"22,117,072",0.5,"61,463","6,713,380",0.3
S41_JUNB_KD_rep2,"44,580,282","1,205,626",0.03,"21,094,674",0.49,"21,802,676",0.49,"63,744","7,507,671",0.34
S42_JUNB_KD_rep3,"47,247,252","1,135,344",0.02,"22,463,458",0.49,"23,127,406",0.49,"68,897","8,569,514",0.37
S43_JUND_KD_rep1,"46,231,250","2,686,770",0.06,"20,616,118",0.47,"22,442,634",0.49,"60,914","6,540,577",0.29
S44_JUND_KD_rep2,"48,137,500","1,901,196",0.04,"22,138,016",0.48,"23,569,536",0.49,"64,676","7,762,565",0.33
S45_JUND_KD_rep3,"46,478,564","1,551,584",0.03,"22,072,100",0.49,"22,364,360",0.48,"65,785","7,806,005",0.35
S46_Cas9_and_NTC_KD_Control_rep1,"45,830,088","2,665,514",0.06,"21,571,824",0.5,"21,207,956",0.46,"71,666","7,595,708",0.36
S47_Cas9_and_NTC_KD_Control_rep2,"45,800,882","2,468,148",0.05,"21,744,946",0.5,"21,185,992",0.46,"66,948","6,676,940",0.32
S48_Cas9_and_NTC_KD_Control_rep3,"39,420,732","1,646,512",0.04,"18,226,618",0.48,"19,218,116",0.49,"65,558","6,541,546",0.34
S49_FOSL2_KO_and_JUN_KD_rep1,"44,436,910","1,457,824",0.03,"20,188,468",0.47,"22,335,454",0.5,"64,959","7,110,773",0.32
S50_FOSL2_KO_and_JUN_KD_rep2,"42,178,588","1,076,024",0.03,"19,127,036",0.47,"21,574,434",0.51,"67,593","6,367,557",0.3
S51_FOSL2_KO_and_JUN_KD_rep3,"44,970,334","1,095,770",0.02,"21,879,466",0.5,"21,576,952",0.48,"67,656","6,253,812",0.29
S52_FOS_KO_and_JUNB_KD_rep1,"64,565,920","1,244,008",0.02,"25,046,970",0.4,"37,645,256",0.58,"84,289","12,641,989",0.34
S53_FOS_KO_and_JUNB_KD_rep2,"42,866,422","887,260",0.02,"20,247,748",0.48,"21,176,954",0.49,"65,722","6,138,518",0.29
S54_FOS_KO_and_JUNB_KD_rep3,"44,538,864","964,118",0.02,"21,242,856",0.49,"21,953,276",0.49,"73,286","7,450,463",0.34
S55_FOS_KO_and_JUND_KD_rep1,"41,425,528","1,327,938",0.03,"18,833,152",0.47,"20,891,386",0.5,"65,077","6,535,838",0.31
S56_FOS_KO_and_JUND_KD_rep2,"47,362,706","1,248,166",0.03,"22,780,026",0.49,"22,939,544",0.48,"72,312","7,420,160",0.32
S57_FOS_KO_and_JUND_KD_rep3,"43,333,428","1,788,914",0.04,"20,186,858",0.49,"20,999,366",0.48,"65,750","6,764,513",0.32
S58_JUN_KO_and_JUND_KD_rep1,"45,073,282","2,755,946",0.06,"20,922,588",0.49,"20,992,224",0.47,"59,962","6,368,858",0.3
S59_JUN_KO_and_JUND_KD_rep2,"44,503,318","2,420,928",0.05,"20,604,216",0.49,"21,069,266",0.47,"60,129","5,991,091",0.28
S60_JUN_KO_and_JUND_KD_rep3,"40,955,208","1,352,562",0.03,"18,988,082",0.48,"20,177,756",0.49,"43,704","3,051,269",0.15peak_calling_metrics.csv…]()








