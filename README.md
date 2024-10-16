# lcWGS GREEN SEA TURTLE

---

# UPSTREAM ANALYSIS

```bash
RAW=/home/ralvarezv/raw-genomes  #Directory where raw data is located 
TRIM=/home/ralvarezv/new_workflow/trimmed_reads  #Directory where trimmed reads are located
MITOREF=/home/ralvarezv/reference  #Directory where the reference mitogenome is located 
SPLIT=/home/ralvarezv/new_workflow/bbsplit  #Directory where splitted reads are located
REF=/home/ralvarezv/new_workflow/reference #Directory where the reference genome is located
BAM=/home/ralvarezv/new_workflow/bam_files #Directory where .bam files (aligned) are located
METALIGN=/home/ralvarezv/new_workflow/metrics/bwa_align #Directory where metrics about alignment (using samtools-flagstat) are located
SORTBAM=/home/ralvarezv/new_workflow/sorted_bam		#Directory where sorted bam files are located.
```
# Step 1a: Quality Control of Raw reads - fastQC

FastQC performs a series of analyses, each represented by a module, to evaluate different aspects of data quality.
This tool takes as input one or more sequence files in various formats, such as FASTQ or BAM. In this case, we are going to use a fastq file.
- Input: fq.gz/fastq.gz files
- Output: .html and .zip files per sample (forward and reverse, separately)

Since we have numerous samples, we used a loop that will help us analyze all the samples with the following script:

```bash
#------Set relative path--------------------
RAW=/home/ralvarezv/raw_reads #Directory where raw reads are located
#------Command------------------------------
cat list.txt | while read a
do
  fastqc –o $RAW/${a} –t 1 –f fastq $RAW/${a}.gz
done
```
Details:
- The names of all the samples to be analyzed must be in the list.txt
- -o: the name that the output .html and .zip files will have
- -t: number of threads (cpu) to use to perform the analysis
- -f: to indicate that the input has fastq format, then you have provide the file path 

  
# Step 1b: Quality Control of Raw reads - multiQC
MultiQC is designed to aggregate and summarize results from multiple analysis tools, including FastQC, into a single, easy-to-read report.
To run multiqc we must be in the folder where our fastq files are located.
- Input: .fastqc files from all samples (forward and reverse, separately)
- Output: one .html file (all samples together)

Command:
``` bash
multiqc .
```
Details:
- By placing the '.' we are indicating that our input corresponds to all the fastqc files that are in our current directory

# Step 2: Trimming Raw reads - Trimmomatic
We used Trimmomatic to trim quality, remove adapter sequences, and filter out low-quality reads from raw sequencing data.
- Input: fq.gz files
- Output: four files per sample, two unpaired (U) and two paired (P)
- Add list of adapters (Illumina)

Script:
``` bash
# Loop through the list of samples to perform trimming
cat list.txt | while read sample
do
 java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred33 \
  $RAW/${sample}_R1.fq.gz $RAW/${sample}_R2.fq.gz \
  $TRIM/${sample}_R1_paired.fq $TRIM/${sample}_R1_unpaired.fq \
  $TRIM/${sample}_R2_paired.fq $TRIM/${sample}_R2_unpaired.fq \
  ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:75
done
```
Details:
- -phred33: the encoding scheme used for representing base quality scores in the input sequencing data files
- -baseout: the path and name of the output files
- ILLUMINACLIP: the path to the adapters files
- LEADING: any bases at the start of a read (5') with a quality score less than 3 will be trimmed off
- TRAILING: Bases at the end of a read (3') with a quality score less than 3 will be trimmed off
- MINLEN: reads with a length below 36 bases will be removed from the dataset

## Run fastQC and multiQC post-trimming using paired files to re-evaluate the content of adapters

# Step 3: Separation of the mitogenome from the nuclear genome | BBsplit
To separate the nuclear genome from the mitogenome we used BBSplit. 
This tool is designed to separate sequencing reads into different bins based on the reference databases provided by the user.

Script:
``` bash
# BBSplit command to separate mitochondrial reads from nuclear reads
cat list.txt | while read sample
do
  bbsplit.sh ref=$MITOREF/mitogenome_NC_000886.fasta \
  in1=$TRIM/${sample}_R1_paired.fq in2=$TRIM/${sample}_R2_paired.fq \
  basename=${sample}mt%.fq out1=$SPLIT/${sample}_nomt_R1.fq out2=$SPLIT/${sample}_nomt_R2.fq
done
```
Details:
- ref=$MITOREF: path to the reference mitogenome
- in1= path to reads (forward)
- in2= path to reads (reverse)
- basename= path and name of the file with the mitogenome (output name)
- out1= path and name of the nuclear reads (forward)
- out2= path and name of the nuclear reads (reverse)

Move the the splited mitognomes and reference mitogenome to a new folder called "mitogen_03".

# Step 4: Alignment against the reference genome | BWA  (Burrows-Wheeler Aligner)
BWA is a software package for mapping low-divergent sequences against a large reference genome.

### Genome indexing
Before aligning the reads to the reference genome, it is necessary to index the reference genome. 
The reference genome was used: GCF_015237465.2

#### BWA Indexing of Reference Genome

Script:
``` bash
# BWA index the reference genome
bwa index -p $REF/reference_genome.fasta
```
Details:
- p: Output database prefix
- The results should be in the “reference” folder and must be 5 files:
    - *.amb
    - *.ann
    - *.bwt
    - *.pac
    - *.sa

#### Samtools Indexing of Reference Genome
The samtools faidx command is used to create an index for a FASTA-formatted reference genome. 
This index allows for efficient retrieval of sequences from specific genomic regions, providing quick access to the underlying nucleotide information. 
The index file created by samtools faidx has the extension ".fai."

``` bash
# Index the reference genome using Samtools
samtools faidx $REF/reference_genome.fasta
```

#### Create Sequence Dictionary for Reference Genome
The CreateSequenceDictionary step generates a sequence dictionary for the reference genome. This dictionary provides information about the reference sequence, including the name, length, and order of each chromosome or contig.

```bash
# GATK CreateSequenceDictionary command
gatk CreateSequenceDictionary -R $REF/reference_genome.fasta --MAX_RECORDS_IN_RAM 17000
```

### Genome aligment
To align our reads against the reference genome we used bwa mem.
It is the primary and most commonly used algorithm in BWA for aligning high-quality sequencing reads, particularly from next-generation sequencing (NGS) platforms like Illumina. Samtools provides a set of utilities that allow users to manipulate and analyze data in SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) formats.

Script:
```
# BWA read alignment
cat list.txt | while read sample
do
  bwa mem -t 16 $REF/reference_genome \
  $SPLIT/${sample}_nomt_R1.fq $SPLIT/${sample}_nomt_R2.fq \
  | samtools view -bS -@ 16 - > $BAM/${sample}.bam
done
```
Details:
- This script aligns paired-end reads to the reference genome using BWA-MEM, and outputs the results in BAM format using Samtools.
- -b: output is delivered in .sam format
- -S: Automatically detect input

# Step 5: Watch statistics and Sort by genome coordinates
### Samtools Flagstat - Alignment Statistics
The samtools flagstat command is used to generate simple statistics from a BAM file. 
These statistics provide information about the number and types of reads in the alignment file based on their alignment flags. 
The output includes details about how many reads are properly paired, how many are singletons, and various other categories.

Script:
```
# Generate alignment statistics using Samtools flagstat
cat list.txt | while read sample
do
  samtools flagstat $BAM/${sample}.bam > $METALIGN/${sample}_stats.txt
done
```
Details: This is an example of what we can see in the stats.txt file of a sample:
  ```
  25781849 + 0 in total (QC-passed reads + QC-failed reads)
  22603 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  25781849 + 0 mapped (100.00% : N/A)
  25759246 + 0 paired in sequencing
  12883130 + 0 read1
  12876116 + 0 read2
  25759246 + 0 properly paired (100.00% : N/A)
  25759246 + 0 with itself and mate mapped
  0 + 0 singletons (0.00% : N/A)
  0 + 0 with mate mapped to a different chr
  0 + 0 with mate mapped to a different chr (mapQ>=5)
  ```

### Sorting BAM Files
The samtools sort command is used to sort a BAM file by genomic coordinates. 
Sorting is a necessary step before many downstream analyses, as it allows for efficient and ordered access to the aligned reads.

Script:
```
# Samtools sort command
cat list.txt | while read a
do
  samtools sort -o $SORTBAM/${a}_sorted.bam -@ 5 $BAM/${a}.bam
done
```
Details:
- -o: leave the results in the sorted folder with this name
- -@ 5: the number of cpu to use

# Step 6: Cleaning bam files

### Identify and mark duplicate reads in BAM files
Picard MarkDuplicates is a command from the Picard Tools suite, which is a collection of command-line tools for manipulating high-throughput sequencing data. 
The MarkDuplicates tool is specifically designed to identify and mark duplicate reads in BAM files, that may have originated from the same original DNA fragment. 
These duplicates can arise during library preparation or sequencing and may impact downstream analyses, particularly in variant calling or other applications where unique reads are essential.
The metrics file (metrics.txt) provides information about the number and types of duplicates found during the process.

Script:
```
# Picard MarkDuplicates command
cat list.txt | while read a
do
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  -I $SORTBAM/${a}_sorted.bam \
  -O dedup/${a}_sorted.dedup.bam \
  -METRICS_FILE metrics/dedup/${a}_metrics.txt \
  -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true \
  -ASSUME_SORT_ORDER coordinate -MAX_RECORDS_IN_RAM 78000
done
```
Details:
- METRICS_FILE: .txt file with the deduplication metrics
- VALIDATION_STRINGENCY=LENIENT: Picard is less strict when validating data, will be permissive in terms of input validation.
- CREATE_INDEX true: creates an index file for the deduplicated BAM file.
- CREATE_MD5_FILE true: creates an MD5 checksum file for the deduplicated BAM file.
- TAGGING_POLICY All: Specifies the tagging policy for duplicate reads. In this case, it's set to "All," which means all duplicate reads will be marked.
- ASSUME_SORT_ORDER coordinate: Assumes that the input BAM file is sorted in coordinate order. It is important to correctly specify the sorting order to ensure accurate duplicate marking.

### Cleaning the BAM files
The CleanSam tool of Picard, cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads.
```
#------Set relative path--------------------
cat list.txt | while read sample
do
  java -jar $EBROOTPICARD/picard.jar CleanSam \
    -I ${sample}_sorted.dedup.bam \
    -O ${sample}_sorted.dedup_clean.bam
done
```
### Add or Replace Read Groups with Picard
Read groups are essential for tracking sequencing data, and this step ensures each BAM file has the correct metadata associated with it.
```bash
cat list_RG.txt | while read sample
do
  java -jar picard.jar AddOrReplaceReadGroups \
    I=${sample}_sorted.dedup_clean.bam \
    O=${sample}_sorted.dedup_clean_RG.bam \
    RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${sample} \
    CREATE_INDEX=True
done
```

# Step 7: Indel Realignment with GATK and Samtools Indenxing

The Indel Realignment step is important because it helps correct misalignments caused by insertions and deletions (indels) in sequencing data.

### Create Indel Realignment Targets with GATK
This step identifies regions where indels (insertions/deletions) might cause misalignment and creates a list of target regions for realignment.
```bash
cat list_RG.txt | while read sample
do
  gatk3 -T RealignerTargetCreator \
    -R reference_genome.fasta \
    -I ${sample}_sorted.dedup_clean_RG.bam \
    -o ${sample}_realignment_targets.intervals \
    -drf BadMate
done
```
### Perform Indel Realignment with GATK
Using the list of target regions, this step realigns the reads in the BAM file to correct misalignments caused by indels.
```bash
cat list_RG.txt | while read sample
do
  gatk3 -T IndelRealigner \
    -I ${sample}_sorted.dedup_clean_RG.bam \
    -R reference_genome.fasta \
    -targetIntervals ${sample}_realignment_targets.intervals \
    -o ${sample}_realigned.bam \
    -consensusDeterminationModel USE_READS
done
```
### Index the BAM Files with SAMtools
This step creates an index (.bai) for the BAM files, which is necessary for efficient access and downstream processing.
```bash
cat list.txt | while read sample
do
  samtools index ${sample}_realigned.bam
done
```

# Step 8: Calculate coverage with Mosdepth
The calculation of coverage is crucial for the analyses that will be performed later in the downstream pipeline, especially in ANGSD. Coverage depth directly impacts the quality and accuracy of variant calling, population genetics statistics, and other genomic analyses. Ensuring proper coverage helps to minimize biases, avoid missing data, and improve the overall reliability of the results obtained from ANGSD and other downstream tools.

### Calculate Coverage with Mosdepth (Pre-Deduplication)
This step calculates the depth of coverage for each BAM file before removing duplicate reads, which helps assess the overall sequencing quality.
```bash
cat list.txt | while read sample
do
  singularity run mosdepth.sif mosdepth -n --fast-mode -t 2 -Q 20 \
    ${sample}_prededup_coverage \
    ${sample}_sorted.bam
done
```
### Calculate Coverage with Mosdepth (Post-Deduplication)
After indel realignment and other cleaning steps, this command calculates the coverage again to verify sequencing depth in the cleaned files.
```bash
cat list_RG.txt | while read sample
do
  singularity run mosdepth.sif mosdepth -n --fast-mode -t 2 -Q 20 \
    ${sample}_postdedup_coverage \
    ${sample}_realigned.bam
done
```

