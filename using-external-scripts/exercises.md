# Exercises

## R script or Python script

The file `gene_ex.csv` contains gene expression data of yeast. Your job is to come up with a gene expression heatmap visualization of this file by running a snakemake script using an external R or external Python script.

Hint: the heatmap needs a matrix. In R you can use something like this:

```r
d <- as.matrix(read.csv('gene_ex.csv', header=FALSE, sep=",")[-1,-1])
rownames(d) <- read.csv('gene_ex.csv', header=FALSE, sep=",")[-1,1]
colnames(d) <- read.csv('gene_ex.csv', header=FALSE, sep=",")[1,-1]
```

Hint: first try to create the heatmap with a separate script in R or Python. Then try to incorporate this in the snakemake script. Data for this exercise is to be found at&#x20;

```bash
/commons/Themas/Thema11/Dataprocessing/WC04/data
```

## End-to-end pipeline

The scripts’ job is to take a genomic reference (`reference.fa`) and an input file (`reads.txt`) and output a pileup (`out.vcf`). The pipeline goes through the following stages:

1. Create a BWA index in the genomic reference
2. Align the reads in the input file against the genomic reference
3. Convert the alignment into a `.sam` file
4. Convert the .sam file into a `.bam` file and sort it
5. Detect and remove duplicates
6. Index the results
7. Create the pileup and convert it into a `.vcf` file

The Picard-tools needs to be installed.

* First, clone the repo:

```bash
    git clone https://github.com/broadinstitute/picard.git
    cd picard/
```

* Picard is now built using [gradle](http://gradle.org/). A wrapper script (`gradlew`) is included which will download the appropriate version of gradle on the first invocation.
* To build a fully-packaged, runnable Picard jar with all dependencies included, run:

```bash
    ./gradlew shadowJar
```

* The resulting jar will be in `build/libs`.

The reference sequence and the reads to align can be found in

```bash
    /commons/Themas/Thema11/Dataprocessing/WC04/data
```

named `reference.fa`. and `reads.txt`. The bash commands are listed below.

```bash
#!/bin/bash


#Create a BWA index in the genomic reference
bwa index reference.fa
#Align the reads in the input file against the genomic reference
bwa aln -I -t 8 reference.fa reads.txt > temp/out.sai
bwa samse reference.fa temp/out.sai reads.txt > alligned/out.sam
#Convert the .sam file into a .bam file and sort it
samtools view -S -b alligned/out.sam > temp/out.bam
samtools sort temp/out.bam -o sorted/out.sorted.bam
#Detect and remove duplicate
java -jar picard/build/libs/picard.jar MarkDuplicates \
                            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
                            METRICS_FILE=out.metrics \
                            REMOVE_DUPLICATES=true \
                            ASSUME_SORTED=true  \
                            VALIDATION_STRINGENCY=LENIENT \
                            INPUT=sorted/out.sorted.bam \
                            OUTPUT=filtered/out.dedupe.bam
#Index the results
samtools index filtered/out.dedupe.bam
#Create the pileup and convert it into a .bcf file
samtools mpileup -uf reference.fa filtered/out.dedupe.bam | bcftools view -> result/out.vcf
```

Make a snakemake pipeline that generates the steps from the bash script to be used for several examples. Test first the script with the `reference.fa` and `reads.txt`. If the script does not produce any error you can use any other reference genome and reads to test the multiple samples.

Hints:

* run the individual commands in the command line to see what for each step goes in and comes out
* create the snakemake script step by step. For example first output the index files

```python
rule all:
    input:
        'bwa_index.done'



rule bwa_index:
    input:
        'reference.fa'
    output:
        touch('bwa_index.done')
    shell:
        "bwa index {input}"

    
```

Then add the next step

```python
rule all:
    input:
        'temp/out.sai'



rule bwa_index:
    input:
        'reference.fa'
    output:
        touch('bwa_index.done')
    shell:
        "bwa index {input}"


rule bwa_allign1:
    input:
        check = "bwa_index.done",
        gen = 'reference.fa',
        reads = "reads.txt"
    output:
        "temp/out.sai"
    shell:
        "bwa aln -I -t 8 {input.gen} {input.reads} > {output}"
```

And so on. Then try to figure out how to make the filenames flexible and use multiple reads
