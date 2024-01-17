# Exercise

## Create a fake workspace with FASTQ files

For demonstration purposes, we create some fake test data. With the fake data, we will run the pipeline script. To create the data type in your terminal, do the following:

```bash
# make your workdir
mkdir snakemake-test
mkdir snakemake-test/data

# make fake genome
touch snakemake-test/data/genome.fa

# make fake data
touch snakemake-test/data/Sample1.R1.fastq.gz snakemake-test/data/Sample1.R2.fastq.gz
touch snakemake-test/data/Sample2.R1.fastq.gz snakemake-test/data/Sample2.R2.fastq.gz
```

## Create the snakefile

Open your preferred text editor, type the code below, and save it into your work directory  (NB. In the example above our work directory is snakemake-test and the pipeline script is called Snakefile so we store the script in `snakemake-test/Snakefile`

```python
SAMPLES = ['Sample1', 'Sample2']

rule all:
    input:
        expand('results/{sample}.txt', sample = SAMPLES)

rule quantify_genes:
    input:
        genome = 'genome.fa',
        r1 = 'data/{sample}.R1.fastq.gz',
        r2 = 'data/{sample}.R2.fastq.gz'
    output:
        'results/{sample}.txt'
    shell:
        'echo {input.genome} {input.r1} {input.r2} > {output}'
```

### Running the pipeline

We can run the pipeline by invoking snakemake. It knows to look for a file called Snakefile. Otherwise, you can specify a file to use with the --snakefile option.

```python
snakemake --snakefile Snakefile
```

```bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	2	quantify_genes
	3

rule quantify_genes:
    input: genome.fa, data/Sample1.R2.fastq.gz, data/Sample1.R1.fastq.gz
    output: results/Sample1.txt
    jobid: 1
    wildcards: sample=Sample1

Finished job 1.
1 of 3 steps (33%) done

rule quantify_genes:
    input: genome.fa, data/Sample2.R2.fastq.gz, data/Sample2.R1.fastq.gz
    output: results/Sample2.txt
    jobid: 2
    wildcards: sample=Sample2

Finished job 2.
2 of 3 steps (67%) done

localrule all:
    input: results/Sample1.txt, results/Sample2.txt
    jobid: 0

Finished job 0.
3 of 3 steps (100%) done
```

As you can see first the rule quantify\_genes is executed for Sample1, then quantify\_genes is executed for Sample2 and then rule all is executed

Here are the output files that were created:

```bash
head Sample?.txt
```

```bash
==> Sample1.txt <==
genome.fa data/Sample1.R1.fastq.gz data/Sample1.R2.fastq.gz

==> Sample2.txt <==
genome.fa data/Sample2.R1.fastq.gz data/Sample2.R2.fastq.gz

```

## Run with different options

Run the same snakefile with the following options.&#x20;

```bash
snakemake --snakefile Snakefile
snakemake --cores 4
snakemake -n 
snakemake -p 
snakemake --reason 
snakemake --list 
snakemake --verbose 
snakemake --summary 
```

## Add rule clean

If you try to rerun you will see that `snakemake` does not rerun the jobs. Instead it says:

```
nothing to be done
```

This is because the files that need to be generated already exist. It is one of the features of `snakemake` that it only processes changed or new files. A bit annoying in this development and testing phase. A way to deal with this is adding a rule clean:

```python
rule clean:
    shell:
        'rm *.txt'
```

In the terminal I now can call the rule from my command line

```python
snakemake clean
```

The same way I can test a single rule as well.

## --delete-all-output

If you want to delete all output you can use `snakemake --delete-all-output`

## Add third sample

* Try to extend your sample script with Sample3
* Try to extend your sample script with a final output file test.txt containing the text "Sample1.txt Sample2.txt and Sample3.txt are successfully processed" by the use placeholders {sample}.

## Try different command line options

Read the Run snakemake section. Try the different printing options&#x20;
