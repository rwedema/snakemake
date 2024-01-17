# Data paths

So far we used our working directory with the Snakefile for input and output files. We defined input files along with their subdirectory

```python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
        
```

All paths in the snakefile are interpreted relative to the directory snakemake is executed in. This behaviour can be overridden by specifying a workdir in the snakefile:

```python
workdir: "path/to/workdir"
```

So far we only used the path of our working directory. But if we decide to move our data we can change the data paths in the script anywhere. Better practise is to use variables declared in the first lines of the script for the data paths. Since snakemake is python based we can work with string concatenations or joins.

```python
#concatenation example
wdir = "/data/storix2/student/thema11/jdoe/"
FASTQ_DIR = 'data/'

rule quantify_genes:
    input:
        genome = wdir + 'genome.fa'
        r1 = wdir + FASTQ_DIR + '{sample}.R1.fastq.gz',
        r2 = wdir + FASTQ_DIR + '{sample}.R2.fastq.gz'
    output:
        wdir + '{sample}.txt'
    shell:
        'echo {input.genome} {input.r1} {input.r2} > {output}'
```

```python
#join example
from os.path import join

wdir = "/data/storix2/student/thema11/jdoe/"
FASTQ_DIR = "data/"

rule quantify_genes:
    input:
        genome = join(wdir,'genome.fa')
        r1 = join(wdir, FASTQ_DIR, '{sample}.R1.fastq.gz'),
        r2 = join(wdir, FASTQ_DIR, '{sample}.R2.fastq.gz')
    output:
        '{sample}.txt'
    shell:
        'echo {input.genome} {input.r1} {input.r2} > {output}'
```

Best practise of course is not to use a variable wdir for our working directory but the `workdir:` feature

```python
#concatenation example with workdir
workdir: "/data/storix2/student/thema11/jdoe/"
FASTQ_DIR = 'data/'

rule quantify_genes:
    input:
        genome = 'genome.fa'
        r1 = FASTQ_DIR + '{sample}.R1.fastq.gz',
        r2 = FASTQ_DIR + '{sample}.R2.fastq.gz'
    output:
        '{sample}.txt'
    shell:
        'echo {input.genome} {input.r1} {input.r2} > {output}'
```
