# FAQ

### General tip

* give each rule a different output directory so snakemake cannot become confused

### Tips for testing

* use snakemake dry runs

```bash
snakemake -np
```

* Test your script with small data. For instance the yeast genome. Only **after** the optimization of your script (use of threads and cores and fast machines like assemblix) you are ready to test with the actual data.

### The system has no data storage available, what should I do?

First test with only a small amount of samples, for example, 2 samples. Try to use a small test genome. If your pipeline works well, create files that are 'intermediate' (for example a bam file) temporarily with `temp ()`. Only **then** perform the pipeline with all samples. If you have finished your assignment and you still have many 'intermediate files' on the system, please **clear** them up. **All that is needed for the assessment is the final output (often only one count file) and the log files**.

In case of emergency you can go to `assemblix2019:/local-fs/bachelor-students/2019-2020`. But this is not backed up. If you do please notify Marcel.

### Where should I store my data?

* Please read @assemblix2019:/students LEESMIJ file
* Please read @assemblix2019:/data LEESMIJ-local-fs
* other possibility is @assemblix2019:/data/storix2/students/2019-2020/Thema11
* storage in commons is not allowed!

### Help I got an error

* find out if it is a rule error or a shell command error
* Need more help? [Click here](http://bfy.tw/HU3K)

### Rule error or shell command error?

If your rule generates an error it is always best practice to find out the source of the error. Is it your rule or is it the shell command itself? To find out try to run the shell command directly in your terminal. For instance, if your script is like this:

```python
rule bowtie2_index_bacteria:
    input:
        wdir + "Fastafiles/all_bacteria.fna" 
    params:
        bacteria = wdir + "Indices/bacteria/all_bacteria"
    output:
        expand(wdir + "Indices/bacteria/all_bacteria.{index}.bt2", index=range(1,5)),
        expand(wdir + "Indices/bacteria/all_bacteria.rev.{index}.bt2", index=range(1,3))
    log:
        wdir + "logs/index_bacteria.log"
    shell:
        "(bowtie2-build {input} {params.bacteria}) 2> {log}"
```

Then run in your terminal the command with the actual files you want to use, for example

```
bowtie2-build /data/storage/dataprocessing/Fastafiles/all_bacteria.fna
/data/storage/dataprocessing/Indices/bacteria/all_bacteria  2> logs/index_bacteria.log
```

If the terminal command line code does not return an error then you need to dig into the snakemake script. Otherwise, you need to fix the shell command.

### I fail to fix the shell command. Can I use a similar tool?

First, try to read the manual of the tool carefully. Did you use the tool properly? If you cannot find the error you are allowed to use another similar tool

### Workaround for dependable input-output data

In the case of dependable input-output like this:

```bash
awk '$4>=24' 1.txt > a.txt
awk '$4>=24' 2.txt > b.txt
```

We need a smart solution to only retrieve 1.txt from the input list in the case of a.txt as output. Best practice is to use `config files` and `wildcards` and a `lambda` function.

```python
lambda wildcards: config["sample"][wildcards.sample]
```

The `wildcards.sample` is either a or b and depending on the wildcard it will retrieve the value from the config file.

My config file is as follows:

```
sample:
  a: data/1.txt
  b: data/2.txt
```

My snakefile script is as follow:

```
configfile: "config.yaml"
SAMPLES = ['a', 'b']

rule all:
    input:
        expand('processed/{sample}.txt', sample = SAMPLES)

rule quantify_genes:
    input:
        lambda wildcards: config["sample"][wildcards.sample]
    output:
        "processed/{sample}.txt"
    shell:
        "awk '$4>=24' {input} > {output}"
```

### bowtie index build

Building an index with bowtie needs a workaround because bowtie generates several output files. It needs information to build the indexes. Here is an overview of solutions:

```python
rule create_index:
    """creates a bowtie2 index of the required fasta from config.yaml.
       the wildcard reference_name refers to the key under which the file path is stored in config.yaml
       given fasta(in config.yaml) shoulf be .fa or .fna
    """
    input:
        #wildcards are not accessible like config[{wildcard}]. this lambda is a workaround
        lambda wildcards: config["refgenomes"][wildcards.reference_name]
    output:
        "indices/{reference_name}"
    shell:
        "bowtie2-build {input} {output}/{wildcards.reference_name}  --large-index"
```

```
rule build_indexes:
    input:
        lambda wildcards: config["refgenomes"][wildcards.indexName]
    output:
        temp(expand("indexes/{{indexName}}.{ind}.bt2", ind=range(1,5))),
        temp(expand("indexes/{{indexName}}.rev.{ind}.bt2", ind=range(1,3)))
    shell:
        "bowtie2-build {input} indexes/{wildcards.indexName}"
```

### What does wildcard mean?

The `[wildcards.reference_name]` is basically the same as `[{reference_name}]`

### I need a list of names in my input

You can do this is a function that reads all the files in a certain directory, or you could use a parameter

```python
params:
		filelist = ",".join([file for file in config["files"]])
	shell:
		"({config[tool]}"
		"-thread {threads} "
		"-f {params.filelist} "
```
