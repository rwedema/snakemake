# Step 1: mapping

Our first Snakemake rule maps reads of a given sample to a given reference genome. For this, we will use the tool `bwa`, specifically the subcommand `bwa mem`:

```
bwa mem data/genome.fa data/samples/A.fastq | samtools view -Sb - > mapped_reads/A.bam
```

We will not do this in the terminal but writing a Snakefile script. In the working directory, create a new file called Snakefile

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

Snakemake executes a rule if the output file needs to be created. Therefor we can test the script by calling the outputfile

```python
snakemake mapped_reads/A.bam
```

### Generalising

Obviously, the rule will only work for a single sample with reads in the file `data/samples/A.fastq`. However, Snakemake allows to generalise rules by using named wildcards. Simply replace the A in the second input file and in the output file with the wildcard `{sample}`, leading to

```python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    message: "executing bwa mem on the following {input} to generate the following {output}"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

When Snakemake determines that this rule can be applied to generate a target file by replacing the wildcard `{sample}` in the output file with an appropriate value, it will propagate that value to all occurrences of `{sample}` in the input files and thereby determine the necessary input for the resulting job. Note that you can have multiple wildcards in your file paths, however, to avoid conflicts with other jobs of the same rule, all output files of a rule have to contain exactly the same wildcards. When executing

```python
snakemake mapped_reads/{A,B}.bam
```

this will execute the flow for Sample A and B.
