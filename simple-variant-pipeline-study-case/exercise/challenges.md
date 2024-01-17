# Challenges

Once you finished these parts you can continue with two additional challenges. These are optional

**Challenge 1**: Change your script to use the data directly from

```bash
/commons/Themas/Thema11/Dataprocessing/data 
```

without typing the full path every time you need the path it in the script. Upload the script to the repository

**Challenge 2**: Change `bwa mem` to `bowtie2`. Search on the internet for bowtie manuals how to do this. Compare the outcomes. What is there a difference? How about the performance? Measure this by adding a benchmark. With the benchmark directive, Snakemake can be instructed to measure the wall clock time of a job. Activate benchmarking for the rule like below:

```python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```
