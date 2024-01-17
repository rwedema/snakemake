# Benchmark

## Benchmarking

With the `benchmark` directive, Snakemake can be instructed to measure the **wall clock time of a job**. In theory our job should run faster if we use multiple threads.  In the script below we use 20 threads.&#x20;

```python
rule bwa_map:
    input:
        join(FDIR, "genome.fa"),
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    threads: 20
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

If we now run snakemake it will generate benchmark files for each sample

```
A.bwa.benchmark.txt
B.bwa.benchmark.txt
C.bwa.benchmark.txt
```

the content of these files looks like this(depending on your machine):

```
s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load
1.0691	0:00:01	28.57	271.92	20.00	21.50	0.00	1.00	0.00
```

## Report benchmark

You can include the benchmark files in your report for more access convenience. Just add a T2 variable.

```python
rule report:
    input:
        T1="calls/all.vcf",
	    T2=expand("benchmarks/{sample}.bwa.benchmark.txt", sample=SAMPLES)
    output:
        "out.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as f:
            n_calls = sum(1 for line in f if not line.startswith("#"))

        report("""
        An example workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        Benchmark results for BWA can be found in the tables T2_.
        """, output[0], **input)
```
