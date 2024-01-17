# Step 2: sorting reads

For later steps, we need the mapped read alignments in the BAM files to be sorted. So the bam file in the mapped\_reads need to be sorted and stored somewhere else. This can be achieved with the `samtools` command, for instance:

```bash
samtools sort -f mapped_reads/A.bam sorted_reads/A.bam 
```

To generalise this we can add a rule to our script.

```python
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```

If we now want to run the rule `samtools_sort` we can call

```python
snakemake sorted_reads/A.bam
```

Because the rule samtools\_sort need as an input mapped\_reads/A.bam it will first execute the rule bwa\_map if no bam file is available to generate the mapped\_reads/A.bam file. To test the workflow we can delete the bam files in the mapped\_reads directory and the bam files in the sorted\_reads directory and run

```python
snakemake sorted_reads/{A,B,C}.bam
```

The output should be like this (if we removed all the files from mapped\_reads and sorted\_reads first):

```bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	3	bwa_map
	3	samtools_sort
	6

rule bwa_map:
    input: data/genome.fa, data/samples/C.fastq
    output: mapped_reads/C.bam
    jobid: 4
    wildcards: sample=C

[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 25000 sequences (2525000 bp)...
[M::mem_process_seqs] Processed 25000 reads in 0.978 CPU sec, 0.980 real sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa mem data/genome.fa data/samples/C.fastq
[main] Real time: 1.281 sec; CPU: 1.019 sec
Finished job 4.
1 of 6 steps (17%) done

rule bwa_map:
    input: data/genome.fa, data/samples/B.fastq
    output: mapped_reads/B.bam
    jobid: 5
    wildcards: sample=B

[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 25000 sequences (2525000 bp)...
[M::mem_process_seqs] Processed 25000 reads in 0.954 CPU sec, 0.955 real sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa mem data/genome.fa data/samples/B.fastq
[main] Real time: 1.252 sec; CPU: 0.992 sec
Finished job 5.
2 of 6 steps (33%) done

rule bwa_map:
    input: data/genome.fa, data/samples/A.fastq
    output: mapped_reads/A.bam
    jobid: 3
    wildcards: sample=A

[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 25000 sequences (2525000 bp)...
[M::mem_process_seqs] Processed 25000 reads in 0.947 CPU sec, 0.948 real sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa mem data/genome.fa data/samples/A.fastq
[main] Real time: 1.257 sec; CPU: 0.991 sec
Finished job 3.
3 of 6 steps (50%) done

rule samtools_sort:
    input: mapped_reads/C.bam
    output: sorted_reads/C.bam
    jobid: 1
    wildcards: sample=C

Finished job 1.
4 of 6 steps (67%) done

rule samtools_sort:
    input: mapped_reads/B.bam
    output: sorted_reads/B.bam
    jobid: 2
    wildcards: sample=B

Finished job 2.
5 of 6 steps (83%) done

rule samtools_sort:
    input: mapped_reads/A.bam
    output: sorted_reads/A.bam
    jobid: 0
    wildcards: sample=A

Finished job 0.
6 of 6 steps (100%) done
```

Let us take a closer look at the resulting DAG of jobs. By executing

```bash
snakemake --dag sorted_reads/{A,B,C}.bam | dot -Tsvg > dag.svg
```

we create a visualisation of the DAG.
