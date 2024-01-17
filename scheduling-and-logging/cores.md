# Cores

<pre><code><strong>fennaf@idefix:~$ lscpu | egrep 'Thread|Core|Socket|^CPU\('
</strong>CPU(s):                24
Thread(s) per core:    2
Core(s) per socket:    6
Socket(s):             2
</code></pre>

When a workflow is executed, the number of threads the jobs need is considered by the Snakemake scheduler. In particular, the scheduler ensures that the sum of the threads of all running jobs does not exceed a given number of available CPU cores. This number can be given with the –cores command line argument **(per default, Snakemake uses only 1 CPU core).**

```
snakemake --cores 12
```

Executes the workflow with **12 cores** if 12 cores are available.

\
Since the rule bwa\_map needs 8 threads, only one job of the rule can run at a time, and the Snakemake scheduler will try to saturate the remaining cores with other jobs like, e.g., samtools\_sort. The threads directive in a rule is interpreted as a maximum: **when fewer cores than threads are provided, the number of threads a rule uses will be reduced to the number of given cores**\


```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    threads: 6
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

```
snakemake --cores 12
```

Rule bwa\_map needs 6 threads now. So with 12 cores _two_ jobs of rule bwa\_map can run parallel in case of a single thread per core.

```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    threads: 6
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

```
snakemake --cores 2
```

\
Rule bwa\_map needs 6 threads now. But only 2 cores are provided. It will execute only 1 job with using two threads.
