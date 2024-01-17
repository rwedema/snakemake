# Step 4: variants

The next step in our workflow will aggregate the mapped reads from all samples and jointly call genomic variants on them (see Background). For the variant calling, we will combine the two utilities `samtools` and `bcftools`.

```bash
samtools mpileup -g -f data/genome.fa sorted_reads/A.bam | \
bcftools call -mv - > calls/all.vcf
```

where the genome and bam file is input and the all.vcf is the output. We can translate this in the following rule:

```python
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```

With multiple input or output files, it is sometimes handy to refer them separately in the shell command. This can be done by specifying names for input or output files (here, e.g., fa=...). The files can then be referred in the shell command via, e.g., `{input.fa}`. For long shell commands like this one, it is advisable to split the string over multiple indented lines. Python will automatically merge it into one. Furthermore the special function `expand()` is used. It takes a string like `"sorted_reads/{sample}.bam"` and expands it into a list like `["sorted_reads/A.bam","sorted_reads/B.bam", "sorted_reads/C.bam"]`. Remember to put `SAMPLES = ["A", "B", "C"]` on top of the Snakefile.&#x20;

**Exercise** Obtain the updated DAG of jobs for the target file `calls/all.vcf`.
