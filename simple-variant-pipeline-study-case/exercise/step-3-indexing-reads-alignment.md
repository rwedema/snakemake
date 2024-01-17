# Step 3: indexing reads alignment

Next, we need to use `samtools` again to index the sorted read alignments for random access. The samtools command looks like this:

```bash
samtools index sorted_reads/A.bam
```

This can be done with the following rule:

```python
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
```

Now test the rule by the snakemake `sorted_reads/{sample}.bam.bai` command like in the previous example and make a new visualisation.
