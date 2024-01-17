# Modularisation

Snakemake is python based and therefor we can import python libraries using `import`. If we want to import other snakefile scripts however we need the `include` directive.&#x20;

```
include:
```

It works as follow:

```python
'''
Main script to be executed.
Includes all necessary Snakefiles needed
to execute the workflow.
'''


rule all:
    input:
        "results/summary.vchk"
        
        
configfile: "config.yaml"
include: "gatk.smk"
include: "fasta.snk"
include: "bam.smk"
include: "bcftools.smk"

```

Best practice is to group rules into logical groups and create a snakefile for each group of rules. More reading see:

{% embed url="https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#" %}
