# Configuration files

Configuration files specify the parameters of an application, such as data paths, thresholds, or filter values. These parameters typically are written in a key-value pair format. Placing the parameter values in a configuration file gives a central and well-organized location to specify what arguments ought to be parsed to a program without hard coding such as in a script. This makes code more user-friendly and flexible.&#x20;

Snakemake provides a configuration file mechanism. Configuration files can be written in JSON or YAML, and loaded with the configuration file directive.&#x20;

```python
configfile: "config.yaml"
```

This mechanism is more flexible since you can use data from several sources.

## YAML

“YAML Ain’t Markup Language” (abbreviated YAML) is a human-readable data serialisation standard that can be used in conjunction with all programming languages and is often used to write configuration files. YAML is designed to be useful and friendly to people working with data. Its features are derived from Perl, C, HTML, and other languages. YAML is a superset of JSON that comes with multiple built-in advantages such as including comments, self-referencing, and support for complex data types. Full documentation for YAML can be found on its [official site](https://yaml.org/spec/1.2/spec.html).

## PyYAML <a href="#pyyaml" id="pyyaml"></a>

PyYAML is a YAML parser and emitter for Python. It can be installed with pip

```
pip install pyyaml
```

## YAML file <a href="#yaml-file" id="yaml-file"></a>

A YAML file is just a flat file value with key : value pairs. YAML supports strings, integers, floats, lists, and associative&#x20;

For our Snakefile we need a separate configuration file with some genome and sample data.

First, we need to create a `config.yaml` file in the working directory listing the files. Instead of expanding a list with the sample names we now expand the samples from the configuration file

```python
vcf=expand("calls/{sample}.g.vcf", sample=config["samples"])
```

Data (A.bam .. J.bam) for this example below is to be found at bioinf.nl/\~ronald/snakemake/WC03/data

```yaml
#configfile config.yaml
genome: data/A
ext : .fa
samples:
    A : /data/A.bam
    B : /data/B.bam
    C : /data/C.bam
    D : /data/D.bam
    E : /data/E.bam
    F : /data/F.bam
    G : /data/G.bam
    H : /data/H.bam
    I : /data/I.bam
```

```python
configfile: "config.yaml"
    
rule merge_variants:
    input:
         fa=config["genome"] + config["ext"],
         fai=config["genome"] + config["ext"] +  ".fai",
         dict=config["genome"] + ".dict",
         vcf=expand("calls/{sample}.g.vcf", sample=config["samples"]),
    output:
        temp("calls/merged_results.vcf")
    message: "Executing GATK CombineGVCFs with {threads} threads on the following files {input}."
    shell:
        "java -jar ./GATK/GenomeAnalysisTK.jar -T CombineGVCFs -R {input.fa} {vcf2} -o {output} "
    
```

The config principle can be used for other wildcard variables as well.

```python
config["genome"]
config["ext"]
```

Using the config file allows you to develop a pipeline without hardcoded data paths.&#x20;

Instead of using the system avail
