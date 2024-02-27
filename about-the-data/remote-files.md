# Remote files

Snakemake can access web resources via several remote services. These are available through [plugins](https://snakemake.github.io/snakemake-plugin-catalog/) that can be installed using `conda install.`

Besides these plugins, there are wrappers for common programs from snakemake rules and workflows. See the [wrapper documentation](https://snakemake-wrappers.readthedocs.io/en/stable/) for an extensive overview of the available wrappers.&#x20;

Using the wrapper functions allows us to download files from remote sites or directly access remote files as input files. Using wrapper functions can be used for remote webservices access like public like NCBI. Also, cloud services like Amazon Simple Storage Service, google, and dropbox can be accessed via plugins or wrappers. Most of the bioinformatics research uses NCBI data. But, access to public data through HTTP becomes more relevant as well. Consumers participate in large-scale research projects, such as Lifelines in the northern three provinces. These large-scale datasets are stored in databases and made available to researchers and healthcare professionals.&#x20;



### HTTP

Below you find an example in which the WGET is used to download files from the site https://bioinf.nl/\~ronald/snakemake/WC03/data into a current work directory.

```python
url = "bioinf.nl/~ronald/snakemake/WC03/data/"
SAMPLES = [x for x in 'ABDEFGHIJ']

rule all:
    input:
       expand("{sample}.bam", sample=SAMPLES)

rule download:
    output:
       "{sample}.bam"
    shell:
        "wget {url}{output}"
```

### NCBI

To retrieve files from NCBI we can make use of the [efetch](https://www.ncbi.nlm.nih.gov/books/NBK179288/) program that is available using the [efetch wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/entrez/efetch.html).

The efetch wrapper will retrieve files from the NCBI source database. For example, `id="AAH32336.1",db="protein", format="fasta")`will download a FASTA file while format=`"gb")` will download a GenBank-format file. Below you find an example of the NCBI provider downloading a fasta file.

To start the workflow we have to specify to use Conda using the `--use-conda` argument. In doing so, snakemake will create a conda environment downloading the efetch tool.&#x20;

```bash
// snakemake -c 1 --use-conda --conda-frontend conda
```

The `--conda-frontend conda` in the previous command specifies to use conda instead of mamba.

```python
rule get_fasta:
    output:
        "test.fasta",
    log:
        "logs/get_fasta.log",
    params:
        id="AAH32336.1",
        db="protein",
        format="fasta",
        # optional mode
        mode=None,
    wrapper:
        "v3.4.0/bio/entrez/efetch"
```

### remote access

Downloading files is not always needed for pipeline processing. One might set the file output using the `temporary` function. This will remove the file once it is no longer needed. An example is given below. In the example, the size of the file is counted by the shell command `wc -c`

```python
rule all:
    input:
        "size.txt"

rule download_and_count:
    output:
        temporary("test.fa"),
    params:
        id="AAH32336.1",
        db="protein",
        format="fasta",
        # optional mode
        mode=None,
    wrapper:
        "v3.4.0/bio/entrez/efetch"

rule size:
    input:
        "test.fa",
    output:
        "size.txt"
    run:
        shell("wc -c {input} > {output}"
```
