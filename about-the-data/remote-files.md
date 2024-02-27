# Remote files

Snakemake can access web resources via a several remote services. The Snakefile supports a wrapper function, remote(), indicating a file is on a remote storage provider. Using the remote() function allows us to download files from remote sites or directly access remote files as input files. The remote() wrapper function can be used for remote webservices access like public like NCBI. Also cloud services like Amazon Simple Storage Service, google and dropbox can be accessed via the remote wrapper. Most of the bioinformatics research us NCBI data. But access to public data through HTTP becomes more relevant as well. Consumers participate in large-scale research projects, such as Lifelines in the northern three provinces. These large-scale datasets are stored in databases and made available to researchers and healthcare professionals. To use the remote wrapper we need to import the snakemake.remote library applicable.rt

### HTTP

Below you find an example in which the WGET is used to download files from to the site https://bioinf.nl/\~ronald/snakemake/WC03/data into a current work directory.

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

The `NCBI.remote("file.fasta", db="source database")` will retrieve the fasta file from the NCBI source database. For example, `NCBI.RemoteProvider().remote("file.fasta", db="source database")` will download a FASTA file while `NCBI.RemoteProvider().remote("file.gb", db="source database")` will download a GenBank-format file. Below you find an example of the NCBI provider

```python
import os
from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
NCBI = NCBIRemoteProvider(email="f.feenstra@pl.hanze.nl") # email required by NCBI

rule all:
    input:
        NCBI.remote("AAH32336.1.fasta", db="protein")
    run:
        outputName = os.path.basename("test.fasta")
        shell("mv {input} {outputName}")
```

### remote access

Downloading files is not always needed for the pipeline processing. One might be simply read the file remotely by using the file as input and store only the outputfile. An example is given below. In the example the size of the file is count by the shell command `wc -c`

```python
from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
NCBI = NCBIRemoteProvider(email="f.feenstra@pl.hanze.nl") # email required by NCBI to prevent abuse

rule all:
    input:
        "size.txt"

rule download_and_count:
    input:
        NCBI.remote("AAH32336.1.fasta", db="protein")
    output:
        "size.txt"
    run:
        shell("wc -c {input} > {output}")
```

One can even query NCBI by

```
NCBI.search()
```

Below you find an example of a query searching for Zika virus genomes. Searched is the organism Zika Virus with a Sequence length \[SLEN] between 9000 and 20000, with a publication date \[PDAT] between 2017/03/20 and 2017/03/24

```python
from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
NCBI = NCBIRemoteProvider(email="f.feenstra@pl.hanze.nl") # email required by NCBI to prevent abuse

# get accessions for the first 3 results in a search for full-length Zika virus genomes
# the query parameter accepts standard GenBank search syntax
query = '"Zika virus"[Organism] AND (("9000"[SLEN] : "20000"[SLEN]) AND ("2017/03/20"[PDAT] : "2017/03/24"[PDAT])) '
accessions = NCBI.search(query, retmax=4)

# give the accessions a file extension to help the RemoteProvider determine the
# proper output type.
input_files = expand("{acc}.fasta", acc=accessions)

rule all:
    input:
        "mf.txt"

rule download_and_count:
    input:
        # Since *.fasta files could come from several different databases, specify the database here.
        # if the input files are ambiguous, the provider will alert the user with possible options
        # standard options like "seq_start" are supported
        NCBI.remote(input_files, db="nuccore", seq_start=5000)

    output:
        "mf.txt"
    run:
        shell("cat {input} >> mf.txt")
```

Mind you that the `input_files` variable is defined on the top of the script. The rule download\_and\_count uses this variable. If you dryrun the script above with -np you will see that the content of the input\_files are:

```
KY785484.1.fasta, KY785481.1.fasta, KY785480.1.fasta
```

This content is filled by the expand function which uses the accessions, which queries the three files (retmax = 3). When you run the snakemake (without the dryrun) you see a message like this

```
rule download_and_count:
    input: KY785484.1.fasta, KY785481.1.fasta, KY785480.1.fasta
    output: mf.txt
    jobid: 1

Removing local output file: KY785484.1.fasta
Removing local output file: KY785481.1.fasta
Removing local output file: KY785480.1.fasta
```

So it is downloaded from the database but removed after usage, since it is not defined as an output file.

### More reading

{% embed url="https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html" %}
