# R scripts

Snakemake utils has a module R which you can use for small R scripts. We have to use `run:`

```python
 run:
    from snakemake.utils import R
    R("""
    some R script
    """)
```

The R script

```r
    data <-read.scv(file = "out.csv", header=FALSE, sep=";")
    jpeg("result/histogram.jpg")
    hist(x=data$V1, 
         main="expression values of gene CCND3 Cyclin D3",
         ylab="count",
         xlab="gene CCND3 Cyclin D3 expression")
    dev.off()
```

Can be implemented as follows:

```python
datadir= '/commons/Themas/Thema11/Dataprocessing/WC04/data/'

rule all:
    """ final rule """
    input: 'result/histogram.jpg'


rule make_histogram:
    """ rule that creates histogram from gene expression counts"""
    input:
        datadir + 'out.csv'
    output:
         'result/histogram.jpg'
    run:
        from snakemake.utils import R
        R("""
        data <-read.csv(file = "{input}", header=FALSE, sep=";")
        jpeg("{output}")
        hist(x=data$V1,
             main="expression values of gene CCND3 Cyclin D3",
             ylab="count",
             xlab="gene CCND3 Cyclin D3 expression")
        dev.off()
        """)
```

The R code is pasted between the

```python
       R("""
       """)
```

Furthermore the `file = "out.csv"` is replaced by `file = "{input}"` and the line `jpeg("result/histogram.jpg")` is replaced by `jpeg("{output}")`

Maybe even more flexible is to make a separate R script:

```r
args = commandArgs(trailingOnly=TRUE)
data <-read.csv(file = args[1], header=FALSE, sep=";")
jpeg(args[2])
hist(x=data$V1,
     main="expression values of gene CCND3 Cyclin D3",
     ylab="count",
     xlab="gene CCND3 Cyclin D3 expression")
dev.off()
```

and calling the script from the Snakefile shell:

```python
datadir = '/commons/Themas/Thema11/Dataprocessing/WC04/data/'

rule all:
    """ final rule """
    input: 'result/histogram.jpg'


rule make_histogram:
    """ rule that creates histogram from gene expression counts"""
    input:
        datadir + 'out.csv'
    output:
         'result/histogram.jpg'
    shell:
        "Rscript scripthist.R {input} {output}"
```

Mind you that the shell command parses `args[1]` and `args[2]` which needs to be fetched from the command line by the R script

Want to test yourself? Use the data from `/commons/Themas/Thema11/Dataprocessing/WC04/data`
