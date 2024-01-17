# Step 5: report

Later on in this course we will use R for analysing and reporting the data. For now we will use a built in reporting function called `report`. It generates a html file.

```python
rule report:
    input:
        "calls/all.vcf"
    output:
        "out.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as f:
            n_calls = sum(1 for line in f if not line.startswith("#"))

        report("""
        An example workflow
        ===================================

        Reads were mapped to the Yeas reference genome 
        and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], metadata="Author: Mr Pipeline", T1=input[0])
      
```

In the example above, we notice that this rule does not entail a shell command. Instead, we use the run directive, which is followed by Python code.

First, we import the report function from `snakemake.utils`. Second, we open the VCF file by accessing it via its index in the input files (i.e. input\[0]), and count the number of non-header lines (which is equivalent to the number of variant calls). Third, we create the report using the report function. The function takes a string that contains RestructuredText markup. In addition, we can use the familiar braces notation to access any Python variables (here the samples and `n_calls` variables we have defined before). The second argument of the report function is the path were the report will be stored (the function creates a single HTML file). Then, report expects any number of keyword arguments referring to files that shall be embedded into the report. Importantly, you can refer to the files from within the report via the given keywords followed by an underscore (here T1\_).

## NB

It might be the case that you have to install the some additional packages. You can install these with pip3 install
