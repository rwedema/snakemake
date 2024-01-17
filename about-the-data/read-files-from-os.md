# Read files from os

### Retrieve files function

Using config files is one way to work with flexible input, but a long list of files might be not easy to maintain in a config file. One way to deal with that is to simply retrieve the files from the os. In snakemake we can define functions and use such a function as input:

```python
def myfunc(wildcards):
    return [... a list of input files depending on given wildcards ...]

rule:
    input: myfunc
    output: "someoutput.{somewildcard}.txt"
    shell: "..."
```

### glob.glob

The glob module finds all the pathnames matching a specified pattern according to the rules used by the Unix shell. Using `glob.glob` is useful in the list comprehension to retrieve files.

```python
import glob
import os

SAMPLES = [os.path.basename(x).rstrip(".fastq") for x in glob.glob("data/samples/*")]
print(SAMPLES)
```

###
