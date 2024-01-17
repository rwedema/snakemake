# Python scripts

A rule can also point to an external script instead of a shell command or inline Python code, e.g

```python
rule NAME:
    input:
        "path/to/inputfile",
        "path/to/other/inputfile"
    output:
        "path/to/outputfile",
        "path/to/another/outputfile"
    script:
        "path/to/script.py"
```

The script path is always relative to the Snakefile (in contrast to the input and output file paths, which are relative to the working directory). Inside the script, you have access to an object snakemake that provides access to the same objects that are available in the run and shell directives (input, output, params, wildcards, log, threads, resources, config), e.g. you can use `snakemake.input[0]` to access the first input file of the above rule. An example is given below:

### Run directive

So far we used little bits and pieces of Python code in our snakemake script. We used for instance the run directive to run a piece of Python code.

```python
rule collate_outputs:
    input:
        expand('{sample}.txt', sample=SAMPLES)
    output:
        'test.txt'
    run:
        with open(output[0], 'w') as out:
            for i in input:
                sample = i.split('.')[0]
                for line in open(i):
                    out.write(sample + ' ' + line)
```

### Script directive

Using the run directive as above is only reasonable for short Python scripts. As soon as your script becomes larger, it is reasonable to separate it from the workflow definition. For this purpose, Snakemake offers the script directive. Using this principle we can change the original run script of `collate_output` in a script command

```python
rule collate_outputs:
    input:
        expand('{sample}.txt', sample=SAMPLES)
    output:
        'test.txt'
    script:
        'scripts/collate_output.py'
```

```python
#!/usr/bin/env python3

"""collate_output.py"""

with open(snakemake.output[0], 'w') as out:
    for i in snakemake.input:
        sample = i.split('.')[0]
        for line in open(i):
            out.write(sample + ':' + line)
```

The Python script itself makes use of the `snakemake.output` variable and the `snakemake.input` variable. Make sure that you make your `.py` script executable.
