# Editor tips

### Editor tips

Both Pycharm and Visual Code can be used for convenient syntax highlighting.  Search for `snakemake` in the plugin section of the IDE.

In snakemake you can simply use the

```
"""docstrings""" and #comments 
```

to make your Snakefile more readable. Furthermore, you can use messages to inform you a bit more about the process.

```
message: "Executing Samtools on the following files {input}"
```

Later on we will use the message attribute.
