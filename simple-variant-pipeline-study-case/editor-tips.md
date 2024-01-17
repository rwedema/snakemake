# Editor tips

### Editor tips

The Emacs editor can be used for convenient syntax highlighting. Start the editor, and enter on top of the Snakefile:

```python
# -*- python -*-
```

Now save the file, and restart the editor. You will see how the syntax highlights. If you do not work remote you can install editors like atom, visual studio code or pycharm on your system to edit the scripts.&#x20;

In snakemake you can simply use the

```
"""docstrings""" and #comments 
```

to make your Snakefile more readable. Furthermore, you can use messages to inform you a bit more about the process.

```
message: "Executing Samtools on the following files {input}"
```

Later on we will use the message attribute.
