# Temporary files

Output file marked as `temporary` is deleted after all rules that use it as input are completed. The advantage of temporary files is that no unnecessary storage space is used. The disadvantage is that when the workflow needs to be re-run the execution of the file creation consumes time, since the previous run did not store the file.&#x20;

```python
rule NAME:
    input:
        "path/to/inputfile"
    output:
        temporary("path/to/outputfile")
    shell:
        "somecommand {input} {output}"
```

