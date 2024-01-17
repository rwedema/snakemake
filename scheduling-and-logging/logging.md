# Logging

## Log directive

In the script, you saw a `log` directive. When executing a large workflow, it is usually desirable to store the output of each job persistently in files instead of just printing it to the terminal. For this purpose, Snakemake allows to specify log files for rules. Log files are defined via the `log` directive and handled similarly to output files, but they are not subject to rule matching and are not cleaned up when a job fails. We modify our rule `bwa_map` as follows:

```python
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem  -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

The shell command is modified to collect `STDERR` output of both bwa and samtools and pipe it into the file referred by `{log}`. Log files must contain exactly the same wildcards as the output files to avoid clashes. It is best practice to store all log files in a subdirectory `logs/`, prefixed by the rule or tool name.

This logging is especially useful to check afterwards if your script runned as expected. After execution you can check all the updated files output files and logfiles with the command

```bash
snakemake --summary
```
