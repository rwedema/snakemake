# Threads

## Threads

Snakemake has a threads directive. Within a rule, you can specify the number of threads to be used for the processing part. The number of threads can be propagated to the shell command with the familiar braces notation (i.e. `{threads}`). If no threads directive is given, a rule is assumed to need **1 thread**.

Some tools can in theory use multiple threads. You can make Snakemake aware of this, such that the information can be used for scheduling. Add a directive `threads:` to the rule and alter the shell command to `{threads}` to use multiple threads.

```
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

Be aware that if we define threads we also need to define cores to allow multithreading. The example below demonstrates this.

{% hint style="danger" %}
Threads needs --cores to schedule for multithreading!
{% endhint %}

If we re-run snakemake it does not seem to increase performance.  Looking in the log file we see the following:

```
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
```

If 1 core (default) is provided the threads will be scaled down. We have to tell snakemake to use more cores!

```
snakemake --cores 8
```

And now it seems to run much faster!&#x20;
