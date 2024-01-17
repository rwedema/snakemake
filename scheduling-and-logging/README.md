# Scheduling and logging

\
**Scheduling**

#### CPU's

Before we optimize our code it is important to get an idea about the capability of the computer to execute jobs in parallel. To get an idea of what is available on threads, cores, sockets and CPU's can type:

```bash
lscpu | egrep 'Thread|Core|Socket|^CPU\('
```

Try this for different machines

* assemblix
* helix
* asterix
* mordor
* your bin machine
* another bin machine
* nuc410

What is the difference between bin201 and nuc410?

\


#### Processor versus core id's

To get a more detailed idea about the cores type

```bash
 grep -E 'processor|core id' /proc/cpuinfo | sed 'N;s/\n/; /'
```

Draw a picture of your own bin machines with sockets, processors and cores like this \
How many threads are available on your machine?

\


#### Benchmarking

With the `benchmark` directive, Snakemake can be instructed to measure the **wall clock time of a job**. In theory, our job should run faster if we use multiple threads. Let us try this. Copy the complete script of Tutorial 2 to a new directory. Incorporate the following in rule `bwa_map`:

```python
rule bwa_map:
    input:
        join(FDIR, "genome.fa"),
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

If we now run snakemake it will generate benchmark files for each sample

```bash
A.bwa.benchmark.txt
B.bwa.benchmark.txt
C.bwa.benchmark.txt
```

the content of these files looks like this(depending on your machine):

```bash
s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load
1.0691	0:00:01	28.57	271.92	20.00	21.50	0.00	1.00	0.00
```

\


#### Threads

Now let's us increase the number of threads.

```python
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

If we rerun snakemake it does not seem to increase performance. The benchmark file has similar values.

```bash
s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load
1.0734	0:00:01	28.55	271.92	19.91	21.41	0.00	1.00	0.00
```

What happened? Looking in the log file we see the following:

```
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
```

We have to tell snakemake to use more cores!

```bash
snakemake --cores 4
```

And now it seems to run much faster! Look ath your benchmark files. These confirm this.

```
s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load
0.4817	0:00:00	2.94	10.93	0.27	0.81	0.00	0.00	0.00
```

\


#### Report benchmark

You can include the benchmark files in your report for more access convenience. Just add a T2 variable.

```python
rule report:
    input:
        T1="calls/all.vcf",
	    T2=expand("benchmarks/{sample}.bwa.benchmark.txt", sample=SAMPLES)
    output:
        "out.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as f:
            n_calls = sum(1 for line in f if not line.startswith("#"))

        report("""
        An example workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        Benchmark results for BWA can be found in the tables T2_.
        """, output[0], **input)
```

\


#### Exercise increase performance

Choose your slowest Snakemake script which you made so far. Try to increase the performance with threads and cores. Check your performance with benchmark. First try to increase your performance on your bin machine. What will happen if you run snakemake with more cores than available? What will happen if you run snakemake with more threads than available? Now run the script on idefix. What is the maximum threads and cores you can use? How will it affect other jobs if you ask for maximum CPU sources? How can you avoid overload of memory usage?

#### Running a script on SLURM

You can run a snakemake script on a slurm cluster

```
snakemake --cluster "sbatch -A {cluster.account} -t {cluster.time} -p {cluster.partition} -N {cluster.nodes}" --cluster-config cluster_config.yml —jobs NUM_JOBS_TO_SUBMIT

```

The information within the `{}` are the parameters that snakemake will read from the `cluster_config.yml` file.

```
# cluster_config.yml - cluster configuration
__default__:
    account: ACCOUNT
    partition: PARTITION # for example assemblix
    time: 01:00:00 # time limit for each job
    nodes: 1
    ntasks-per-node: 14 #Request n cores be allocated per node.
    chdir: /directory/to/change/to #working dir
    output: a_name_for_my_job-%j.out
    error: a_name_for_my_job-%j.err
```

{% embed url="https://hackmd.io/@bluegenes/BJPrrj7WB" %}

#### Running a script on Condor (obsolete for BIN grid)

You can run a snakememake script on the condor cluster.&#x20;

```
condor_status -master
```

gives an overview of condor machines. First you need to make a make a bash file in which you execute snakemake. The bash file is a wrapper in which you activate and deactivate the virtual environment.

```
#!/bin/bash
source venv/bin/activate
snakemake --cores 16
deactivate
```

Lastly you configure a `condor.cfg` file

```
####################
# Condor config
####################

Executable     = task.sh
Universe = vanilla
Log     = /homes/fennaf/tmp/condor.log
Error   = /homes/fennaf/tmp/condor.error
Output =  /homes/fennaf/tmp/condor.out
Queue
```

To execute the condor script type in your terminal

```
condor_submit condor.cfg
```

#### Exercise condor

Try to run a snakemake script on Slurm. Commit the scripts and the log and error file on your repository. If you encounter error related to paths and tools try a more simple script.

#### Review exercise

Check the [the following script](https://bitbucket.org/nmkeur/data-processing)

* Why do you think he used different threads for different rules
* What would you recommend the writer to improve his script

\


### Logging

In the script you saw a `log` directive. When executing a large workflow, it is usually desirable to store the output of each job persistently in files instead of just printing it to the terminal. For this purpose, Snakemake allows to specify log files for rules. Log files are defined via the `log` directive and handled similarly to output files, but they are not subject of rule matching and are not cleaned up when a job fails. We modify our rule bwa\_map as follows:

```
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem  -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

The shell command is modified to collect STDERR output of both bwa and samtools and pipe it into the file referred by {log}. Log files must contain exactly the same wildcards as the output files to avoid clashes. It is best practice to store all log files in a subdirectory logs/, prefixed by the rule or tool name.

This logging is especially useful to check afterwards if your script runned as expected. After execution you can check all the updated files output files and logfiles with the command

```
snakemake --summary
```

#### Exercise

Incorporate for every data processing rule an error log with the `log` directive. Run your script. When it is finished run:

```bash
snakemake --summary
```
