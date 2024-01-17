# Exercises

### Threads, cores, sockets, and CPU's

Get an idea of what is available on threads, cores, sockets, and CPUs. Try this for different machines

* assemblix
* helix
* asterix
* mordor
* your bin machine
* another bin machine
* nuc410

What is the difference between bin201 and nuc410?

### Increase performance

Choose the slowest Snakemake script that you made so far. Try to increase the performance with threads and cores. Check your performance with a benchmark. First, try to increase your performance on your bin machine. What will happen if you run snakemake with more cores than available? What will happen if you run snakemake with more threads than available? Now run the script on assemblix. What are the maximum threads and cores you can use? How will it affect other jobs if you ask for maximum CPU sources? How can you avoid overload of memory usage?&#x20;

### Log&#x20;

Incorporate for every data processing rule an error log with the `log` directive. Run your script. Inspect the log files.

### Slurm

Try to run a snakemake script on Slurm. Commit the scripts and log files on your repository. If you encounter errors related to paths and tools try a more simple script.

### Code review

Check [the following script](https://bitbucket.org/nmkeur/data-processing)

* Why do you think he used different threads for different rules?
* What would you recommend the writer to improve his script?
