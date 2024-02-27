# Exercises

Before you start: Make sure that you activate your conda environment.&#x20;

### HTTP

Download the https://bioinf.nl/\~ronald/snakemake/test.txt file via a Snakefile.&#x20;

### NCBI

Download the **KY785484.1.fasta** file from the NCBI **nuccore** database via a Snakefile.

### Work directory

Move your data from studycase 02 to a data directory on your commons directory. It is not needed to do this using a snakemake script, you can copy the data. Now rewrite the script from tutorial 02 using the workdir directive.&#x20;

### Configuration file

Write a snakemake script that runs a merge variants workflow merging the data from bioinf.nl/\~ronald/snakemake/WC03/data. Use a configuration file.&#x20;

### Include modules

Split the snakefile into a main snakefile and another snakefile. Import the other rules using includes.

### Databases

Review [https://github.com/veerdonk/data\_processing/blob/master/Snakefile](https://github.com/veerdonk/data\_processing/blob/master/Snakefile)

Answer the following questions:

1. How is the input output construction handled in case of a database?
2. How is de Snakefile documented?
3. What would you improve in this script?
