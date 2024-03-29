# Installing SnakeMake

<pre><code>Snakemake works best when installed through [<a data-footnote-ref href="#user-content-fn-1">Conda</a>](https://docs.conda.io/en/latest/). 
First, install Conda and create a new environment. 
It is best practice to create a new environment for each new project.
</code></pre>

```bash

# Create a new conda environment
conda create -n new_snake_env

# activate this environment
conda activate new_snake_env

# install snakemake from the bioconda channel
conda install -c conda-forge -c bioconda snakemake
```

To deactivate the virtual environment you can type

```bash
conda deactivate
```

To see a list of available environments, type:

```bash
conda env list
```

See the [Conda cheatsheet](https://docs.conda.io/projects/conda/en/stable/\_downloads/843d9e0198f2a193a3484886fa28163c/conda-cheatsheet.pdf) for the most used commands.

[^1]: 
