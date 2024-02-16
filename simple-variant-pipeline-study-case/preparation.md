# Preparation

### Environment Preparation

First, you need to activate your work environment again.

```bash
conda activate <name_of_created_conda_env>
```

### Data Preparation

Next, we need the data. It is important to organize our data in subdirectories to avoid overriding original data with processed data. Fetch the data from

```bash
/commons/Themas/Thema11/Dataprocessing/WC02/data
/commons/Themas/Thema11/Dataprocessing/WC02/data/samples
```

and organise as follows in your work directory

```python
data/:
genome.fa
genome.fa.amb
genome.fa.ann
genome.fa.bwt
genome.fa.fai
genome.fa.pac
genome.fa.sa

data/samples:
A.fastq
B.fastq
C.fastq
```
