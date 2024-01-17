# Preparation

### Environment Preparation

First, you need to activate your work environment again.

```python
source {name}/bin/activate
```

### Data Preparation

Next, we need the data. It is important to organize our data in subdirectories so that we do not override original data with processed data. Fetch the data from

```bash
/commons/Themas/Thema11/Dataprocessing/data
/commons/Themas/Thema11/Dataprocessing/data/samples
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
