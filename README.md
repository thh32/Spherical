# Spherical

Spherical is an iterative approach to assembling metagenomic datasets. Spherical has two uses; firstly the reducion of RAM usage as to allow metagenomic assemblies to be produced on computational infrastructure that would otherwise be insufficient, secondly Spherical allows for the assembly and study of low abundance varients otherwise hidden by the more common species.



##Requirements

- Python 2.7
- Velvet
- Bowtie2

Python modules;
- Numpy
- HTSeq



##Usage

Once all dependancies are install and Spherical has been downloaded you can begin using Spherical from command line.

The most basic use of Spherical is;
```
python Spherical.py -fasta -i $INPUT -velvet -bowtie2 -o $OUTPUT

```

###Options

Full customisation of Spherical is possible using the commands below.

| Option command| Description                                                                                              | Default |
| ------------- | -------------------------------------------------------------------------------------------------------- | ------- |
|-align $INT    | Identifies an alignment rate that must be obtained before Spherical stops producing iterations           | 70      |
| -iter $INT    | States the number of iterations that must be produced                                                    | 5       |
| -k $INT       | Sets the Kmer size which should be used for each assembly                                                | 31      |
| -R $INT       | States the fraction of the input which should be used in each iteration of assembly                      | NA      |
| -limit        | A switch which prevents contigs smaller from 300bp from being produced                                   | NA      |
| -m            | A switch which activates the combination of all iterations assemblies at the end to produce a final file | NA      |
| -fasta        | States that the input file is in FASTA input                                                             | NA      |
| -fastq        | States that the input file is in FASTQ input                                                             | NA      |
| -i  $FILE     | Identifies the input file                                                                                | NA      |
| -o $FILE      | Identifies string to be used for output file                                                             | NA      |





