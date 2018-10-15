# Spherical

Spherical is an iterative approach to assembling metagenomic datasets. Spherical has been designed to produce a more complete assembly from deep sequenced metagenomic data. Utilization of multiple iterations of assembly allows for regions which otherwise would be missed to be assembled without a reduction in contig accuracy. Another use for Spherical is its ability to produce metagenomic assemblies using a subset of the initial input file, allowing for assembly of a metagenome whilst using a fraction of the RAM that would otherwise be required. 

## Requirements

- Python 2.7
- Velvet (tested using version 1.2.10)
- Bowtie2 (tested using version 2.2.3)

Python modules;
- Numpy
- HTSeq



## Usage

Once all dependancies are install and Spherical has been downloaded you can begin using Spherical from command line.

The most basic use of Spherical is;
```
python Spherical.py -fasta -i $INPUT -velvet -bowtie2 -o $OUTPUT

```

## Tutorial

To ensure all the dependancies for Spherical are correctly installed and get you used to the various options of Spherical there is a test dataset.

The test dataset consists of a Chicken cecum microbiome dataset obtained from MG-RAST (entry;101).

```
python Spherical.py -fasta -m -k 21 -R 1 -align 99 -iter 5 -i test_data.fa -velvet -bowtie2 -o test_assembly

```
In this command we tell Spherical to produce an assembly using Velvet as the assembler and Bowtie2 as the aligner. The subsample size given of the  input data to Velvet is 1 which means the entire file is to be used. Velvet will use a kmer size of 21 and only finish once either the alignment rate reaches 99% or Spherical has completed 5 iterations of assembly. Finally we have also used the `-m` command to tell Spherical to combine the finished assemblies into a single file.

This will produce an output assembly titled `test_assembly.combined.fa` which when compared using the code below the value `0` should be returned.
```
grep -F -x -v -f test_assembly.combined.fa provided_assembly.fa |grep '>'| wc -l
```

A value other than `0` inicates that an error has occured during the assembly process leading to a difference between the assemblies. This is most likely due to the user using a different version of Velvet or Bowtie2.


### Options


Full customisation of Spherical is possible using the commands below.

| Option command| Description                                                                                              | Default |
| ------------- | -------------------------------------------------------------------------------------------------------- | ------- |
|-align $INT    | Identifies an alignment rate that must be obtained before Spherical stops producing iterations           | 70      |
| -iter $INT    | States the number of iterations that must be produced                                                    | 5       |
| -k $INT       | Sets the Kmer size which should be used for each assembly                                                | 31      |
| -R $FLOAT     | States the input size which should be used in each iteration of assembly e.g. -R 0.25 = 25% of reads     | NA      |
| -limit        | A switch which prevents contigs smaller from 300bp from being produced                                   | NA      |
| -m            | A switch which activates the combination of all iterations assemblies at the end to produce a final file | NA      |
| -fasta        | States that the input file is in FASTA input                                                             | NA      |
| -fastq        | States that the input file is in FASTQ input                                                             | NA      |
| -i  $FILE     | Identifies the input file                                                                                | NA      |
| -o $FILE      | Identifies string to be used for output file                                                             | NA      |
| -bowtie2      | Identifies Bowtie2 as the aligner of choice                                                              | NA      |
| -velvet       | Identifies Velvet as the assembler of choice                                                             | NA      | 
| -x1  $STRING     | Allows additional options for assembly to be used in Velveth or ABYSS steps, used for options starting with '-'                                                                                | NA      |
| -x2  $STRING     | Allows additional options for assembly to be used in Velveth or ABYSS steps, used for options starting with '--'                                                                                | NA      |





