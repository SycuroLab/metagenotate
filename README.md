# metagenotate

Snakemake pipeline for metagenomic analysis, annotation and classification. assembly, binning, bin refinement

## Overview

Input: 

* Filtered and cleaned paired-end reads from shotgun metagenome sequencing.

Output: 

Metagenome assembly with metaspades and megahit (metaWRAP). Assembles metagenome using metaspades with filtered and cleaned reads as input. Reads that were not used in the assembly are then assembled using megahit. Both assemblies are then combined into one metagenome assembly file.

Metagenome bins (MAGs)

Metagenome assembly and metagenome bin annotation with Prokka and MetaErg.



## Pipeline summary

### Steps


## Installation

To use this pipeline, navigate to your project directory and clone this repository into that directory using the following command:

```
git clone https://github.com/SycuroLab/metagenotate.git
```

Note: you need to have **conda** and **snakemake** installed in order to run this. To install conda, see the instructions [here](https://github.com/ucvm/synergy/wiki). 

To install snakemake using conda, run the following line:

```
conda install -c bioconda -c conda-forge snakemake
```

See the snakemake installation [webpage](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for further details.

## Config file

All the parameters required to run this pipeline are specified in a config file, written in yaml. See/modify the provided example file with your custom parameters, called `config.yaml`. This is the only file that should be modified before running the pipeline. Make sure to follow the syntax in the example file in terms of when to use quotations around parameters.

## Data and list of files

Specify the full path to the directory that contains your data files in the config file. You also need to have a list of sample names which contains the names of the samples to run the pipeline on, one sample per line. You can run this pipeline on any number or subset of your samples. Sample names should include everything up to the R1/R2 (or 1/2) part of the file names of the raw fastq files. Specify the path and name of your list in the config file.

## Description of parameters
| Parameter | Description | Example |
| -------------- | --------------- | ------------ |
| list_files | Full path and name of your sample list. | `"/home/aschick/project/list_files.txt"` |
| path | Location of input files. | `"/home/aschick/project/data/filtered/"` |
| output_dir | Location of output files. | `"/home/aschick/project/metaphlan"` |
| metaphlan_database | Location of the Metaphlan bowtie2 database. | `"/bulk/IMCshared_bulk/shared/dbs/metaphlan3"` |
| threads | Number of threads for Metaphlan | `8` |
| paired | Are reads paired? Set to `TRUE` if the reads are paired. | `FALSE` |
| for | If paired, suffix of forward reads. | `"_R1_filtered.fastq"` |
| rev | If paired, suffix of reverse reads. | `"_R2_filtered.fastq"` |
| suff | If unpaired, suffix of reads. | `"_merged_trimmed_filtered.fastq"` |

## Running the pipeline on ARC (SLURM cluster)

Test the pipeline by running `snakemake -np`. This command prints out the commands to be run without actually running them. 

To run the pipeline on the ARC compute cluster, enter the following command from the project directory:

```
sbatch < metagenotate_sbatch.sh
```

The above command submits jobs to ARC, one for each sample and step of the metqc pipeline.

Note: the file `cluster.json` contains the parameters for the SLURM job submission system that ARC uses. In most cases, this file should not be modified. Use the `cluster.json` file in the `cluster_files/slurm_files/` folder. 

The ARC Cluster Guide can be found here:
https://rcs.ucalgary.ca/index.php/ARC_Cluster_Guide

The General Guidelines and Policies can be found here:
https://rcs.ucalgary.ca/index.php/General_Cluster_Guidelines_and_Policies



## Results and log files

Snakemake will create a directory for the results of the pipeline as well as a directory for log files. Log files of each step of the pipeline will be written to the `logs` directory.

## Databases and Location

Databases used for metagenotate were downloaded and installed previously for ease of use.

Location: 
`arc.ucalgary.ca`

Directory Path: 
``

If there is a newer version of the database that you want to use for your project you can download the newer version;

```

Place this path in the `` parameter in the `config.yaml` file.


## Known Issues

1) The metaerg_refined_bin rule exits with an error that all the bins failed. Looks like a filehandle error. Solution: Delete all metaerg directories and re-run the pipeline.
