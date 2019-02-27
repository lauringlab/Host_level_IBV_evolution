# Within and between-host evolution of influenza B virus in acute infections

This repository holds code for the data analysis and relies heavily on our other repository [variant_pipeline](https://github.com/lauringlab/variant_pipeline).

# Overview
--------

    project
    |- README          # the top level description of content
    |
    |- data            # raw and primary data, are not changed once created
    |  |- metadata/  # data files containing sample metadata from HIVE
    |  |- reference/  # reference fasta files to be used in in sequence alignment
    |  |- raw/         # raw data, will not be altered. 
    |  |- processed/     # cleaned data.
    |- scripts/           # all analysis code
    |  |- primary_analysis/    # The pbs scripts and option files used to process the sequencing data from fastq format to variant calls
    |  |- secondary_analysis/  # R scripts used for variant analysis.
    |  |- figures # R scripts used to make the final figures
    |- results         # all output from workflows and analyses
    |  |- figures/     #  final figures
    |  |- plots/     # intermediate plots that may not be included in final figures
    |  |- *ipynb and *py    # Python notebooks used to run a few steps of the analysis - generating trees, etc.
    +- Makefile        # Executable makefile for this study. Secondary analysis and figure generation.
    
  --------
# Dependencies    

The analysis expects the variant_pipeline repository to be your home directory. See the tutorial in the variant_pipeline repository for setup. 

# Reproducing the analysis

Reproducing the analysis involves 4 steps:

1) Downloading the fastq files from the SRA 
2) Primary analysis - Calling iSNV in each sample using the variant_pipeline repository referenced above
3) Secondary analysis - maniputating the data and running the models 
4) Making the figures

All the intermediate files needed to run the secondary analysis are included in this repository. The Makefile can be used to run the secondary analysis. 

Due to space limitations, the raw fastq files and intermediate bam files are not included here but can be remade using the commands below.

Please note that in many places we refer to the genomic segment "NA" as "NR". This is to avoid complications in R as "NA" is a special term. 
