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


The analysis requires the helper repository, variant_pipeline for the primary analysis of variant identification. For the secondary analysis involved in the paper, the R and python scripts are available in this repository. 

# Reproducing the analysis

For Adam as we go through the analysis (will be changed for publication!):
