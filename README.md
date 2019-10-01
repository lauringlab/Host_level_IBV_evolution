# Within and between-host diversity of influenza B virus in natural infections

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
    |  |- processed/     # cleaned and intermediate data.
    |- scripts/           # all analysis code
    |  |- primary_analysis/    # The pbs scripts and option files used to process the sequencing data from fastq format to variant calls
    |  |- secondary_analysis/  # R scripts used for variant analysis.
    |  |- figures # R scripts used to make the final figures
    |- results         # output from workflows and analyses
    |  |- figures/     #  raw panels for figures
    
  --------

The analysis requires the helper repository, variant_pipeline for the primary analysis of variant identification. For the secondary analysis involved in the paper, the R and python scripts are available in this repository. 

Raw sequencing data is available through the NCBI SRA with BioProject PRJNA561158.

After downloading data into an appropriate directory, the primary analysis can be replicated with the commands in IBV_setup.sh.

After downloading the variants (all.variants.csv), coverage data (all.coverage.csv), and consensus sequences (parsed_fa), the secondary analysis can be replicated with the scripts listed in analysis_script_order.txt. Due to file size limitations, the raw all.variants.csv file is not included in this repository.
