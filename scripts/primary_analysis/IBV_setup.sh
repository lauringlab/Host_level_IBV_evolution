
### Author: ALV 6/4/2019
### Purpose: A log of all the commands to go from data on /nfs/turbo to running the pipeline for IBV.
### Working directory: IBV_HIVE

# ========== Set up directory and copy files from permanent storage ===============
cd /scratch/alauring_fluxm/avalesan
mkdir IBV_HIVE
cd IBV_HIVE
mkdir Run_2445
mkdir Run_2446
mkdir Run_2445/data
mkdir Run_2446/data
mkdir Run_2445/data/fastq
mkdir Run_2446/data/fastq

qsub copy_files.pbs

# ============= Merge, rename, and unzip fastq files =============
cd Run_2445/data/fastq
qsub merge_lanes_2445.pbs # need merge_lanes_2445.sh
rm *-*
python ~/variant_pipeline/scripts/change_names_miseq.py -s . -f ../fastq_renamed -run
cd ../fastq_renamed
gunzip -v *gz

cd ../../../Run_2446/data/fastq
qsub merge_lanes_2446.pbs # need merge_lanes_2446.sh
rm *-*
python ~/variant_pipeline/scripts/change_names_miseq.py -s . -f ../fastq_renamed -run
cd ../fastq_renamed
gunzip -v *gz

# =========== Get fastq files into format specified by options files =============
cd ../../../
qsub pipeline_setup.pbs

# =========== Prepare modules for pipeline ===============
cd ~/variant_pipeline
module load muscle bowtie2 fastqc R/3.5.0 python-anaconda2/201704
R
packrat::restore()
quit() # don't need to save workspace
module load R/3.5.0

# ================== Run the pipeline ==============
cd /scratch/alauring_fluxm/avalesan/IBV_HIVE/Run_2445
qsub IBV_2445_VIC.pbs # Will need to change path names in these job files
qsub IBV_2445_YAM.pbs

cd ../Run_2446
qsub IBV_2446_VIC.pbs
qsub IBV_2446_YAM.pbs
