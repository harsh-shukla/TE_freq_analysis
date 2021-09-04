#!/bin/bash

## Specify the job name
#SBATCH --job-name=Trimming

## account to charge
#SBATCH -A GRYLEE_LAB

# Define nodes tasks and cpu_per_task
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=1    ## number of cores the job needs

## Specify which queues you want to submit the jobs too
#SBATCH -p standard

# Specify where standard output and error are stored
#SBATCH --error=trim.err
#SBATCH --output=trim.out

# Pass the current environment variables
#SBATCH --export=ALL

# Go to current working directory 
#SBATCH --chdir=.

######## Load all modules needed ###############
module load cutadapt/2.10
################################################

SRR_ID=$1  # Pass the SRR id here

FASTQ_R1=$SRR_ID\_1.fastq.gz
FASTQ_R2=$SRR_ID\_2.fastq.gz

mkdir Trimmed/
OUTDIR=Trimmed/


/data/homezvol1/hgshukla/Softwares/TrimGalore-0.6.0/trim_galore \
--quality 20 --phred33 --stringency 3 --length 70 --trim-n --max_n 2 --gzip \
--output_dir $OUTDIR \
--paired \
$FASTQ_R1 $FASTQ_R2



######## Unload all modules ###################
module unload cutadapt/2.10
################################################





