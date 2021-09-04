#!/bin/bash

## Specify the job name
#SBATCH --job-name=Download_Process

## account to charge
#SBATCH -A GRYLEE_LAB

# Define nodes tasks and cpu_per_task
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=1    ## number of cores the job needs

## Specify which queues you want to submit the jobs too
#SBATCH -p standard

# Specify where standard output and error are stored
#SBATCH --error=download_SRR_process.log.err
#SBATCH --output=download_SRR_process.log.out

# Pass the current environment variables
#SBATCH --export=ALL

# Go to current working directory 
#SBATCH --chdir=.


while IFS=, read -r SRR_ID STRAIN
do


    if [ -d "$STRAIN" ] 
    then
        #echo "Directory /path/to/dir exists."
        cd $STRAIN

        prefetch $SRR_ID -O .
        fastq-dump --outdir $SRR_ID\_F --gzip --split-files $SRR_ID.sra 

    else
        #echo "Error: Directory /path/to/dir does not exists."
        mkdir $STRAIN
        cd $STRAIN         
 
        prefetch $SRR_ID -O .
        fastq-dump --outdir $SRR_ID\_F --gzip --split-files $SRR_ID.sra

    fi

    cd ..

done < $1
