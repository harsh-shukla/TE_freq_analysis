#!/bin/bash


## Specify the job name
#SBATCH --job-name=TE_Calling

## account to charge
#SBATCH -A GRYLEE_LAB

# Define nodes tasks and cpu_per_task
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=1    ## number of cores the job needs

## Specify which queues you want to submit the jobs too
#SBATCH -p standard

# Specify where standard output and error are stored
#SBATCH --error=callTE_Annv3_HPC3.log.err
#SBATCH --output=callTE_Annv3_HPC3.log.out

# Pass the current environment variables
#SBATCH --export=ALL

# Go to current working directory 
#SBATCH --chdir=.

## LOAD MODULES or ENVIRONMENTS  ##
source ~/miniconda3/etc/profile.d/conda.sh
conda activate TE_grace

## ALL THE INPUTs REQUIRED FOR RUNNING THE PIPELINE ##

STRAIN=$1
#STRAIN=MD105

REF=/dfs3b/grylee/hgshukla/D_MEL/Ref_annotations/W1118/w1118_scaffold.fasta
TE_FILE=/dfs3b/grylee/hgshukla/D_MEL/Illumia_data/SVTOOLS/Dmel_w1118_updated_TE_annotation.txt

GENERATE_INP_PHRAP=generate_input_phrap.py
PARSE_SAM=parse_sam.v4.py

BAM_PATH=`pwd`/$STRAIN.sorted.merged.50q.bam

# Get the sra ids
declare -a SRA_IDS
i=0

for f in *.sra
do
	#echo "$f"
        IFS='.' read -ra SRAN <<< "$f"
        SRA_IDS[$i]=${SRAN[0]}

        i=`expr $i + 1`

done

# Get the fastq files name
j=0
for CURR_SRA in ${SRA_IDS[@]}
do

    cd $CURR_SRA\_F

    #FASTQ_R1="$(echo $CURR_SRA\_1.fastq.gz)"
    #FASTQ_R2="$(echo $CURR_SRA\_2.fastq.gz)"
    #echo $FASTQ_R1
    #echo $FASTQ_R2
    #echo $NAME
 
    ## Code to handle trimming

    TRMDIR="Trimmed"
    if [ -d "$TRMDIR" ]; then
         cd $TRMDIR
         FASTQ_R1="$(echo $TRMDIR/$CURR_SRA\_1_val_1.fq.gz)"
         FASTQ_R2="$(echo $TRMDIR/$CURR_SRA\_2_val_2.fq.gz)"
         cd ..
    else
         FASTQ_R1="$(echo $CURR_SRA\_1.fastq.gz)"
         FASTQ_R2="$(echo $CURR_SRA\_2.fastq.gz)"
    fi    


    CURR_FASTQ_PATH_1=`pwd`/$FASTQ_R1
    CURR_FASTQ_PATH_2=`pwd`/$FASTQ_R2    

    FASTQ_FILES[$j]=$CURR_FASTQ_PATH_1
    j=`expr $j + 1`
    FASTQ_FILES[$j]=$CURR_FASTQ_PATH_2
    j=`expr $j + 1`


    cd ..

done


FASTQ_A_1=${FASTQ_FILES[0]}
FASTQ_A_2=${FASTQ_FILES[1]}

#CHR_IN=$1
mkdir Run_pipeline
cd Run_pipeline

### INITIALIZE THE RUN ##

rm Run.Finished      # Delete the older log
rm -r Contigs
rm -r SAMs

mkdir Contigs   # Directory to store all the contigs 
mkdir SAMs      # Directory to store all the alignemnts of contig


echo CHR$'\t'START$'\t'END$'\t'SUP_FAM$'\t'FAM$'\t'LEN$'\t'PRESENT$'\t'NUM_CONTIGS$'\t'SPANNING$'\t'WITHIN > $STRAIN.TE.calls.txt
echo CHR$'\t'START$'\t'END$'\t'SUP_FAM$'\t'FAM$'\t'LEN$'\t'PRESENT$'\t'NUM_CONTIGS$'\t'SPANNING$'\t'WITHIN > $STRAIN.err.1.txt
echo CHR$'\t'START$'\t'END$'\t'SUP_FAM$'\t'FAM$'\t'LEN$'\t'PRESENT$'\t'NUM_CONTIGS$'\t'SPANNING$'\t'WITHIN > $STRAIN.err.2.txt


### PROCESS EACH TE/RECORD IN THE THE FILE ###


skip_header=1

while IFS= read -r line
do

    if [ $skip_header -eq 1 ]    # Skip the header
    then        
        skip_header=0
        continue 
    fi

    # Parse the TE file to extract values required for generating bed file

    CHR="$(echo $line | cut -d' ' -f2)"

    #if [ "$CHR" != "$CHR_IN" ]
    #then
    #    continue
    #fi
    
    TE_START="$(echo $line | cut -d' ' -f3)"
    TE_END="$(echo $line | cut -d' ' -f4)"
    SUP_FAM="$(echo $line | cut -d' ' -f5)"
    FAM="$(echo $line | cut -d' ' -f1)"
    TE_id="$(echo $CHR\_$TE_START\_$TE_END)"
  
    #echo $line
    #echo $CHR
    #echo $TE_START
    #echo $TE_END
    #echo $TE_id

    # Generate the bed files and get the (unique) reads mapping within -500 to +500 of the TE [region from (TE start-500)  to  (TE_end + 500)]

    BED_START=`expr $TE_START - 500`
    BED_END=`expr $TE_END + 500`
    echo $CHR$'\t'$BED_START$'\t'$BED_END$'\t'$SUP_FAM$'\t'$FAM > TE_current.bed
    REGION="$(echo $CHR:$BED_START-$BED_END)"

    #bedtools intersect -abam $BAM_PATH -b TE_current.bed > TE_current.bam
    samtools view -h $BAM_PATH $REGION > TE_current.bam
    samtools view TE_current.bam | awk -F"\t" '{print $1}' | sort | uniq > TE_current.reads.UNIQ.list

    # Get the reads from the all fastq files # Rate limiting step no way to index fastq files and cannot extarct sequences from SAM beacuse of hardclip

    seqtk subseq $FASTQ_A_1 TE_current.reads.UNIQ.list > TE_current_1.fastq
    seqtk subseq $FASTQ_A_2 TE_current.reads.UNIQ.list > TE_current_2.fastq

    # assemble the reads using phrap 
     
    python $GENERATE_INP_PHRAP TE_current_1.fastq TE_current_1.fastq # This script will convert the fastq files into format phrap can understand (TE_current.fasta and TE_current.qual)

    mv TE_current.fasta $TE_id.fasta
    mv TE_current.fasta.qual $TE_id.fasta.qual
    phrap $TE_id.fasta -vector_bound 0 -forcelevel 5 -minscore 30 -minmatch 10 > TE_current.phrap.log  # Contigs in $TE_id.fasta.contigs

    /data/homezvol1/hgshukla/Softwares/bwa/bwa-0.7.17/bwa mem $REF $TE_id.fasta.contigs | samtools sort - -o $TE_id.contigs_align.sorted.sam

    python $PARSE_SAM $TE_id.contigs_align.sorted.sam TE_current.bed $STRAIN

    mv $TE_id.fasta.contigs ./Contigs/
    mv $TE_id.contigs_align.sorted.sam ./SAMs/

    rm TE_current*
    rm $TE_id*

    echo $TE_id >> Run.Finished

done < $TE_FILE

#mv $STRAIN.TE.calls.txt $STRAIN.$CHR_IN.TE.calls.txt
