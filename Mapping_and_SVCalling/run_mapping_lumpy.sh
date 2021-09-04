#!/bin/bash

## Specify the job name
#SBATCH --job-name=Mapping

## account to charge
#SBATCH -A GRYLEE_LAB

# Define nodes tasks and cpu_per_task
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=32   ## number of cores the job needs

## Specify which queues you want to submit the jobs too
#SBATCH -p standard

# Specify where standard output and error are stored
#SBATCH --error=pipeline.err
#SBATCH --output=pipeline.out

# Pass the current environment variables
#SBATCH --export=ALL

# Go to current working directory 
#SBATCH --chdir=.

######## Load all modules needed ###############
source ~/miniconda3/etc/profile.d/conda.sh
conda activate TE_grace
################################################

### Run as ###
# sbatch run_mapping_lumpy.sh <STRAIN_NAME>

THREADS=32
REF=/dfs3b/grylee/hgshukla/D_MEL/Ref_annotations/W1118/w1118_scaffold.fasta
STRAIN=$1

# Index the genome if first run
# /data/homezvol1/hgshukla/Softwares/bwa/bwa-0.7.17/bwa index $REF 

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

# Process each SRA files
for CURR_SRA in ${SRA_IDS[@]}
do

    cd $CURR_SRA\_F

    TRMDIR="Trimmed"
    if [ -d "$TRMDIR" ]; then
         cd $TRMDIR
         FASTQ_R1="$(echo ./$TRMDIR/$CURR_SRA\_1_val_1.fq.gz)"
         FASTQ_R2="$(echo ./$TRMDIR/$CURR_SRA\_2_val_2.fq.gz)"
         cd ..
    else
         FASTQ_R1="$(echo $CURR_SRA\_1.fastq.gz)"
         FASTQ_R2="$(echo $CURR_SRA\_2.fastq.gz)"
    fi    

    NAME="$(echo $CURR_SRA)"

    R_STR="$(echo @RG\\tID:$CURR_SRA\\tSM:$STRAIN\\tLB:$CURR_SRA\\tPL:Illumina)" 

    /data/homezvol1/hgshukla/Softwares/bwa/bwa-0.7.17/bwa mem -R $R_STR -t $THREADS $REF $FASTQ_R1 $FASTQ_R2 | samtools sort -@ $THREADS -m 2G -O bam - > $NAME.sorted.bam
    samtools index -@ $THREADS $NAME.sorted.bam
    samtools flagstat -@ $THREADS  $NAME.sorted.bam > $NAME.flagstat.txt

    CURR_BAM_PATH=./$CURR_SRA\_F/$NAME.sorted.bam     
    CURR_BAI_PATH=./$CURR_SRA\_F/$NAME.sorted.bam.bai 
    CURR_FLST_PATH=./$CURR_SRA\_F/$NAME.flagstat.txt

    cd ..


done

# Merge all the bams (for 2 or more SRR ids per strain)
samtools merge -@ $THREADS -b bam_sorted_list.txt $STRAIN.sorted.merged.bam
samtools index -@ $THREADS $STRAIN.sorted.merged.bam
samtools flagstat -@ $THREADS $STRAIN.sorted.merged.bam > $STRAIN.sorted.merged.flagstat.txt

# Link the only bam generated (for 1 SRR ids per strain)
#ln -s $CURR_BAM_PATH $STRAIN.sorted.merged.bam
#ln -s $CURR_BAI_PATH $STRAIN.sorted.merged.bam.bai
#ln -s $CURR_FLST_PATH $STRAIN.sorted.merged.flagstat.txt

# Filter bam Q50 for calling absence/presence of TE's
samtools view -@ $THREADS -b -q 50 $STRAIN.sorted.merged.bam > $STRAIN.sorted.merged.50q.bam
samtools index -@ $THREADS $STRAIN.sorted.merged.50q.bam

################################ Preprocess files for lumpy #########################################

samtools sort -n -@ $THREADS -m 2G -O sam $STRAIN.sorted.merged.bam | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | samtools sort -@ $THREADS -m 2G -O bam - > $STRAIN.samblaster.bam
samtools index -@ $THREADS $STRAIN.samblaster.bam


# Extract the discordant paired-end alignments.
samtools view -b -F 1294 $STRAIN.samblaster.bam  > $STRAIN.discordants.sorted.bam
samtools index $STRAIN.discordants.sorted.bam


# Extract the split-read alignments
samtools view -h $STRAIN.samblaster.bam \
    | /data/homezvol1/hgshukla/Softwares/LUMPY/lumpy_sv/lumpy-sv-v0.2.13/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools sort -@ $THREADS -O bam - \
    > $STRAIN.splitters.sorted.bam
samtools index $STRAIN.splitters.sorted.bam

### Call SVs using Lumpy

/data/homezvol1/hgshukla/Softwares/LUMPY/lumpy_sv/lumpy-sv-v0.2.13/bin/lumpyexpress \
 -B $STRAIN.samblaster.bam \
 -S $STRAIN.splitters.sorted.bam \
 -D $STRAIN.discordants.sorted.bam \
 -o $STRAIN.vcf \
 -P -v \
 -K /data/homezvol1/hgshukla/Softwares/LUMPY/lumpy_sv/lumpy-sv-v0.2.13/bin/lumpyexpress.config

### Filter the vcf file for deletions we are interested in and genotype using svtyper 

bcftools filter -i'(ALT=="<DEL>" && SVLEN<=-95)' $STRAIN.vcf | grep -v "U_" | grep -v "mitochondrion_genome" > $STRAIN.filtered.vcf # Only extratc records for major chromsomes..
/data/homezvol1/hgshukla/Softwares/svtyper/svtyper-0.0.4/svtyper -B $STRAIN.samblaster.bam -S $STRAIN.splitters.sorted.bam -i $STRAIN.filtered.vcf > $STRAIN.filtered.GT.vcf
bcftools filter -i'QUAL>100' $STRAIN.filtered.GT.vcf > $STRAIN.filtered.GT.QUAL100.vcf

