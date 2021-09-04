#!/bin/bash


###### Generate the sorted merged vcf file ######
find ../* | grep "GT.QUAL100.vcf" > lumpyvcfs_file.txt
svtools lsort -f lumpyvcfs_file.txt > DMEL.sorted.vcf                ## replace DMEL with correspoding species you are analysing
svtools lmerge -i DSIM.sorted.vcf -f 20 > DMEL.sorted.merged.vcf     ## replace DMEL with correspoding species you are analysing


###### Update the TE boundaries ######
# The below python files generates a file "DMEL_TE_Annotation.UPDATED.txt" which has the updated boundaries (wherever approriate);final_Script_all_recs_TE.FINAL.txt is a log like file for manual inspection
python get_TE_recs_vcf.update.py         DMEL.sorted.merged.vcf DMEL_TE_Annotation.txt > final_Script_all_recs_TE.FINAL.txt


###### Check statistics for updated TE records 
python get_TE_recs_vcf.update.CHANGED.py DMEL.sorted.merged.vcf DMEL_TE_Annotation.txt # Outputs a file similar to DMEL_TE_Annotation.UPDATED.txt but has an extra column that represents whether that particular TE was updated significantyly or not
python get_significantaly_updated.py DMEL_TE_Annotation.UPDATED.CHANGED.txt DMEL_TE_Annotation.txt # Calculate the statistics -- see below

#TOTAL NUMBER OF TE RECS           :  1030
#
#TOTAL NUMBER OF TE UPDATED        :  701
#      SIGNIFCANTLY UPDATED        :  35
#  NOT-SIGNIFCANTLY UPDATED        :  666
#
#TOTAL NUMBER OF TE NOT UPDATED    :  329


# here NOT-SIGNIFCANTLY UPDATED is if updated boudaries(both) are within +-20 bp of the orignal boundaries and SIGNIFCANTLY UPDATED is rest of the cases  

