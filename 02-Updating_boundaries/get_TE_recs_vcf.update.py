import sys
import re


vcffile = open(sys.argv[1],"r")
TEfile  = open(sys.argv[2],"r")

fout = open("DMEL_TE_Annotation.UPDATED.txt","w")

vcf_recs = {}

for line in vcffile:


    if line[0:2]=="##" or line[0:6]=="#CHROM":
        continue


    cols = re.split("\t|\n",line)

    curr_chr = cols[0]

    if curr_chr not in vcf_recs:
        vcf_recs[curr_chr] = []
        vcf_recs[curr_chr].append(line)
    else:
        vcf_recs[curr_chr].append(line)


for line in TEfile: ## Skip the header
    fout.write(line)
    break


for line in TEfile:

    col1 = re.split("\t|\n",line)

    TE_chr   = col1[1]
    TE_start = int(col1[2])
    TE_end   = int(col1[3])

    TE_len = TE_end - TE_start

    sup_fam = col1[0]
    fam     = col1[4]
    TE_col6    = col1[5]
    TE_col7    = col1[6]

    updated_TE_start = TE_start
    updated_TE_end   = TE_end

    print TE_chr,"\t",TE_start,"\t",TE_end,"\t",TE_len 
    changed=0

    if TE_len > 750:
      
        possible_rec = []
        curr_rec=""

        for every_rec in vcf_recs[TE_chr]:

            col2 = re.split("\t|\n",every_rec)
            curr_vcfrec_start = int(col2[1])

            col3 = re.split(";",col2[7])
            col4 = re.split("=",col3[2]) 
            curr_vcfrec_end   = int(col4[1])
            
            svlen = curr_vcfrec_end - curr_vcfrec_start
                     

            if ( TE_start>(curr_vcfrec_start-1000) and TE_start<(curr_vcfrec_start+1000) ) and ( TE_end>(curr_vcfrec_end-1000) and TE_end<(curr_vcfrec_end+1000) ):
                 if svlen > 750:   # this will get rid of any small svs that got captured while filtering in +- 1000 bp regions
                     possible_rec.append(every_rec)
                     
           
        if len(possible_rec)==1:

            col2 = re.split("\t|\n",possible_rec[0])
            curr_vcfrec_start = int(col2[1])

            col3 = re.split(";",col2[7])
            col4 = re.split("=",col3[2]) 
            curr_vcfrec_end   = int(col4[1])

            poss_imprecise = col3[4]

            if poss_imprecise == "IMPRECISE":
                col6 = re.split("=",col3[9])
                curr_SU = int(col6[1])
            else:
                col6 = re.split("=",col3[8])
                curr_SU = int(col6[1])

            if curr_SU > 500:   # Minmum SU for updating even for single records
                curr_rec=possible_rec[0]            
                updated_TE_start =  curr_vcfrec_start
                updated_TE_end   =  curr_vcfrec_end
                changed=1

        elif len(possible_rec) > 1:

            SU=500   # Minmum SU for updating in case of multiple records
                 
            for every_rec in possible_rec:

                col2 = re.split("\t|\n",every_rec)
                curr_vcfrec_start = int(col2[1])

                col3 = re.split(";",col2[7])
                col4 = re.split("=",col3[2]) 
                curr_vcfrec_end   = int(col4[1])

                poss_imprecise = col3[4]

                if poss_imprecise == "IMPRECISE":
                    col6 = re.split("=",col3[9])
                    curr_SU = int(col6[1])
                else:
                    col6 = re.split("=",col3[8])
                    curr_SU = int(col6[1])


                if curr_SU > SU:   # Get the record with highest SU genrally thats the one that is proper
                    curr_rec = "CHOICE:"  + every_rec
                    SU       = curr_SU
                    updated_TE_start =  curr_vcfrec_start
                    updated_TE_end   =  curr_vcfrec_end                        
                    changed=1

            if  changed==0:
                pass                                       
                print "FLAG1 : There were more than 1 records found having svlen > 750 but none of them met the condtion of SU > 500" # 

        else:
            pass 
            print "FLAG2 : There were more than 1 record found but none of them met the condtion of svlen > 750" # 
                         
    
        print curr_rec
        del possible_rec
        del curr_rec


    else:  # TE < 750bp

        SU=800   # Minmum SU for updating any small SV record
        curr_rec=""
        num_of_recs=0

        for every_rec in vcf_recs[TE_chr]:

            col2 = re.split("\t|\n",every_rec)
            curr_vcfrec_start = int(col2[1])

            col3 = re.split(";",col2[7])
            col4 = re.split("=",col3[2]) 
            curr_vcfrec_end   = int(col4[1])
            
            svlen = curr_vcfrec_end - curr_vcfrec_start
                     
            poss_imprecise = col3[4]

            if poss_imprecise == "IMPRECISE":
                col6 = re.split("=",col3[9])
                curr_SU = int(col6[1])
            else:
                col6 = re.split("=",col3[8])
                curr_SU = int(col6[1])
       
            
            if ( TE_start>(curr_vcfrec_start-1000) and TE_start<(curr_vcfrec_start+1000) ) and ( TE_end>(curr_vcfrec_end-1000) and TE_end<(curr_vcfrec_end+1000) ):
                num_of_recs = num_of_recs + 1

                if curr_SU > SU:   # Get the record with highest SU genrally thats the one that is proper and here minimum SU is 800
                    curr_rec = "SMALL:"  + every_rec
                    SU       = curr_SU
                    updated_TE_start =  curr_vcfrec_start
                    updated_TE_end   =  curr_vcfrec_end                        
                    changed=1
 
        if changed==0 and num_of_recs>0:
            pass                                       
            print "FLAG3 : There were more than 1 record found for this small TE but none of them met the condtion of SU > 800"

        print curr_rec
        del curr_rec




    print "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"


    str_print = sup_fam + "\t" + TE_chr + "\t" + str(updated_TE_start) + "\t" + str(updated_TE_end) + "\t" + fam + "\t" + TE_col6 + "\t" + TE_col7
  
    fout.write(str_print)
    fout.write("\n")
