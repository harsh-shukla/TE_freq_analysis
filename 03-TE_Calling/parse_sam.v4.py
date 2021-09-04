import re
import sys
from cigar import Cigar


sam_file = sys.argv[1]
bed_file = sys.argv[2]
strain_name = sys.argv[3]


fbed = open(bed_file,"r")
fsam = open(sam_file,"r")

for line in fbed:

    cols = re.split("\t|\n",line)

    chr_name  = cols[0]
    TE_start  = int(cols[1]) + 500
    TE_end    = int(cols[2]) - 500
    TE_supfam = cols[3]
    TE_fam    = cols[4]
    TE_len    = TE_end - TE_start 


# The below code parses the alignemts in the <sorted> sam file

alignment_dic = {}  # keeps track of alignments

for line in fsam:

    if line[0]=="@":  # Skip the header lines
        continue

    #print line
 
    col1 = re.split("\t|\n",line)

    curr_contig_name = col1[0]
    curr_contig_chr  = col1[2]

    if curr_contig_chr != chr_name:  # dont even bother with alignmets not of same chromosome 
        continue

    align_start = int(col1[3])
    align_end = align_start

    cigar_str = Cigar(col1[5])

    #print(list(cigar_str.items())[0][1])
    #print(list(cigar_str.items()))

    for each_item in cigar_str.items():

        len_aln = int(each_item[0])
        type_aln    = each_item[1]   

        if type_aln == "M":
            align_end = align_end + len_aln

        elif type_aln == "I":
            pass
 
        elif type_aln == "D":
            align_end = align_end + len_aln

        elif type_aln == "N":
            align_end = align_end + len_aln

        elif type_aln == "S":
            pass

        elif type_aln == "H":
            pass

        elif type_aln == "P":
            pass

        elif type_aln == "=":
            align_end = align_end + len_aln

        elif type_aln == "X":
            align_end = align_end + len_aln

        else:
            print "Something has gone wrong in reading CIGAR"

    if align_start >= (TE_start-2000) and align_end <= (TE_end+2000): # This code was added in v3 for considering alignments only in the TE surrouding(+-2500bp) region (From where the reads were extracted)
        cur_aln_ss = []
        cur_aln_ss.append(align_start)
        cur_aln_ss.append(align_end) 
        
        if curr_contig_name not in alignment_dic:
            alignment_dic[curr_contig_name] = []
            alignment_dic[curr_contig_name].append(cur_aln_ss)

        else:
            alignment_dic[curr_contig_name].append(cur_aln_ss)

#print alignment_dic
#print chr_name,TE_start,TE_end

######################  Make a call if TE is present absent or cannot be determined

present = -1  # 0 for absent , 1 for present and -1 for cannot determined,  2 for conflict

fout_name  = strain_name + ".TE.calls.txt"
ferr1_name = strain_name + ".err.1.txt"
ferr2_name = strain_name + ".err.2.txt"

fout = open(fout_name,"a")
ferr1 = open(ferr1_name,"a")   # stores my sanity check and doublet error
ferr2 = open(ferr2_name,"a")   # stores record having conflict

flagged = 0

num_of_contigs = len(alignment_dic)

### Make a call if the TE is absent  --- Note: This code might fail if the deletion size reported in a single alignemnt is comparable to size of TE  ###

for curr_contig in  alignment_dic:

    num_of_alns = len(alignment_dic[curr_contig])
          
    if num_of_alns >= 2: # check if it has atleast 2 alignments in total

        conditions_true = 0 
        alignment_vec=[]

        for each_aln in alignment_dic[curr_contig]:

            curr_aln_start = each_aln[0]
            curr_aln_end   = each_aln[1]
            curr_aln_len   = curr_aln_end - curr_aln_start

            alignment_vec.append(curr_aln_start)
            alignment_vec.append(curr_aln_end)
 
            if (curr_aln_end >= (TE_start-20) and curr_aln_end <= (TE_start+20)) and curr_aln_len >=50:  # check if one of alignments end within 20bp of TE start and that alignemt is > 50 bp
                conditions_true = conditions_true + 1
           
            if (curr_aln_start >= (TE_end-20) and curr_aln_start <= (TE_end+20)) and curr_aln_len >=50:  #  check if one of alignments start within 20bp of TE end and that alignemt is > 50 bp
                conditions_true = conditions_true + 1


        if conditions_true == 2:
            present = 0   # The TE is definitely absent

        else:             # This code was added in v2 to cover cases when the contig alignment went a little off or strain specifc alignment/short read mapping issues but there is a clear signal of a deletion being present
            alignment_vec.sort()
            if len(alignment_vec) == 4:

                possible_del_len =  alignment_vec[2] - alignment_vec[1]            
                TE_len_err  = (TE_len*10)/100

                if (possible_del_len > (TE_len - TE_len_err)) and (possible_del_len < (TE_len + TE_len_err)) and (alignment_vec[1] <= (TE_start+20) and alignment_vec[2] >= (TE_end-20)):  # this will check if possible deletion is with +- 10% of TE length - the aligmnent vec condition was added in v3 this makes sure that alignment is not going into TE much (20bp error range) to account for cases when phrap aggresivly assembles into a single contig even if the contig is present (to avoid mis labeling TEs))
                    present=0

        if conditions_true == 1:
            flagged = 1   # this might also catch the doublets (record)

        del alignment_vec


if (num_of_contigs == 1 and num_of_alns == 1) and present==-1:  # this is my sanity check to see - this might also catch the doublets (record)
    flagged = flagged + 1

#print present

### Make a call if the TE is present   ###

spanning = 0
within   = 0


    
for curr_contig in  alignment_dic:

    num_of_alns = len(alignment_dic[curr_contig])

    if num_of_alns == 1:   # I am assuming here there is a single alignemnt because i am looking for only two contigs spanning the boundary (or a single contig spanning the entire TE).

        alignment_vec=[] 
        for each_aln in alignment_dic[curr_contig]:
            
            #print each_aln 
            curr_aln_start = each_aln[0]
            curr_aln_end   = each_aln[1]
            curr_aln_len   = curr_aln_end - curr_aln_start

            cond1 = (curr_aln_end >= (TE_start + 20)) and (curr_aln_start <= (TE_start - 30)) # 30bp to 20bp changed in the v4 according to observations
            cond2 = (curr_aln_start <= (TE_end - 20)) and (curr_aln_end >= (TE_end + 30))     # 30bp to 20bo changed in the v4 according to observations   

            if (cond1) and (not(cond2)):       # The contig is spanning only the start insertion site of the TE and not the ending
                spanning = spanning + 1

            elif (not(cond1)) and (cond2):     # The contig is spanning only the end   insertion site of the TE and not the starting 
                spanning = spanning + 1

            elif (cond1) and (cond2):          # The same contig is spanning the entire TE including the insertion site
                spanning = spanning + 2

            elif (curr_aln_start > TE_start) and (curr_aln_end < TE_end):
                within = within + 1

            else:
                pass   # this contig entirely lies outside TE region probably

    if num_of_alns == 2: # This code was added in v4 to account for cases when phrap aggressively assembles two seperate contigs (supposingly each one spanning the TE boundary) into one instead of two seperate Contigs and bwa reports one as primary and other as supplementary alignment beacasue of the gap of bases in TE

        alignment_vec=[]                          
        for each_aln in alignment_dic[curr_contig]:# Above Comment continued : This is technically phrap issue but since i dont have much time to optimize phrap parametesrs i am just gonna add this piece of code here and hopefully it works out

            curr_aln_start = each_aln[0]
            curr_aln_end   = each_aln[1]
            curr_aln_len   = curr_aln_end - curr_aln_start

            alignment_vec.append(curr_aln_start)
            alignment_vec.append(curr_aln_end)

        alignment_vec.sort()
        if len(alignment_vec) == 4:

            first_start  = alignment_vec[0]
            first_end    = alignment_vec[1]
            second_start = alignment_vec[2]
            second_end   = alignment_vec[3]

            if (first_start <= (TE_start-30) and first_end >= (TE_start+50)) and (second_start <= (TE_end - 50) and second_end >= (TE_end+30)): # Boundaries were made more stringent of 50 bp inside just to be sure
                spanning = spanning + 2

        del alignment_vec



if (spanning==2) or (spanning==1 and within>=1):

    if present==0:
        present=2
    else:
        present=1


str_print_rec = chr_name + "\t" + str(TE_start) + "\t" + str(TE_end) + "\t" + TE_supfam + "\t" + TE_fam + "\t" + str(TE_len) + "\t" + str(present) + "\t" + str(num_of_contigs) + "\t" + str(spanning) + "\t" + str(within)

fout.write(str_print_rec)
fout.write("\n")


if flagged==1 and present==-1:
    ferr1.write(str_print_rec)
    ferr1.write("\n")

if present==2:
    ferr2.write(str_print_rec)
    ferr2.write("\n")

fout.close()
ferr1.close()
ferr2.close()

