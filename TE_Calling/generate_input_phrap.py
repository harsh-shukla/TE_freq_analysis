import re
import sys

col1 = re.split('\.|_',sys.argv[1])
fasta_name = col1[0] + "_" + col1[1] + ".fasta" 
qual_name  = col1[0] + "_" + col1[1] + ".fasta.qual"

f_out = open(fasta_name,"w")
q_out = open(qual_name,"w")


f1 = open(sys.argv[1],"r")
line_no=1

for line in f1:

    if line_no%4 == 1:

        col2 = re.split("@| |\.",line)
        curr_rname = ">"+ col2[1] + "_" + col2[2] + ".z"
        #print curr_rname
        #break
        f_out.write(curr_rname)
        f_out.write("\n")
        q_out.write(curr_rname)
        q_out.write("\n")


    elif line_no%4 == 2:
        f_out.write(line)

    elif line_no%4 == 0:
        col3 = re.split("\n",line)
        curr_qual = col3[0]
        qual_str = "" 

        for each_val in curr_qual:           
            conv_val = ord(each_val)-33
            #print(chr(conv_val))
            qual_str = qual_str + str(conv_val) + " "
        
        q_out.write(qual_str)
        q_out.write("\n")

    else:
        pass

    line_no = line_no + 1


f2 = open(sys.argv[2],"r")
line_no=1

for line in f2:

    if line_no%4 == 1:

        col2 = re.split("@| |\.",line)
        curr_rname = ">"+ col2[1] + "_" + col2[2]  + ".y"
        #print curr_rname
        #break
        f_out.write(curr_rname)
        f_out.write("\n")
        q_out.write(curr_rname)
        q_out.write("\n")


    elif line_no%4 == 2:
        f_out.write(line)

    elif line_no%4 == 0:
        col3 = re.split("\n",line)
        curr_qual = col3[0]
        qual_str = "" 

        for each_val in curr_qual:           
            conv_val = ord(each_val)-33
            #print(chr(conv_val))
            qual_str = qual_str + str(conv_val) + " "
        
        q_out.write(qual_str)
        q_out.write("\n")

    else:
        pass

    line_no = line_no + 1


f_out.close()
q_out.close()
    
