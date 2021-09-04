import sys
import re


updated_TE_ch  = open(sys.argv[1],"r")
orignal_TE     = open(sys.argv[2],"r")

orignal = []
updated = []

for line in updated_TE_ch: # Skip the header
    break

for line in updated_TE_ch:
     col1 = re.split("\n",line)
     updated.append(col1[0])

#print updated

for line in orignal_TE: # Skip the header
    break

for line in orignal_TE:
     col1 = re.split("\n",line)
     orignal.append(col1[0])

#print orignal[0]

upd_len = len(updated)
org_len = len(orignal)

not_updated     = 0
updated_t       = 0
sig_updated     = 0
not_sig_updated = 0

if upd_len == org_len:

    for i in range(0,org_len):

        #org_rec = orignal[i]
        #upd_rec = updated[i]

        col1 = re.split("\t",orignal[i])
        orig_start = int(col1[2])
        orig_end   = int(col1[3])

        col2 = re.split("\t",updated[i])
        upd_start = int(col2[2])
        upd_end   = int(col2[3])
        changed   = int(col2[7])


        if changed == 0:   
            not_updated = not_updated + 1  # Not updated

        else:              
            updated_t = updated_t + 1         # Updated

            cond1 = (upd_start >= (orig_start-20)) and (upd_start <= (orig_start+20))
            cond2 = (upd_end   >= (orig_end  -20)) and (upd_end   <= (upd_end   +20))

            if (cond1) and (cond2):
                not_sig_updated = not_sig_updated + 1   # Not significantly updated
            else:
                sig_updated = sig_updated + 1 # Boundaries significantly updated 



print "TOTAL NUMBER OF TE RECS           : ",org_len
print "\n"
print "TOTAL NUMBER OF TE UPDATED        : ",updated_t
print "      SIGNIFCANTLY UPDATED        : ",sig_updated
print "  NOT-SIGNIFCANTLY UPDATED        : ",not_sig_updated 
print "\n"
print "TOTAL NUMBER OF TE NOT UPDATED    : ",not_updated







