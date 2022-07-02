#!/usr/bin/python/3.8.2

######################################################
## filter_blast.py
## run after blast.sh
## authors: Georgina Samaha and Mitchell O'Brien
## version 1
######################################################


# run as filter_blast.py <blastout>

# Filter BLAST output
# Based on E value (needs input threshold), determine SNP position relative to flank side and length (taken from header)
# Exclude alignments that fell short of the SNP
# Exclude alignments that fell short of E threshold
# Exclude alignments that failed to return any result
# Print to stdout: total alignments, how many were successful, how many failed each category

import os
import sys
import re
import pandas as pd

if sys.version_info[0] < 3:
        raise Exception("Must be using Python 3")
else:
        print(sys.version)

max_e=1e-10
ref=sys.argv[1].split('.')[0]
print(ref)
outDir="./Output"

blast=pd.read_csv(sys.argv[1], sep='\t')
probes=pd.read_csv('Output/probes.seq', sep='\t', header=None)
#new=open('%s_newpos.txt’ % ref, "wb")
new=open('%s_newpos.txt'% ref, "w+")
unmap=open('%s_failedBlast.txt'% ref, "w+")
multimap=open('%s_multipleBLAST.txt'% ref, "w+")

print("Filtering blast results ...\n\n")

#create dict to save new snp positions
blasthash = {}
#define count variables
total=0
short=0
good=0
flank_error=0
multi=0
filt_e=0

#create header
unmap.write("SNP_ID\tmappedTo\t%ID\tAlnLen\tMis Gap\tqStart\tqEnd\tsStart\tsEnd evalue\tbit_score\treason\n")
multimap.write("SNP_ID\tmappedTo\t%ID\tAlnLen\tMis\tGap\tqStart\tqEnd\tsStart\tsEnd\tevalue\tbit_score\n")

for line, value in blast['evalue'].iteritems():
        if value < max_e:
                total+=1
                SNP=blast.loc[line,'#SNP_ID']
                id_info=SNP.split("_")
                flank=id_info[-2]
                length=int(id_info[-1])
                id=SNP.replace("_"+flank+"_"+str(length),"")
                if flank == "left":
                        if length == blast.loc[line,'qEnd']:
                                if blast.loc[line,'sStart'] == blast.loc[line,'sEnd']:
                                        snp_pos=int(blast.loc[line,'sEnd']+1)
                                else:
                                        snp_pos=int(blast.loc[line,'sEnd']-1)
                        else:
                                #print("Alignment fell short of SNP: "+str(id))
                                result=(str(id)+'\t'+"NULL"+'\t'+"NULL"+'\t'+str(1))
                                blasthash[id]=result
                                #print(blasthash[id])
                                short+=1
                                unmap.write(blast.loc[[line]].to_string(index=False,header=False)+"\tfell short"+'\n')
                elif flank == "right":
                        if length == blast.loc[line,'qStart']:
                                if blast.loc[line,'sStart'] == blast.loc[line,'sEnd']:
                                        snp_pos=int(blast.loc[line,'sStart']-1)
                                else:
                                        snp_pos=int(blast.loc[line,'sStart']+1)
                        else:
                                #print("Alignment fell short of SNP: "+str(id))
                                result=(str(id)+'\t'+"NULL"+'\t'+"NULL"+'\t'+str(1))
                                blasthash[id]=result
                                #print(blasthash[id])
                                short+=1
                                unmap.write(blast.loc[[line]].to_string(index=False,header=False)+"\tfell short"+'\n')
                else:
                        #print("Unrecognised flank: "+str(line))
                        result=(str(id)+'\t'+"NULL"+'\t'+"NULL"+'\t'+str(1))
                        blasthash[id]=result
                        flank_error+=1
                        unmap.write(blast.loc[[line]].to_string(index=False,header=False)+"\tunrecognised flank"+'\n')
                if snp_pos>0:
                        if id in blasthash:
                                #print(str(id)+" mapped to multiple locations picking best evalue")
                                multi+=1
                                currenthit=blasthash[id].split("\t")
                                current_eval=currenthit[3]
                                multimap.write(blast.loc[[line]].to_string(index=False,header=False)+'\n')
                                if blast.loc[line,'evalue'] < float(current_eval):
                                        chrom=str(blast.loc[line,'mappedTo'])
                                        evalue=str(blast.loc[line,'evalue'])
                                        result=(str(id)+'\t'+chrom+'\t'+str(snp_pos)+'\t'+str(evalue))
                                        blasthash[id]=result
                        else:
                                chrom=str(blast.loc[line,'mappedTo'])
                                evalue=str(blast.loc[line,'evalue'])
                                result=(str(id)+'\t'+chrom+'\t'+str(snp_pos)+'\t'+str(evalue))
                                #print(result)
                                blasthash[id]=result
                                good+=1
        else:
                total+=1
                print("No significant hit for "+str(id))
                result=(str(id)+'\t'+"NULL"+'\t'+"NULL"+'\t'+str(0))
                blasthash[id]=result
                filt_e+=1
                unmap.write(blast.loc[[line]].to_string(index=False,header=False)+"\tfailed e"+'\n')

print("Alignments in BLAST output: "+str(total)+"\n"+"Failed E value cut-off: "+str(filt_e)+"\n"+"Alignment missed SNP: "+str(short)+"\n"+"Error reading flank designation (should be zero): "+str(flank_error)+"\n"+str(good)+" alignments passed filtering"+"\n")
tot=filt_e+short+flank_error+good
if tot != total:
        print("WARNING: "+str(total)+" SNPs in BLAST output, "+str(tot)+" SNPs accounted for during filtering\n")

found=0
missed=0
tot=0
for line, snp in probes[0].iteritems():
        tot+=1
        if snp in blasthash:
                new_pos=blasthash[snp].split('\t')
                if str(new_pos[1]).startswith("chrUn"):
                        new_pos[1]=0
                result=(str(new_pos[0])+'\t'+str(new_pos[1])+'\t'+str(new_pos[2]))
                found+=1
                new.write(result+'\n')
        else:
                result=(snp+'\t'+"0"+'\t'+"0")
                missed+=1
                new.write(result+'\n')
print("Of "+str(tot)+" input SNPs, "+str(found)+" accounted for in BLAST, "+str(missed)+" failed BLAST altogether.")

print ("canfam4.newpos.txt successfully created\n\n")
