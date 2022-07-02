#!/usr/bin/python/3.8.2

#######################################################
## create_probes.py
## authors: Georgina Samaha and Mitchell O'Brien
## run as: python create_probes.py <manifest.csv> --illumina
## version 1
#######################################################

import sys
import os
import re
import pandas as pd

if sys.version_info[0] < 3:
        raise Exception("Must be using Python 3")
else:
        print(sys.version)

print("creating probes.fasta from "+sys.argv[1]+" using "+sys.argv[2]+"...\n\n")

path = "Output"
os.makedirs(path, exist_ok=True)
outfile = open("Output/probes.fasta", "w+")
probeseq = open("Output/probes.seq", "w+")

if sys.argv[2] == "--illumina":

        manifest = pd.read_csv(sys.argv[1],skiprows=7)

        for probes, seq in manifest['SourceSeq'].iteritems():
                if pd.notnull(seq):
                        snpid = manifest.loc[probes,'Name']
                        seq_split = re.split('\\[|\\]',seq)
                        left=seq_split[0]
                        right=seq_split[2]
                        left_length = len(left)
                        right_length = len(right)
                        if left_length >= right_length:
                                outfile.write(">"+snpid+"_left_"+str(left_length)+"\n"+left+"\n")
                                probeseq.write(snpid+"\t"+seq+"\n")
                        else:
                                outfile.write(">"+snpid+"_right_"+str(right_length)+"\n"+right+"\n")
                                probeseq.write(snpid+"\t"+seq+"\n")

elif sys.argv[2] == "--plink":

        manifest = open(sys.argv[1], 'r')

        #def fasta2dict(reffasta):
        #       """
        #       Read a fasta file and create dictionary(name,value)
        #       where name will be the unique chrom/scaff identifier in header
        #       and the value is the corresponding sequence
        #       """
        #       fastadic = {}
        #       chr_seq = []
        #       fasta = open(reffasta)
        #       for line in fasta:
        #               if line.startswith(">"):
        #                       chr_id = line.split(' ')[0].replace('>','')
        #                       chr_seq = ""
        #               else:
        #                       chr_seq += line.rstrip()
        #                       fastadic[chr_id] = ''.join(chr_seq)
        #               return fastadic

        fasta = {}
        chr_seq = []
        reffasta = open(sys.argv[3])
        for line in reffasta:
                if line.startswith(">"):
                        chr_id = line.split(' ')[0].replace('>','')
                        chr_seq = ""
                else:
                        chr_seq += line.rstrip()
                        fasta[chr_id] = ''.join(chr_seq)

        #fasta = fasta2dict(str(sys.argv[3]))
        flanksize = 60

        for line in manifest:
                # print(line)
                fields = line.split("\t")
                snpid = fields[1]
                snp_chr = fields[0]
                #print(snp_chr)
                snp_pos = int(fields[3])-1
                # print(snp_pos)
                ref_al = fields[4]
                # print(ref_al)
                alt_al = fields[5]
                # print (alt_al)
                if snp_chr not in fasta:
                        left = str("N"*60)
                        right = str("N"*60)
                else:
                        left = fasta[snp_chr][snp_pos-flanksize:snp_pos]
                        right = fasta[snp_chr][snp_pos+1:snp_pos+flanksize+1]
                        outfile.write(">"+snpid+"_left_"+flanksize+"\n"+left+"\n")
                # print(">"+snp_chr+"_"+str(snp_pos+1)+"\t"+left+"["+ref_al+"/"+alt_al+"]"+right)
                # outfile.write(">"+snpid+"_left_"+flanksize+"\n"+left+"\n")
                probeseq.write(snpid+"\t"+left+"["+ref_al+"/"+alt_al+"]"+right+"\n")

print ("probes.fasta and probes.seq successfully created\n\n")
