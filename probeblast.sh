#!/bin/bash

#######################################################
## probeblast.sh
## run after create_probes.py
## authors: Georgina Samaha and Mitchell O'Brien
## run as: bash probeblast.sh <reference.fasta>
## version 1
#######################################################

#PBS -P Georgie
#PBS -l select=1:ncpus=2:mem=10GB
#PBS -l walltime=01:00:00
#PBS -N blast
#PBS -m e
#PBS -M georgina.samaha@sydney.edu.au
#PBS -q defaultQ
#PBS -W umask=022

# written for USyd Artemis HPC which uses PBSpro job scheduler
# run as bash blast.sh $ref

probeseq=Output/probes.fasta
ref=$1 #../../Ref/canfam4.fasta

# run as an interactive session on Artemis
# when writing pbs profile config make sure to qsub

module load blast+/2.9.0

NCPUS=24

#blast+ to convert fasta to indexed and searchable version of the same information .nin, .nsq and .nhr
makeblastdb -in ${ref} \
       -dbtype nucl \
       -blastdb_version 5 \
       -parse_seqids \
       -out ${ref}.db

#blastn to search subject fasta with query sequence(s) and return filtered hits
blastn -query ${probeseq} \
       -db ${ref}.db \
       -evalue 1e-20 \
       -outfmt 7 \
       -num_threads ${NCPUS} \
       -out ${ref}.BLASTn

#extract blast output to readable format
echo "#SNP_ID   mappedTo        %ID     AlnLen  Mis     Gap     qStart  qEnd    sStartsEnd     evalue  bit_score" >${ref}.BLASTn_Output
grep -v '#' ${ref}.BLASTn >> ${ref}.BLASTn_Output

mv ${ref}.BLASTn ./Output
mv ${ref}.BLASTn_Output ./Output
