#!/bin/bash

######################################################################
# findMatchesPDB.sh                                                  #
# Author:  Dario Ghersi                                              #
# Version: 20130924                                                  #
# Goal:    extract the chain sequences that are available in         #
#          BioLip, and perform a BLAST search against them for each  #
#          reference sequence                                        #
######################################################################

source ../variables.sh

######################################################################
# CONSTANTS                                                          #
######################################################################

workingDir="."
BLASTDB=${workingDir}/blastdb

######################################################################
# CONCATENATE THE REFERENCE SEQUENCES INTO ONE FILE                  #
######################################################################

echo "1. Concatenating the reference sequences into one file"
echo

if [ ! -e ${workingDir}/allSeqs.fasta ]; then
  allRefSeqs=`ls ../step.01/refSeqs`
  ../code/python/buildSeqDatabase.py ${workingDir}/../step.01/mutDir \
                             ${workingDir}/../step.01/refSeqs \
                             ${workingDir}/allSeqs.fasta
fi

######################################################################
# EXTRACT THE CHAIN SEQUENCES FROM BIOLIP                            #
######################################################################

echo "2. Extracting the chain sequences from BioLip"
echo

biolipDir=${workingDir}/../data/biolip
biolipSeqs=${biolipDir}/biolipSeqs.txt
if [ ! -e $biolipSeqs ]; then
  ../code/python/extractSequencesBioLip.py ${biolipDir}/biolipAnn.txt \
                                   $biolipSeqs
fi

######################################################################
# BUILD THE BLAST DATABASE                                           #
######################################################################

echo "3. Making a BLAST database of the BioLip sequences"

## make a BLAST database with the BioLip sequences
if [ ! -d $BLASTDB ]; then
  mkdir $BLASTDB
fi

if [ ! -e ${BLASTDB}/biolipSeqs.phr ] || \
   [ ! -e ${BLASTDB}/biolipSeqs.pin ] || \
   [ ! -e ${BLASTDB}/biolipSeqs.psq ]; then

  ${BLAST_DIR}/makeblastdb -in ${biolipDir}/biolipSeqs.txt \
                           -dbtype prot -out ${BLASTDB}/biolipSeqs
fi

######################################################################
# RUN THE BLAST SEARCH                                               #
######################################################################

echo "4. Running a BLAST search against the BioLip database..."
echo "   this will take a while"
echo

if [ ! -e refSeqsAgainstBioLip.tab ]; then
  export BLASTDB
  ${BLAST_DIR}/blastp -query allSeqs.fasta -db biolipSeqs \
                      -evalue $BLAST_EVALUE \
                      -outfmt 6 -num_threads 4 \
                      -out refSeqsAgainstBioLip.tab
fi
