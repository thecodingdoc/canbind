#!/bin/bash

######################################################################
# assignBindingToRefSeqs.sh                                          #
# Author:  Dario Ghersi                                              #
# Version: 20130924                                                  #
# Goal:    parse the BLAST search and extract the PDB chains that    #
#          match the reference sequences, then align these against   #
#          the PDB chains using CLUSTAL-OMEGA                        #
######################################################################

source ../variables.sh

######################################################################
# CONSTANTS                                                          #
######################################################################

workingDir="."

BLAST_OUTPUT="${workingDir}/../step.02/refSeqsAgainstBioLip.tab"
HIT_SEQS="${workingDir}/../data/biolip/biolipSeqs.txt"
QUERY_SEQS="${workingDir}/../step.02/allSeqs.fasta"

######################################################################
# PARSE THE BLAST FILE AND WRITE THE RESULTS TO FILE                 #
######################################################################

echo "1. Parsing the blast output"
echo

if [ ! -e ${BLAST_OUTPUT}.parsed ]; then
  ../code/python/parseBLASTTAB.py $BLAST_OUTPUT $QUERY_SEQS $HIT_SEQS \
                                  $BLAST_EVALUE ${BLAST_OUTPUT}.parsed
fi

######################################################################
# ALIGN THE REFERENCE SEQUENCES ONTO THE BIOLIP CHAINS USING         #
# CLUSTAL-OMEGA                                                      #
######################################################################

echo "2. Aligning the COSMIC sequences onto the BioLip chains"
echo

alnDir=${workingDir}/refSeqAlnBioLip
if [ ! -d $alnDir ]; then
  mkdir $alnDir
fi

../code/python/alignRefSeq2BioLip.py ${BLAST_OUTPUT}.parsed \
                                     $QUERY_SEQS $HIT_SEQS $CLUSTALO \
                                     $SEQ_IDENTITY $HIT_COVERAGE \
                                     $alnDir

######################################################################
# ASSIGN BINDING INFORMATION TO THE REFERENCE SEQUENCES              #
######################################################################

echo "3. Assigning binding information to the reference sequences"
echo

bindingDir=${workingDir}/bindingInfoBioLip_ID${SEQ_IDENTITY}_COV${HIT_COVERAGE}_BINDID${BINDING_SEQ_ID}
if [ ! -d $bindingDir ]; then
  mkdir $bindingDir
fi

../code/python/assignBindingBioLip2RefSeq.py $alnDir \
                         ${workingDir}/../data/biolip/biolipAnn.txt \
                         $BINDING_SEQ_ID $bindingDir

######################################################################
# CLUSTER THE BINDING SITES                                          #
######################################################################

echo "4. Clustering binding sites"
echo

# get the list of binding files
find `pwd`/$bindingDir -name "*.binding" > ${workingDir}/bindingList.txt

../code/c++/clusterBindingSites ${workingDir}/bindingList.txt \
	                        $JACCARD_THRESHOLD

######################################################################
# CALCULATE THE RESIDUE WEIGHTS                                      #
######################################################################

echo "5. Calculating per-residue weights"
echo

weightDir=${bindingDir}/weights

if [ ! -d $weightDir ]; then
  mkdir $weightDir
fi

# build the list of files with the binding information
../code/c++/weightByContacts ${workingDir}/bindingList.txt \
                             ../data/biolip/receptor1 \
                             ../data/biolip/ligand \
                             $DISTANCE_WEIGHT $weightDir
                             
