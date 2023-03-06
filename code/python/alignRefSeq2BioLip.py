#!/usr/bin/python

######################################################################
# alignRefSeq2BioLip.py                                              #
# Author:  Dario Ghersi                                              #
# Version: 20130827                                                  #
# Goal:    use CLUSTAL-OMEGA to build pairwise alignments of         #
#          reference sequences to BioLip sequences                   #
# Usage:   alignRefSeq2PDB.py BLAST_PARSED REF_SEQS BIOLIP_SEQS      #
#                             CLUSTAL_PATH IDENTITY COVERAGE OUT_DIR #
######################################################################

import os
import sys

######################################################################
# FUNCTIONS                                                          #
######################################################################

def storeSeqs(seqFileName):
  """
  return a dictionary with the sequences
  """

  seqDict = {}

  seqFile = open(seqFileName, "r")
  seq = ""
  currentID = ""
  for line in seqFile:
    if line[0] == ">":
      if len(seq) != 0:
        seqDict[currentID] = seq
      currentID = line[1:-1]
      seq = ""
    else:
      seq += line[:-1]
  seqFile.close()

  seqDict[currentID] = seq

  return seqDict

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 8:
  print "Usage: alignRefSeqC2BioLip.py BLAST_PARSED REF_SEQS BIOLIP_SEQS ..."
  print "                             CLUSTAL_PATH IDENTITY COVERAGE OUT_DIR"
  sys.exit(1)
blastParsedFileName, refSeqsFileName, biolipSeqsFileName,\
    clustalPath, identity, coverage, outDir = sys.argv[1:]
identity = float(identity)
coverage = float(coverage)

## store the reference and biolip sequences
refSeqs = storeSeqs(refSeqsFileName)
biolipSeqs = storeSeqs(biolipSeqsFileName)

## process each entry
blastParsedFile = open(blastParsedFileName, "r")
header = blastParsedFile.readline()
for line in blastParsedFile:
  fields = line[:-1].split()
  
  # filter based on identity and coverage on the PDB side
  if float(fields[2]) > identity and float(fields[8]) > coverage:
    refSeqID = fields[0]
    pdbID = fields[1]
    
    if refSeqs.has_key(refSeqID) and biolipSeqs.has_key(pdbID):
      # write the sequences to file
      print "Aligning " + refSeqID + " to " + pdbID
      alnName = outDir + "/" + refSeqID + "_" + pdbID[:4] + \
          pdbID[5] + ".aln"
      fastaFileName = outDir + "/temp.fasta"
      outfile = open(fastaFileName, "w")
      outfile.write(">" + refSeqID + "\n" +\
                    refSeqs[refSeqID] + "\n")
      outfile.write(">" + pdbID + "\n " + biolipSeqs[pdbID] + "\n")
      outfile.close()

      # run clustal-omega
      os.system(clustalPath + " -i " + fastaFileName +\
                " -o " + alnName)
      os.remove(fastaFileName)
      
blastParsedFile.close()
