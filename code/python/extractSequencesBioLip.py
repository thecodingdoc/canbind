#!/usr/bin/python

######################################################################
# extractSequencesBioLip.py                                          #
# Author:  Dario Ghersi                                              #
# Version: 20130827                                                  #
# Goal:    extract the sequences from BioLIP and put them in FASTA   #
#          format                                                    #
# Usage:   extractSequencesBioLip.py BIOLIP_ANNOTATION OUTPUT        #
######################################################################

import sys

######################################################################
# CONSTANTS                                                          #
######################################################################

MAX_SEQ_LINE = 70

######################################################################
# FUNCTIONS                                                          #
######################################################################

def storeBioLipSeqs(biolipAnnFileName):
  """
  store the sequences, removing duplicates
  """

  biolipSeqs = {}
  biolipAnnFile = open(biolipAnnFileName, "r")
  for line in biolipAnnFile:
    fields = line[:-1].split("\t")
    seqID = ">" + fields[0] + "_" + fields[1]
    if not biolipSeqs.has_key(seqID):
      seq = fields[-1]
      biolipSeqs[seqID] = seq
  biolipAnnFile.close()

  return biolipSeqs

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 3:
  print "Usage: extractSequencesBioLip.py BIOLIP_ANNOTATION OUTPUT"
  sys.exit(1)
biolipAnnFileName, outfileName = sys.argv[1:]

## store the sequences
biolipSeqs = storeBioLipSeqs(biolipAnnFileName)

## print the results
sortedKeys = biolipSeqs.keys()
sortedKeys.sort()
outfile = open(outfileName, "w")
for item in sortedKeys:
  seq = biolipSeqs[item]
  outfile.write(item + "\n")
  if len(seq) <= MAX_SEQ_LINE:
    outfile.write(seq + "\n")
  else:
    for i in range(0, len(seq), MAX_SEQ_LINE):
      outfile.write(seq[i:i+MAX_SEQ_LINE] + "\n")
      
outfile.close()
