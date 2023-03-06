#!/usr/bin/python

###############################################################################
# buildSeqDatabase.py                                                         #
# Author:  Dario Ghersi                                                       #
# Version: 20130320                                                           #
# Goal:    build a fasta file with all the sequences that have at least a     #
#          missense mutation                                                  #
# Usage:   buildSeqDatabase.py CANCERS_DIR FASTA_DIR OUTPUT                   #
#                                                                             #
###############################################################################

import glob
import os
import sys

###############################################################################
# FUNCTIONS                                                                   #
###############################################################################

def storeFastaSeqs(fastaDir):
  """
  store the fasta sequences
  """

  seqs = {}
  files = glob.glob(fastaDir + "/*.fasta")
  for item in files:
    infile = open(item, "r")
    infile.readline() # skip the header
    geneName = os.path.basename(item).split(".fasta")[0]
    seqs[geneName] = ""
    for line in infile:
      seqs[geneName] += line[:-1]
    infile.close()

  return seqs

###############################################################################
# MAIN PROGRAM                                                                #
###############################################################################

## parse the parameters
if len(sys.argv) != 4:
  print "Usage: buildSeqDatabase.py CANCERS_DIR FASTA_DIR OUTPUT"
  sys.exit(1)
cancersDir, fastaDir, outfileName = sys.argv[1:]

## process each cancer type
cancerTypes = glob.glob(cancersDir + "/*")
print cancerTypes
toConsider = []
for cancerType in cancerTypes:
  files = glob.glob(cancerType + "/*.mut")
  for item in files:
    size = os.stat(item)[6]
    if size > 0:
      geneName = os.path.basename(item).split(".mut")[0]
      toConsider.append(geneName)

## remove redundancy
toConsider = list(set(toConsider))
toConsider.sort()

## store the fasta sequences
seqs = storeFastaSeqs(fastaDir)

## print the results
outfile = open(outfileName, "w")
for geneName in toConsider:
  seq = seqs[geneName]
  outfile.write(">" + geneName + "\n")
  if len(seq) <= 70:
    outfile.write(seq + "\n")
  else:
    for i in range(0, len(seq), 70):
      outfile.write(seq[i:i+70] + "\n")

outfile.close()
