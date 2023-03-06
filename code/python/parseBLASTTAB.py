#!/usr/bin/python

######################################################################
# parseBLASTTAB.py                                                   #
# Author:  Dario Ghersi                                              #
# Version: 20130124                                                  #
# Goal:    parse the tabular output of a blast run and return        #
#          sequence identity and coverage for both query and hit     #
# Usage:   parseBLASTTAB.py BLASTTAB QUERY_SEQS HIT_SEQS EVAL OUTPUT #
######################################################################

import sys

######################################################################
# CONSTANTS                                                          #
######################################################################

queryIDPos = 0
subjectIDPos = 1
identityPos = 2
queryStartPos = 6
queryEndPos = 7
hitStartPos = 8
hitEndPos = 9
evalPos = 10

######################################################################
# FUNCTIONS                                                          #
######################################################################

def getCoverage(blastTab, hitSeqs, querySeqs):
  """
  calculate the coverage on the query and hit sequence
  """

  coverageHit = {}
  coverageQuery = {}
  for query in blastTab:
    coverageHit[query] = {}
    coverageQuery[query] = {}
    for hit in blastTab[query]:
      # get the alignment length for the query
      alnLengthQuery = blastTab[query][hit]["queryEnd"] -\
      blastTab[query][hit]["queryStart"] + 1

      # get the alignment length for the hit
      alnLengthHit = blastTab[query][hit]["hitEnd"] -\
      blastTab[query][hit]["hitStart"] + 1

      # calculate the coverage
      coverageQuery[query][hit] = alnLengthQuery / \
                                  float(len(querySeqs[query]))
      coverageHit[query][hit] = alnLengthHit / \
                                float(len(hitSeqs[hit]))

  return coverageHit, coverageQuery

######################################################################

def getSequences(fastaFileName, ids, typeSeq):
  """
  store the query sequences into a dictionary
  """

  fastaFile = open(fastaFileName, "r")
  sequences = {}
  toStore = False # flag for selecting the sequences to store
  currentID = ""
  for line in fastaFile:
    if line[0] == ">":
      seqID = line[:-1].split()[0][1:]
      if typeSeq == "hit":
        seqID = seqID.split()[0]
      toStore = False
      if seqID in ids:
        toStore = True
        currentID = seqID
        sequences[seqID] = ""
    elif toStore:
      sequences[currentID] += line[:-1].replace("*", "")
  fastaFile.close()

  return sequences

######################################################################

def printResults(blastTab, coverageHit, coverageQuery, hitSeqs,\
                 querySeqs, output):
  """
  print the query sequence, hit, beginning and end, sequence length,
  identity, coverage and the e-value
  """

  ## open the output file
  outfile = open(output, "w")

  ## print the header
  outfile.write("\t".join(("QUERY", "HIT", "IDENTITY", "START", "END",
                           "QUERY_LEN",  "HIT_LEN", "QUERY_COVERAGE",
                           "HIT_COVERAGE", "EVALUE")) + "\n")

  sortedQueries = blastTab.keys()
  sortedQueries.sort()
  for query in sortedQueries:
    querySeqLen = str(len(querySeqs[query]))
    for hit in blastTab[query]:
      hitSeqLen = str(len(hitSeqs[hit]))
                                  
      # convert to string
      covHitString = "%.3f" % (100 * coverageHit[query][hit])
      covQueryString = "%.3f" % (100 * coverageQuery[query][hit])
      evalue = str(blastTab[query][hit]["eval"])
      start = str(blastTab[query][hit]["queryStart"])
      end = str(blastTab[query][hit]["queryEnd"])
      
      toPrint = [query, hit, blastTab[query][hit]["identity"], start,
                 end, querySeqLen, hitSeqLen, covQueryString,
                 covHitString, evalue]
      outfile.write("\t".join(toPrint) + "\n")

  outfile.close()

######################################################################

def processBLASTTAB(blastTabFileName, evalueCutoff):
  """
  extract the query id, the subject id, the percentage identity, the
  length of the alignment, the beginning and end of the alignment
  """

  blastTabFile = open(blastTabFileName, "r")
  blastTab = {}
  for line in blastTabFile:
    fields = line[:-1].split()
    # get the evalue
    evalue = float(fields[evalPos])
    
    if evalue <= evalueCutoff: # if evalue < threshold
      queryID = fields[queryIDPos]
      subjectID = fields[subjectIDPos]
      
      if not blastTab.has_key(queryID):
        blastTab[queryID] = {}

      # store the values
      blastTab[queryID][subjectID] =\
        {"identity": fields[identityPos],
         "queryStart": int(fields[queryStartPos]),
         "queryEnd": int(fields[queryEndPos]),
         "hitStart": int(fields[hitStartPos]),
         "hitEnd": int(fields[hitEndPos]),
         "eval": float(fields[evalPos])}
      
  blastTabFile.close()

  return blastTab

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 6:
  print "Usage: parseBLASTTAB.py BLASTTAB QUERY_SEQS HIT_SEQS EVAL OUTPUT"
  sys.exit(1)
blastTabFileName, querySeqsFileName, hitSeqsFileName, evalueCutoff,\
output = sys.argv[1:]
evalueCutoff = float(evalueCutoff)

## process the blastTab file
blastTab = processBLASTTAB(blastTabFileName, evalueCutoff)
sys.stderr.write("   blast tab output processed\n")

## store the relevant hit sequences
hits = []
for item in blastTab:
  hits.extend(blastTab[item])
hits = list(set(hits))
hitSeqs = getSequences(hitSeqsFileName, hits, "hit")
sys.stderr.write("   hit sequences stored\n")

## store the relevant query sequences
querySeqs = getSequences(querySeqsFileName, blastTab.keys(), "query")
sys.stderr.write("   query sequences stored\n")

## get the coverage
coverageHits, coverageQuery = getCoverage(blastTab, hitSeqs,\
                                          querySeqs)
sys.stderr.write("   coverage computed\n")

## print the results
printResults(blastTab, coverageHits, coverageQuery, hitSeqs,\
             querySeqs, output)
