#!/usr/bin/python

######################################################################
# assignBindingBioLip2RefSeq.py                                      #
# Author:  Dario Ghersi                                              #
# Version: 20130828                                                  #
# Goal:    use CLUSTAL-OMEGA to build pairwise alignments of         #
#          reference sequences                                       #
#          to BioLip sequences                                       #
# Usage:   assignBindingBioLip2RefSeq.py REF_MAPPED_ONTO_BIOLIP ...  #
#                                        BIOLIP_ANNOTATION           #
#                                        CONSERVATION OUT_DIR        #
# N.B.:    the indices start from 0                                  #
######################################################################

import glob
import os
import sys

######################################################################
# CONSTANTS                                                          #
######################################################################

BIOLIP_PDBCODE = 0
BIOLIP_CHAIN = 1
BIOLIP_SITEID = 3
BIOLIP_LIGAND = 4
BIOLIP_RESIDUES = 6

EXTENSION = "aln"
POSTFIXLEN = 10

######################################################################
# FUNCTIONS                                                          #
######################################################################

def assignBinding(outfile, pdb2RefSeq, refSeqSeq, biolipInfo, pdbCode,
                  bindingInfo, consThreshold):
  """
  assign the binding information obtained from homolog structures
  """

  # assign the binding information
  sortedSites = biolipInfo[pdbCode].keys()
  sortedSites.sort()
  for site in sortedSites:
    mappedSite = []
    for res in biolipInfo[pdbCode][site]["bindingRes"]:
      resNum = int(res[1:])
      if pdb2RefSeq.has_key(resNum -1):
        refSeqRes = refSeqSeq[pdb2RefSeq[resNum - 1]] + \
            str(pdb2RefSeq[resNum - 1] + 1)
          
        mappedSite.append(refSeqRes)
        
    # print the results
    if len(mappedSite) > 0 and\
        len(mappedSite) == len(biolipInfo[pdbCode][site]["bindingRes"]):
      
      # calculate the conservation between mapped and structural site
      lenSite = len(mappedSite)
      matches = sum(map(lambda x:\
                        biolipInfo[pdbCode][site]["bindingRes"][x][0] ==\
                        mappedSite[x][0], range(lenSite)))
      conservation = 100.0 * float(matches) / lenSite

      # consider binding only if conservation of the site is >= threshold
      if conservation >= consThreshold:

        bindingInfo.extend(map(lambda x: int(x[1:]) - 1, mappedSite))
        toPrint = "\t".join([pdbCode, site,\
                             " ".join(biolipInfo[pdbCode][site]["bindingRes"]),\
                             " ".join(mappedSite),\
                             biolipInfo[pdbCode][site]["ligand"]]) + "\n"
        outfile.write(toPrint)

######################################################################

def getBindingInfo(biolipFileName):
  """
  extract the binding information from the biolip annotation file
  """

  biolipInfo = {}
  biolipFile = open(biolipFileName, "r")
  for line in biolipFile:
    fields = line[:-1].split("\t")
    
    # get the pdb id
    pdbID = fields[BIOLIP_PDBCODE] + fields[BIOLIP_CHAIN]

    # get the site id
    siteID = fields[BIOLIP_SITEID]

    # get the name of the ligand
    ligand = fields[BIOLIP_LIGAND]

    # get the binding residues (renumbered from 1)
    bindingRes = fields[BIOLIP_RESIDUES].split()

    # store the information in the dictionary
    if len(bindingRes) > 0:
      if not biolipInfo.has_key(pdbID):
        biolipInfo[pdbID] = {}

      biolipInfo[pdbID][siteID] = {}
      biolipInfo[pdbID][siteID]["ligand"] = ligand
      biolipInfo[pdbID][siteID]["bindingRes"] = bindingRes
    
  biolipFile.close()

  return biolipInfo

######################################################################

def getRefSeqSeqs(dirName, postfixLen, extension):
  """
  get the dictionary of refSeq sequences
  """

  refSeqSeqs = {}
  allFiles = glob.glob(dirName + "/*." + extension)
  for item in allFiles:
    root = os.path.basename(item)
    seqName = root[:-postfixLen]
    if refSeqSeqs.has_key(seqName):
      refSeqSeqs[seqName].append(item)
    else:
      refSeqSeqs[seqName] = [item]

  return refSeqSeqs

######################################################################

def mapPDB2RefSeq(alnFileName):
  """
  map the RefSeq sequence onto the aligned PDB sequence
  N.B. the script assumes that the RefSeq sequence is in first position
  and the PDB sequence in second position
  """

  pdb2RefSeq = {}
  alnFile = open(alnFileName, "r")

  ## get the two sequences
  refSeqSeq = ""
  pdbSeq = ""
  count = 0
  storeRefSeq = False
  storePDB = False
  for line in alnFile:
    if line[0] == ">":
      if count == 0:
        storeRefSeq = True
      elif count == 1:
        storeRefSeq = False
        storePDB = True
      else: # stop, we got the two sequences
        break
      count += 1
    elif storeRefSeq:
      refSeqSeq += line[:-1]
    elif storePDB:
      pdbSeq += line[:-1]
  alnFile.close()

  ## do the mapping
  numPos = len(refSeqSeq)
  countRefSeq = 0
  countPDB = 0
  for i in range(numPos):
    if refSeqSeq[i] != "-" and pdbSeq[i] != "-" and refSeqSeq[i] != "*":
      pdb2RefSeq[countPDB] = countRefSeq
      
    # increase the counts
    if refSeqSeq[i] != "-" and refSeqSeq[i] != "*":
      countRefSeq += 1
    if pdbSeq[i] != "-":
      countPDB += 1

  refSeqSeq = refSeqSeq.replace("-", "")
  refSeqSeq = refSeqSeq.replace("*", "")
  
  return pdb2RefSeq, refSeqSeq

######################################################################

def removeEmptyFiles(outDir):
  """
  remove empty .binding files, and their corresponding .struct files
  """

  bindingFiles = glob.glob(outDir + "/*.binding")

  for item in bindingFiles:
    ## count the number of lines in the file
    numLines = sum(1 for line in open(item, "r"))

    ## remove the file if empty, and its corresponding struct file
    if numLines == 0:
      os.remove(item)
      structFileName = item.split(".binding")[0] + ".struct"
      os.remove(structFileName)

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 5:
  print "Usage: assignBindingRefSeq.py REF_MAPPED_ONTO_BIOLIP ..."
  print "                              BIOLIP_ANNOTATION CONSERVATION OUT_DIR"
  sys.exit(1)
refSeqMappedDir, biolipFileName, conservation, outDir = sys.argv[1:]
conservation = float(conservation)

## get the dictionary of the RefSeq sequences
refSeqSeqs = getRefSeqSeqs(refSeqMappedDir, POSTFIXLEN, EXTENSION)

## get the binding information from BioLip
biolipInfo = getBindingInfo(biolipFileName)

## assign binding information to each RefSeq sequence
sortedSeqIDs = refSeqSeqs.keys()
sortedSeqIDs.sort()

for refSeqID in sortedSeqIDs:
  outfileBinding = open(outDir + "/" + refSeqID + ".binding", "w")
  bindingRes = [] # residues with binding information
  structRes = [] # residues with structural information
  for alnID in refSeqSeqs[refSeqID]:
    temp = os.path.basename(alnID).split("." + EXTENSION)[0]
    pdbCode = temp[len(temp) - 5:-1] + temp[-1]

    # map the RefSeq sequence onto the PDB sequence
    pdb2RefSeq, refSeqSeq = mapPDB2RefSeq(alnID)

    # assign the binding and write to file
    assignBinding(outfileBinding, pdb2RefSeq, refSeqSeq, biolipInfo,
                  pdbCode, bindingRes, conservation)

    # check how many residues have structural information
    structRes.extend(pdb2RefSeq.values())

  # write the sequence of the protein specifying whether the residues
  # have structural information or not
  structRes = list(set(structRes))
  outfileStruct = open(outDir + "/" + refSeqID + ".struct", "w")
  for i in range(0, len(refSeqSeq)):
    code = "U"
    if i in bindingRes:
      code = "B"
    elif i in structRes:
      code = "S"
    
    outfileStruct.write(str(i + 1) + "\t" + refSeqSeq[i] +\
                        "\t" + code + "\n")

  # close the output files
  outfileStruct.close()
  outfileBinding.close()

## remove empty files
removeEmptyFiles(outDir)
