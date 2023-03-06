#!/usr/bin/python

###############################################################################
# weightByContacts.py                                                         #
# Author:  Dario Ghersi                                                       #
# Version: 20130301                                                           #
# Goal:    assign a weight to each binding residue based on the average       #
#          number of contacts it has in the different protein structures      #
# Usage:   weightByContacts.py BINDING_BIOLIP PDB_DIRS LIGAND_DIRS ...        #
#                              DISTANCE OUT_DIR                               #
###############################################################################

import glob
import os
import sys

###############################################################################
# CONSTANTS                                                                   #
###############################################################################

PDB_POS = 0
SITE_POS = 1
PDB_RESIDUES = 2
MAPPED_RESIDUES = 3
LIG_POS = 4

X_POS = [30, 38]
Y_POS = [38, 46]
Z_POS = [46, 54]
RESNUM_POS = [22, 26]

###############################################################################
# FUNCTIONS                                                                   #
###############################################################################

def buildGene2FileDict(bindingBioLipDir):
  """
  build a gene->file dictionary. If a gene doesn't have a cluster file (i.e.,
  it has only one binding site) use the .binding file, otherwise the .clusters
  """

  gene2File = {}
  binding = set(map(lambda x: x.split(".binding")[0],\
                    glob.glob(bindingBioLipDir + "/*.binding")))
  clusters = set(map(lambda x: x.split(".clusters")[0],\
                     glob.glob(bindingBioLipDir + "/*.clusters")))

  onlyBinding = binding.difference(clusters)
  for item in onlyBinding:
    gene2File[item] = item + ".binding"
    
  for item in clusters:
    gene2File[item] = item + ".clusters"

  return gene2File

###############################################################################

def getNumberOfContacts(pdbFileName, ligandFileName, residues, distance):
  """
  get the number of contacts per residue
  """

  ## store the coordinates of the ligand
  ligandFile = open(ligandFileName, "r")
  ligandCoords = []
  for line in ligandFile:
    if line[:4] == "ATOM" or line[:6] == "HETATM":
      x = float(line[X_POS[0]:X_POS[1]])
      y = float(line[Y_POS[0]:Y_POS[1]])
      z = float(line[Z_POS[0]:Z_POS[1]])
      ligandCoords.append([x, y, z])
  ligandFile.close()

  ## get the residue numbers
  resNums = map(lambda x: x[1:], residues)

  ## store the coordinates of the residues
  pdbFile = open(pdbFileName, "r")
  selectedAtoms = []
  currentNum = ""
  currentIndex = -1
  for line in pdbFile:
    if line[:4] == "ATOM":
      resNum = line[RESNUM_POS[0]:RESNUM_POS[1]].replace(" ", "")
      if resNum in resNums:
        # get the coordinates
        x = float(line[X_POS[0]:X_POS[1]])
        y = float(line[Y_POS[0]:Y_POS[1]])
        z = float(line[Z_POS[0]:Z_POS[1]])

        # store the results
        if resNum != currentNum:
          selectedAtoms.append([])
          currentIndex += 1
        currentNum = resNum
        selectedAtoms[currentIndex].append([x, y, z])
  pdbFile.close()

  ## calculate the number of contacts per residue
  numContacts = []
  distance = float(distance)
  threshold = distance * distance
  for i in range(len(selectedAtoms)):
    contacts = 0.0
    for pdbAtom in selectedAtoms[i]:
      for ligandAtom in ligandCoords:
        xDiff = (pdbAtom[0] - ligandAtom[0])
        xDiff *= xDiff
        yDiff = (pdbAtom[1] - ligandAtom[1])
        yDiff *= yDiff
        zDiff = (pdbAtom[2] - ligandAtom[2])
        zDiff *= zDiff

        sqDist = xDiff + yDiff + zDiff

        if sqDist <= threshold:
          contacts += 1
          break

    numContacts.append(float(contacts) / len(selectedAtoms[i]))

  return numContacts

###############################################################################

def processBindingSite(infileName, pdbDir, ligandDir, distance, outDir):
  """
  process individual binding sites
  """

  ## get the binding site line
  infile = open(infileName, "r")
  data = infile.readline()[:-1]
  infile.close()

  ## get the pdb and ligand names
  fields = data.split("\t")
  pdbFileName = pdbDir + "/" + fields[PDB_POS] + ".pdb"
  ligandFileName = ligandDir + "/" + fields[PDB_POS] + "_"+ \
      fields[SITE_POS] + "_" + fields[LIG_POS] + ".pdb"
  
  ## get the number of contacts per residue
  residuesPDB = fields[PDB_RESIDUES].split()
  residuesSeq = fields[MAPPED_RESIDUES].split()
  numOfContacts = getNumberOfContacts(pdbFileName, ligandFileName,\
                                      residuesPDB, distance)
  numOfContacts = map(lambda x: "%.4f" % x, numOfContacts)

  ## print the results to file
  infileName = os.path.basename(infileName)
  outfileName = outDir + "/" + infileName.split(".binding")[0] + ".weighted"
  outfile = open(outfileName, "w")
  toPrint = map(lambda x: " ".join(x), zip(residuesSeq, numOfContacts))
  if len(toPrint) > 0 :
    outfile.write("\t".join(toPrint) + "\n")
  outfile.close()

###############################################################################

def processClusters(infileName, pdbDir, ligandDir, distance, outDir):
  """
  process clusters of binding sites
  """

  ## get the binding site info
  infile = open(infileName, "r")
  data = infile.readlines()
  infile.close()

  # open the output file
  infileName = os.path.basename(infileName)
  outfileName =  outDir + "/" + infileName.split(".clusters")[0] + ".weighted"
  outfile = open(outfileName, "w")

  ## process each binding cluster
  for line in data:
    sites = line[:-1].split(" - ")
    numOfContacts = {}
    for site in sites:
      # get the pdb and ligand file names
      fields = site.split("\t")
      pdbFileName = pdbDir + "/" + fields[PDB_POS] + ".pdb"
      ligandFileName = ligandDir + "/" + fields[PDB_POS] + "_"+ \
          fields[SITE_POS] + "_" + fields[LIG_POS] + ".pdb"

      # get the residues
      residuesPDB = fields[PDB_RESIDUES].split()
      residuesSeq = fields[MAPPED_RESIDUES].split()
      contacts = getNumberOfContacts(pdbFileName, ligandFileName,\
                                     residuesPDB, distance)
      # calculate the normalized number of contacts per residue per structure
      for i in range(len(residuesSeq)):
        if numOfContacts.has_key(residuesSeq[i]):
          numOfContacts[residuesSeq[i]].append(contacts[i])
        else:
          numOfContacts[residuesSeq[i]] = [contacts[i]]

    # calculate the average number of contacts per residue and write to file
    residues = numOfContacts.keys()
    indices = sortByResNum(residues)
    toPrint = []
    for i in indices:
      residue = residues[i]
      meanValue = float(sum(numOfContacts[residue]))\
          / len(numOfContacts[residue])
      toPrint.append(residue + " %.4f" % meanValue)

    if len(toPrint) > 0:
      outfile.write("\t".join(toPrint) + "\n")
    
    """
    # print average and number of times the residue shows up
    for res in numOfContacts:
      print float(len(numOfContacts[res])) / len(sites), float(sum(numOfContacts[res]))\
          / len(numOfContacts[res])
    """
    
  outfile.close()

###############################################################################

def sortByResNum(residues):
  """
  sort by residue number and return the indices
  """

  ## convert from X123 format to 123, discarding the amino acid name
  residues = map(lambda x: int(x[1:]), residues)
  indices = sorted(range(len(residues)), key=residues.__getitem__)

  return indices

###############################################################################
# MAIN PROGRAM                                                                #
###############################################################################

## parse the parameters
if len(sys.argv) != 6:
  print "Usage: weightByContacts.py BINDING_BIOLIP PDB_DIRS LIGAND_DIRS ..."
  print "                           DISTANCE OUT_DIR"
  sys.exit(1)
bindingBioLipDir, pdbDir, ligandDir, distance, outDir = sys.argv[1:]

## get a gene->file dictionary
gene2FileDict = buildGene2FileDict(bindingBioLipDir)

## process each file
for geneID in gene2FileDict:
  print geneID
  # determine whether we are dealing with clusters or single binding sites
  infileName = gene2FileDict[geneID]
  extension = infileName.split(".")[-1]
  if extension == "binding": # single binding site
    processBindingSite(infileName, pdbDir, ligandDir, distance, outDir)
  elif extension == "clusters": # clusters
    processClusters(infileName, pdbDir, ligandDir, distance, outDir)
  else: # unrecognized extension
    print "Unrecognized extension: " + extension
    sys.exit(1)
