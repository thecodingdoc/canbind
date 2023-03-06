#!/usr/bin/python

######################################################################
# processTCGAMut.py                                                  #
# Author:  Dario Ghersi                                              #
# Version: 20131218                                                  #
# Goal:    process a TCGA MAF file and output the FASTA sequence of  #
#          the gene (protein format) and the coordinate of the       #
#          mutations in the protein                                  # 
# Usage:   processTCGAMut.py TCGA_MAF NCBI_HSREF NCBI_GBS NP_SEQS ...#
#                            MUT_DIR OUT_FILE                        #
#                                                                    #
# Note:    the script can be made more efficient in the mapping of   #
#          mutations part                                            #
######################################################################

import glob
import os
import re
import sys

######################################################################
# CONSTANTS                                                          #
######################################################################

COMPLEMENT = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}

GENCODE = {
    'ATA':'I',    # Isoleucine
    'ATC':'I',    # Isoleucine
    'ATT':'I',    # Isoleucine
    'ATG':'M',    # Methionine
    'ACA':'T',    # Threonine
    'ACC':'T',    # Threonine
    'ACG':'T',    # Threonine
    'ACT':'T',    # Threonine
    'AAC':'N',    # Asparagine
    'AAT':'N',    # Asparagine
    'AAA':'K',    # Lysine
    'AAG':'K',    # Lysine
    'AGC':'S',    # Serine
    'AGT':'S',    # Serine
    'AGA':'R',    # Arginine
    'AGG':'R',    # Arginine
    'CTA':'L',    # Leucine
    'CTC':'L',    # Leucine
    'CTG':'L',    # Leucine
    'CTT':'L',    # Leucine
    'CCA':'P',    # Proline
    'CCC':'P',    # Proline
    'CCG':'P',    # Proline
    'CCT':'P',    # Proline
    'CAC':'H',    # Histidine
    'CAT':'H',    # Histidine
    'CAA':'Q',    # Glutamine
    'CAG':'Q',    # Glutamine
    'CGA':'R',    # Arginine
    'CGC':'R',    # Arginine
    'CGG':'R',    # Arginine
    'CGT':'R',    # Arginine
    'GTA':'V',    # Valine
    'GTC':'V',    # Valine
    'GTG':'V',    # Valine
    'GTT':'V',    # Valine
    'GCA':'A',    # Alanine
    'GCC':'A',    # Alanine
    'GCG':'A',    # Alanine
    'GCT':'A',    # Alanine
    'GAC':'D',    # Aspartic Acid
    'GAT':'D',    # Aspartic Acid
    'GAA':'E',    # Glutamic Acid
    'GAG':'E',    # Glutamic Acid
    'GGA':'G',    # Glycine
    'GGC':'G',    # Glycine
    'GGG':'G',    # Glycine
    'GGT':'G',    # Glycine
    'TCA':'S',    # Serine
    'TCC':'S',    # Serine
    'TCG':'S',    # Serine
    'TCT':'S',    # Serine
    'TTC':'F',    # Phenylalanine
    'TTT':'F',    # Phenylalanine
    'TTA':'L',    # Leucine
    'TTG':'L',    # Leucine
    'TAC':'Y',    # Tyrosine
    'TAT':'Y',    # Tyrosine
    'TAA':'_',    # Stop
    'TAG':'_',    # Stop
    'TGC':'C',    # Cysteine
    'TGT':'C',    # Cysteine
    'TGA':'_',    # Stop
    'TGG':'W',    # Tryptophan
}

MAX_SEQ_LINE = 70

NONDECIMAL = re.compile(r'[^\d]+')

VALID_CHROMOSOMES = map(lambda x: str(x), range(1, 23))
VALID_CHROMOSOMES.extend(['X', 'Y'])

######################################################################
# FUNCTIONS                                                          #
######################################################################

def getComplement(seq):
  """
  get the reverse complement of a DNA string
  """

  reverseComplement = ""

  ## reverse the string
  seq = seq[::-1]

  ## get the complement
  for base in seq:
    reverseComplement += COMPLEMENT[base]

  return reverseComplement

######################################################################

def getMutPos(refPos, mutBase, seq, cds, isoforms, mutDir, geneID,
              sampleID):
  """
  get the mutated position in the protein sequence, the reference
  amino acid and the mutated amino acid
  """

  mutPos = []
  for isoform in isoforms:
    # get the corresponding protein sequence
    seqFiles = glob.glob(mutDir + "/" + geneID + \
                         ".%03d" % (isoform + 1) + ".fasta")

    # move on if no protein sequence exists
    if len(seqFiles) == 0:
      mutPos.append(-1)
      continue
      
    currentCDS = cds[isoform][0]
    
    # check if complementation is needed
    complement = False
    currentMutBase = mutBase
    if currentCDS.find("complement") != -1:
      complement = True
      currentCDS = currentCDS.replace("complement(", "")
      try:
        currentMutBase = COMPLEMENT[mutBase]
      except:
        return False

    # get the exons
    currentCDS = currentCDS.replace("join(", "").replace(")", "").split(",")
    if complement:
      currentCDS.reverse()

    # examine each exon
    rnaSeq = ""
    spot = -1
    for exon in currentCDS:
        
      # primary strand
      if exon.find("..") != -1:
        beg, end = exon.split("..")
        beg = int(NONDECIMAL.sub('', beg))
        end = int(NONDECIMAL.sub('', end))
      else:
        beg = end = int(NONDECIMAL.sub('', exon))

      # complement
      if complement:
        if beg <= refPos and end >= refPos:
          spot = len(rnaSeq) + end - refPos
        rnaSeq += getComplement(seq[beg - 1:end])
      else:
        if beg <= refPos and end >= refPos:
          spot = len(rnaSeq) + refPos - beg
        rnaSeq += seq[beg - 1:end]

    # get the position and the mutated amino acid
    seqCount = 1
    storedPos = False
    for i in range(0, len(rnaSeq), 3):
      codon = rnaSeq[i:i+3]
      if len(codon) == 3:
        if i <= spot and (i + 3) > spot:
          mutCodon = ""
          refAA = GENCODE[codon]

          for j in range(i, i+3):
            if j == spot:
              mutCodon += currentMutBase
            else:
              mutCodon += rnaSeq[j]

          try:
            mutAA = GENCODE[mutCodon]
          except:
            return False  
          mutPos.append([seqCount, refAA, mutAA, sampleID])
          storedPos = True
        seqCount += 1

    if not storedPos:
      mutPos.append(False)

  return mutPos
    
######################################################################

def getSeqExons(seq, cds):
  """
  stitch together the exons
  """

  # check if complementation is needed
  complement = False
  if cds.find("complement") != -1:
    complement = True
    cds = cds.replace("complement(", "")

  # get the exons
  cds = cds.replace("join(", "").replace(")", "").split(",")
  if complement:
    cds.reverse()
    
  rnaSeq = ""
  for exon in cds:
    if exon.find("..") != -1:
      beg, end = exon.split("..")
      beg = int(NONDECIMAL.sub('', beg))
      end = int(NONDECIMAL.sub('', end))
    else:
      beg = end = int(NONDECIMAL.sub('', exon))
      
    if complement:
      rnaSeq += getComplement(seq[beg - 1:end])
    else:
      rnaSeq += seq[beg - 1:end]
  
  protSeq = ""
  for i in range(0, len(rnaSeq), 3):
    codon = rnaSeq[i:i+3]
    if len(codon) == 3:
      protSeq += GENCODE[rnaSeq[i:i+3]]

  return protSeq

######################################################################

def mapMuts(tcgaMut, npSeqs, hsrefDir, gbsDir, seqDir):
  """
  get the protein sequences for the genes of interest and then map the
  mutations onto them, printing the results
  """

  mutations = {}
  validIsoforms = {}
  ## process each chromosome in turn
  for chromosome in tcgaMut:

    # store chromosome information
    seq, cds, synonyms = storeChromosome(hsrefDir, gbsDir, chromosome)

    # process each mutation in turn
    for mut in tcgaMut[chromosome]:
      geneID = mut[0]
      refBase = mut[3]
      refPos = mut[2]
      mutBase = mut[4]
      if mut[4] == mut[3]:
        mutBase = mut[5]

      # ignore cases where the reference base doesn't match the
      # chr location
      try:
        seq[refPos - 1]
      except:
        continue
      if seq[refPos - 1] != refBase:
        continue

      # check if the gene information is available in the CDS
      if cds.has_key(geneID):
        geneName = geneID
      elif synonyms.has_key(geneID) and cds.has_key(synonyms[geneID]):
        geneName = synonyms[geneID]
      else:
        continue # nothing to be done, move on to the next

      # write the sequences
      isoforms = writeSeq(seq, cds[geneName], npSeqs, geneName, seqDir)
      validIsoforms[geneName] = isoforms

      # get the mutation info
      if len(isoforms) > 0:
        mutPos = getMutPos(refPos, mutBase, seq, cds[geneName],\
                           isoforms, seqDir, geneName, mut[1])
        if mutPos == False:
          continue

        # store the results
        if mutations.has_key(geneName):
          mutations[geneName].append(mutPos)
        else:
          mutations[geneName] = [mutPos]

  return mutations, validIsoforms

######################################################################

def processHeader(header):
  """
  return the positions for columns of interest
  """

  pos = {}
  fields = header.split("\t")
  for i in range(len(fields)):
    if fields[i] == "Hugo_Symbol":
      pos["genePos"] = i
    elif fields[i] == "Chrom":
      pos["chromPos"] = i
    elif fields[i] == "Start_Position":
      pos["mutPos"] = i
    elif fields[i] == "Variant_Classification":
      pos["classPos"] = i
    elif fields[i] == "Tumor_Sample_Barcode":
      pos["sampleID"] = i
    elif fields[i] == "Reference_Allele":
      pos["refAll"] = i
    elif fields[i] == "Tumor_Seq_Allele1":
      pos["tumAll1"] = i
    elif fields[i] == "Tumor_Seq_Allele2":
      pos["tumAll2"] = i

  return pos

######################################################################

def processTCGAFile(mafFileName):
  """
  process the TCGA MAF file and store the following unique tuples for
  missense mutations organized by chromosome:
  CHROMOSOME = [[GENE1, SAMPLE_ID, POS, REF_ALLELE TUMOR_ALLELE],
                [GENE2, SAMPLE_ID, POS, REF_ALLELE TUMOR_ALLELE]]
  """

  mutData = {}
  mafFile = open(mafFileName)
  header = mafFile.readline() # get the header

  ## process the header
  pos = processHeader(header)

  ## process the entries
  for line in mafFile:
    fields = line[:-1].split("\t")
    if fields[pos["classPos"]] == "Missense_Mutation": # check if missense mut.
      chromosome = fields[pos["chromPos"]]
      if chromosome in VALID_CHROMOSOMES: # check if valid chromosome
          
        # collect info on the mutation
        gene = fields[pos["genePos"]]
        mutPos = int(fields[pos["mutPos"]])
        sampleID = fields[pos["sampleID"]]
        sampleID = "-".join(sampleID.split("-")[:4])
        refAll = fields[pos["refAll"]]
        tumAll1 = fields[pos["tumAll1"]]
        tumAll2 = fields[pos["tumAll2"]]

        # add the mutation if not already present
        mutLine = [gene, sampleID, mutPos, refAll, tumAll1, tumAll2]
        if not mutData.has_key(chromosome):
          mutData[chromosome] = [mutLine]
        elif not mutLine in mutData[chromosome]:
          mutData[chromosome].append(mutLine)
  
  mafFile.close()

  return mutData

######################################################################

def storeChromosome(hsrefDir, gbsDir, chromosome):
  """
  store the sequence of the chromosome and the CDS data for all genes
  """

  ## set up the variables
  seq = ""
  cds = {}
  synonyms = {}

  # get the sequence data
  seqFileName = glob.glob(hsrefDir + "/*chr" + chromosome + ".fa")[0]
  seqFile = open(seqFileName)
  header = seqFile.readline()
  for line in seqFile:
    seq += line[:-1]

  seqFile.close()

  ## get the CDS information
  gbsFileName = glob.glob(gbsDir + "/*chr" + chromosome + ".gbs")[0]
  gbsFile = open(gbsFileName, "r")
  storeInstr = False
  instructions = ""
  fieldHeader = ""
  
  for line in gbsFile:
    if len(line[:21].replace(" ", "")) > 0:
      fieldHeader = line[:21]
      storeInstr = False

    if fieldHeader.find(" CDS ") != -1: # CDS SECTION
      if line.find(" CDS ") != -1:
        storeInstr = True
        instructions = line[:-1].split()[-1]

      elif line.find("/gene=") != -1: # gene name
        storeInstr = False
        geneName = line[:-1].split("/gene=")[1]
        geneName = geneName.replace('"', '')

      elif line.find("/protein_id=") != -1: # protein ID
        proteinID = line[:-1].split("/protein_id=")[1]
        proteinID = proteinID.replace('"', '')

        # store the isoform
        if cds.has_key(geneName):
          cds[geneName].append((instructions, proteinID))
        else:
          cds[geneName] = [(instructions, proteinID)]

      elif storeInstr:
        instructions += line[:-1].split()[-1]

    elif fieldHeader.find(" gene ") != -1: # gene SECTION
      if line.find("/gene=") != -1:
        geneName = line[:-1].split("/gene=")[1]
        geneName = geneName.replace('"', '')
      elif line.find("/gene_synonym=") != -1:
        syns = line[:-1].split("/gene_synonym=")[1].replace('"', '')
        syns = syns.replace(" ", "").split(";")
        for syn in syns:
          synonyms[syn] = geneName

  gbsFile.close()

  ## remove redundancy from cds
  for gene in cds:
    cds[gene] = list(set(cds[gene]))

  return seq, cds, synonyms

######################################################################

def storeNPSeqs(npSeqsFileName):
  """
  store the protein sequences
  """

  npSeqs = {}
  protID = ""
  seq = ""
  npSeqsFile = open(npSeqsFileName, "r")
  for line in npSeqsFile:
    if line[0] == ">":
      if len(seq) > 0: # store sequence
        npSeqs[protID] = seq
      seq = ""
      protID = line.split("|")[3]
    else:
      seq += line[:-1]
  npSeqsFile.close()

  npSeqs[protID] = seq

  return npSeqs

######################################################################

def writeResults(mutations, isoforms, mutDir, outfileName):
  """
  write individual files for each isoforms, with the following format
  POSITION FREQUENCY REF_AA MUTATED_AA
  """

  ## open the output file
  globalOutFile = open(outfileName, "w")

  ## process each gene in turn
  sortedGenes = mutations.keys()
  sortedGenes.sort()
  for gene in sortedGenes:
    # aggregate the mutations per isoform
    numIsoforms = len(mutations[gene][0])
    isoMuts = []
    for isoform in range(numIsoforms):
      isoMuts.append({})
    for mutation in mutations[gene]:
      for i in range(numIsoforms):
        if mutation[i]:
          if mutation[i][1] != mutation[i][2]: # make sure it is really missense
            # write it in the output file
            currIso = isoforms[gene][i] + 1
            fullName = gene + "." + "%03d" % currIso
            toWrite = "\t".join([fullName, mutation[i][3],
                                 str(mutation[i][0]),
                                 mutation[i][1], mutation[i][2]])
            globalOutFile.write(toWrite + "\n")
            
            if isoMuts[i].has_key(mutation[i][0]):
              isoMuts[i][mutation[i][0]] += 1
            else:
              isoMuts[i][mutation[i][0]] = 1

    # print the results to file
    for i in range(numIsoforms):
      isoform = isoforms[gene][i] + 1
      sortedKeys = isoMuts[i].keys()
      if len(sortedKeys) > 0:
        outfileName = mutDir + "/" + gene + ".%03d" % (isoform) + ".mut"
        outfile = open(outfileName, "w")
        sortedKeys.sort()
        for item in sortedKeys:
          outfile.write("%d\t%d\n" % (item, isoMuts[i][item]))
        outfile.close()

  globalOutFile.close()

######################################################################

def writeSeq(seq, cdsInfo, npSeqs, geneID, seqDir):
  """
  write the protein sequence of the gene to file and return the
  isoforms that have been written successfully
  """

  isoforms = []
  ## make sure the file doesn't already exist
  for i in range(len(cdsInfo)): # process each isoform
    protSeq = getSeqExons(seq, cdsInfo[i][0])[:-1]

    # if there is no available protein sequence, move on
    if npSeqs.has_key(cdsInfo[i][1]):
      if protSeq == npSeqs[cdsInfo[i][1]]: # make sure the transcripts match

        # build the sequence file name
        outfileName = seqDir + "/" + geneID + "." +\
            "%03d" % (i+1) + ".fasta"

        # don't overwrite files
        if not os.path.exists(outfileName):
          outfile = open(outfileName, "w")
          outfile.write(">" + cdsInfo[i][1] + "\n")
          
          # write the sequence
          if len(protSeq) <= MAX_SEQ_LINE:
            outfile.write(protSeq + "\n")
          else:
            for j in range(0, len(protSeq), MAX_SEQ_LINE):
              outfile.write(protSeq[j:j + MAX_SEQ_LINE] + "\n")

          outfile.close()

        isoforms.append(i)

  return isoforms

######################################################################
# MAIN PROGRAM                                                       #
######################################################################

## parse the parameters
if len(sys.argv) != 7:
  print "Usage:   processTCGAMut.py TCGA_MAF NCBI_HSREF NCBI_GBS NP_SEQS ... "
  print "                           MUT_DIR OUT_FILE"
  sys.exit(1)
mafFileName, hsrefDir, gbsDir, npSeqsFileName, mutDir,\
outFileName = sys.argv[1:]
seqDir = "./refSeqs"

## process the maf file
tcgaMut = processTCGAFile(mafFileName)

## store NP sequences
npSeqs = storeNPSeqs(npSeqsFileName)

## map the mutations and build protein sequences
mutations, isoforms = mapMuts(tcgaMut, npSeqs, hsrefDir, gbsDir, seqDir)

## write the results
writeResults(mutations, isoforms, mutDir, outFileName)
