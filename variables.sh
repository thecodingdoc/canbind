#!/bin/bash

######################################################################
# variables.sh                                                       #
# bash variables used by CANBIND                                     #
######################################################################

BLAST_DIR="/usr/local/ncbi/blast/bin/"
CLUSTALO="/Users/dghersi/Programs/clustal-omega-1.1.0/src/clustalo"

BLAST_EVALUE=1E-6 # expected value for BLAST searches
HIT_COVERAGE=80 # coverage of the structure
SEQ_IDENTITY=60 # minimum sequence identity
BINDING_SEQ_ID=90 # sequence identity in the binding site
JACCARD_THRESHOLD=0.0 # threshold to cluster binding sites
DISTANCE_WEIGHT=4.0 # cutoff for weighting the binding residues
NUM_ITER=100000 # number of randomizations for p-value calculations
