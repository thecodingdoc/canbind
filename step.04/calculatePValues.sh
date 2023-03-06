#!/bin/bash

######################################################################
# calculatePValues.sh                                                #
# Author:  Dario Ghersi                                              #
# Version: 20130924                                                  #
# Goal:    calculate empirical p-values per-site                     #
# N.B.:    a p-value of -1 means that no mutation in a residue with  #
#          structural information was found. By default a value of 1 #
#          is assigned to binding sites that do not have mutations   #
# Usage:   calculatePValues CANCER_TYPE                              #
######################################################################

source ../variables.sh

cancerType=$1
if [ -z $cancerType ]; then
  echo "Please specify the cancer type"
  exit
fi

## create the structure files and weight files
path=`pwd`
absPath=`dirname $path`
bindingDir=${absPath}/step.03/bindingInfoBioLip_ID${SEQ_IDENTITY}_COV${HIT_COVERAGE}_BINDID${BINDING_SEQ_ID}
find $bindingDir -name "*.struct" > structFiles
find ${bindingDir}/weights -name "*.weighted" > weightFiles

## create the mutation file
cancerType=`basename $cancerType`
mutDir=${absPath}/step.01/mutDir/$cancerType
if [ ! -d $mutDir ]; then
  echo "Mutation directory does not exist"
  exit
fi
find $mutDir -name "*.mut" > mutFiles

## run the p-value calculations
../code/c++/calculatePValues mutFiles structFiles weightFiles \
                             $NUM_ITER ${cancerType}_pvals.txt