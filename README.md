# CANBIND v0.1

This repository contains all the code and datasets necessary to reproduce the results in:

Dario Ghersi, Mona Singh
*"Interaction-based discovery of functionally important genes in cancers"*,
Nucleic Acids Research, 2014, (42):3

----------------------------------------------------------------------

DISCLAIMER:
The programs are distributed in the hope that they will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.

----------------------------------------------------------------------

+Requirements:

1. NCBI BLAST+ (http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

2. Clustal Omega (http://www.clustal.org/omega/)

3. A C++ compiler (GNU gcc recommended) with make

In Ubuntu Linux, both BLAST+ and Clustal Omega are available as "apt"
packages ready to be installed.

----------------------------------------------------------------------

+Setting up the pipeline:

0. Get the data folder

From the main repository folder, type the following commands:

```
wget https://www.dropbox.com/s/5v1xqqerum0f5qj/data.tar.bz2?dl=0
tar xjvf data.tar.bz2
```

1. Setting up the variables

Edit the "variables.sh" file and edit the following two variables:

BLASTP: path to the bin directory of the BLAST+ suite
CLUSTALO: path to the CLUSTAL-OMEGA executable

The other variables can also be changed, but the pipeline has been
tested only with the default values.

2. Compiling the C++ programs

To compile the required C++ programs, type:

```
$ cd code/c++
$ make
```

You should end up with three executables in the same directory
("calculatePValues", "clusterBindingSites", and "weightByContacts")

----------------------------------------------------------------------

The pipeline consists of four sequential steps (contained in separate
directories).

1. The first step parses The Cancer Genome Atlas
raw mutation files, and create reference sequences and mutation files.
The files produced in this step can be replaced with other files
generated directly by the user. For example, the user may have
mutation data that comes from other sources. This data can be used,
but has to be in the same format as the parsed TCGA data.
See the "step.01" directory for more details and an example.

2. The second step finds structures that match the reference
sequences obtained in step 1. Please note that the BLAST search may
take a long time, depending on the number of reference sequences.

3. The third step assigns binding information from BioLip to
the reference sequences and weights the residues by proximity to
the ligands.

4. In the last (optional) step empirical p-values are computed.

Please see the README files in each directory for further information
on how to carry out these steps, and do not hesitate to contact us
if you have further questions at dghersi@unomaha.edu.
