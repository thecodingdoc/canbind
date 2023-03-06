//////////////////////////////////////////////////////////////////////
// calculatePValues.C  -- part of the MAPBIND pipeline              //
// Author:  Dario Ghersi                                            //
// Version: 20130830                                                //
// Goal:    The program calculates empirical p-values               //
//                                                                  //
// MAPBIND is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as          //
// published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.              //
//                                                                  //
// MAPBIND is distributed in the hope that it will be useful,       //
// but WITHOUT ANY WARRANTY; without even the implied warranty of   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    //
// GNU General Public License for more details.                     //
//                                                                  //
// You should have received a copy of the GNU General Public        //
// License along with BLOCKS.                                       //
// If not, see <http://www.gnu.org/licenses/>.                      //
//////////////////////////////////////////////////////////////////////

#include "utils.h"
using namespace std;

//////////////////////////////////////////////////////////////////////
// STRUCTURES                                                       //
//////////////////////////////////////////////////////////////////////

struct Mutations {
  string isofName;
  vector<string> positions;
  vector<unsigned int> frequencies;
};

struct Structs {
  string isofName;
  vector<string> positions;
  vector<string> binding;
};

struct Site {
  vector<string> positions;
  vector<double> scores;
};

struct Weights {
  string isofName;
  vector<Site> sites;
};

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

vector<double> calculateBindingScores(Mutations, Weights);
void calculatePVals(vector<Mutations>, vector<Structs>,
		    vector<Weights>, char *);
unsigned int calculateStructMut(Mutations, Structs);
string getIsoformName(string );
vector<double> getRandomScores(Structs, Weights,
			       unsigned int, unsigned int,
			       vector<double>);
Site getScores(string);
vector<string> storeFileNames(char *);
vector<Mutations> storeMutations(vector<string>);
vector<Structs> storeStruct(vector<string>);
vector<Weights> storeWeights(vector<string>);

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////

void calculatePVals(vector<Mutations> muts, vector<Structs> structs,
		    vector<Weights> weights, unsigned int numIter,
		    char *outfileName)
{

  vector<double> bindingScores;

  // open the output file
  fstream outfile;
  outfile.open(outfileName, fstream::out);

  // initialize the random number generator
  srand(time(NULL));

  // process each isoform in turn
  for (unsigned int i = 0; i < muts.size(); i++) {

    // check if there is structural information
    Weights w;
    for (unsigned int j = 0; j < weights.size(); j++) {
      if (muts[i].isofName == weights[j].isofName) {

	// calculate the binding scores for each site
	w = weights[j];
	bindingScores = calculateBindingScores(muts[i], w);
	break;
      }
    }

    // check if there are mutations in a binding site
    bool mutBind = false;
    for (std::vector<double>::iterator it = bindingScores.begin();
	 it != bindingScores.end(); it++) {
      if (*it > 0.0) {
	mutBind = true;
	break;
      }
    }

    // get the number of binding sites
    unsigned int numSites = w.sites.size();
    unsigned int site;

    // calculate how many mutations fall in a structural part
    unsigned int numMutStruct = 0;
    Structs s;
    vector<double> pvals (numSites, 0.0);
    for (unsigned int j = 0; j < structs.size(); j++) {
      if (muts[i].isofName == structs[j].isofName) {
	s = structs[j];
	numMutStruct = calculateStructMut(muts[i], s);
	break;
      }
    }

    // get the p-values
    if (mutBind) { // if there are mutations in binding sites
      // calculate the p-values

      pvals = getRandomScores(s, w, numMutStruct, numIter,
			      bindingScores);

    }
    // no mutations in binding sites but in structural part
    if (!mutBind && numMutStruct > 0) {
      fill(pvals.begin(), pvals.end(), 1.0); 
    }
    else if (numMutStruct == 0) {
      fill(pvals.begin(), pvals.end(), -1.0);
    }

    // write the results
    for (site = 0; site < numSites; site++) {
      outfile << muts[i].isofName << " " << pvals[site] << endl;
    }
  }

  outfile.close();
}

//////////////////////////////////////////////////////////////////////

vector<double> calculateBindingScores(Mutations m, Weights w) {
  // calculate the binding scores for every binding site

  vector<double> scores (w.sites.size(), 0.0);

  // calculate the score for each site
  for (unsigned int i = 0; i < w.sites.size(); i++) {
    for (unsigned int j = 0; j < w.sites[i].positions.size(); j++) {
      for (unsigned int k = 0; k < m.positions.size(); k++) {
	if (m.positions[k] == w.sites[i].positions[j]) {
	  scores[i] += (double) m.frequencies[k] * w.sites[i].scores[j];
	}
      }
    }
  }

  return scores;
}

//////////////////////////////////////////////////////////////////////

unsigned int calculateStructMut(Mutations m, Structs s)
{
  unsigned int numStructMut = 0;

  for (unsigned int i = 0; i < s.positions.size(); i++) {
    for (unsigned int j = 0; j < m.positions.size(); j++) {
      if (s.positions[i] == m.positions[j]) {
	numStructMut += m.frequencies[j];
      }
    }
  }

  return numStructMut;
}

//////////////////////////////////////////////////////////////////////

string getIsoformName(string s)
{
  // extract the isoform name from a full path

  string temp = basename(s);
  int lastindex = temp.find_last_of(".");
  string isoformName = temp.substr(0, lastindex);

  return isoformName;
}

//////////////////////////////////////////////////////////////////////

vector<double> getRandomScores(Structs s, Weights w,
			       unsigned int numMutStruct,
			       unsigned int numIter,
			       vector<double> bindingScores)
{
  unsigned int numSites = w.sites.size();

  // calculate the p-values
  vector<double> pvals (numSites, 0.0);

  // proceed only if there are mutations in a structural part
  unsigned int sizeStruct = s.positions.size();
  unsigned int randomInt, site;

  // get the random scores
  unsigned int k, q;
  for (unsigned int j = 0; j < numIter; j++) {
    vector<double> randomScores (numSites, 0.0);
    for (k = 0; k < numMutStruct; k++) {
      randomInt = rand() % sizeStruct;

      for (site = 0; site < numSites; site++) {
	for (q = 0; q < w.sites[site].positions.size(); q++) {
	  if (w.sites[site].positions[q] == 
	      s.positions[randomInt]) {
	    randomScores[site] += w.sites[site].scores[q];
	  }
	}
      }
    }

    // calculate if random scores >= real ones
    for (k = 0; k < numSites; k++) {
      if (randomScores[k] >= bindingScores[k]) {
	pvals[k]++;
      }
    }
  }
  for (site = 0; site < numSites; site++) {
    pvals[site] = (pvals[site] + 1) / (numIter + 1);
  }

  return pvals;
}

//////////////////////////////////////////////////////////////////////

Site getScores(string s)
{

  Site site;
  size_t pos = 0;
  string token;
  string position;
  double score;
  const string delimiter = "\t";

  // process each residue
  while ((pos = s.find(delimiter)) != string::npos) {
    token = s.substr(0, pos);
    s.erase(0, pos + delimiter.length());
    
    // get the position and score
    istringstream iss(token);
    iss >> position >> score;
    site.positions.push_back(position.substr(1, string::npos));
    site.scores.push_back(score);
  }

  // process the last residue
  istringstream iss(s);
  iss >> position >> score;
  site.positions.push_back(position.substr(1, string::npos));
  site.scores.push_back(score);

  return site;
}

//////////////////////////////////////////////////////////////////////

vector<string> storeFileNames(char *fileName)
{

  // open the input file
  fstream infile;
  infile.open(fileName, fstream::in);

  // count how many lines are in the file
  string line;
  unsigned int numLines = 0;
  while (getline(infile, line))
    numLines++;

  // store the lines
  vector<string> allFileNames (numLines);
  infile.clear();
  infile.seekg (0, infile.beg);
  int pos = 0;
  while (getline(infile, line)) {
    allFileNames[pos] = line;
    pos++;
  }
  infile.close();

  return allFileNames;
}

//////////////////////////////////////////////////////////////////////

vector<Mutations> storeMutations(vector<string> mutFileNames)
{
  // store the mutations for each isoform
  vector<Mutations> muts (mutFileNames.size());

  for (unsigned int i = 0; i < mutFileNames.size(); i++) {

    // get the isoform name
    muts[i].isofName = getIsoformName(mutFileNames[i]);
    
    // store the file content
    fstream infile;
    string line;
    string pos;
    unsigned int freq;
    infile.open(mutFileNames[i].c_str(), fstream::in);
    if (infile.fail()) {
      cerr << "Cannot open " << mutFileNames[i] << endl;
      exit(1);
    }
    while (getline(infile, line)) {
      istringstream iss(line);
      if (! (iss >> pos >> freq)) {
	cerr << "Problem with mutation file " <<
	  mutFileNames[i] << endl;
	exit(1);
      }
      muts[i].positions.push_back(pos);
      muts[i].frequencies.push_back(freq);
    }

    infile.close();
  }

  return muts;
}

//////////////////////////////////////////////////////////////////////

vector<Structs> storeStruct(vector<string> structFileNames) {
  // store the structural data

  vector<Structs> structs (structFileNames.size());

  for (unsigned int i = 0; i < structFileNames.size(); i++) {

    // get the isoform name
    structs[i].isofName = getIsoformName(structFileNames[i]);

    // store the file content
    fstream infile;
    string line;
    string pos;
    char aa, code;
    infile.open(structFileNames[i].c_str(), fstream::in);
    if (infile.fail()) {
      cerr << "Cannot open " << structFileNames[i] << endl;
      exit(1);
    }
    while (getline(infile, line)) {
      istringstream iss(line);
      if (! (iss >> pos >> aa >> code)) {
	cerr << "Problem with structural file " <<
	  structFileNames[i] << endl;
	exit(1);
      }
      if (code == 'S' or code == 'B') {
        structs[i].positions.push_back(pos);
      }
      if (code == 'B') {
	structs[i].binding.push_back(pos);
      }
    }
  }

  return structs;
}

//////////////////////////////////////////////////////////////////////

vector<Weights> storeWeights(vector<string> weightFileNames)
{

  vector<Weights> weights (weightFileNames.size());

  for (unsigned int i = 0; i < weightFileNames.size(); i++) {

    // get the isoform name
    string isofName = getIsoformName(weightFileNames[i]);
    weights[i].isofName = isofName;

    // store the file content
    fstream infile;
    string line;
    infile.open(weightFileNames[i].c_str(), fstream::in);
    if (infile.fail()) {
      cerr << "Cannot open " << weightFileNames[i] << endl;
      exit(1);
    }
    while (getline(infile, line)) {
      weights[i].sites.push_back(getScores(line));
    }
  }

  return weights;
}

//////////////////////////////////////////////////////////////////////
// MAIN                                                             //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // parse the command-line arguments
  if (argc != 6) {
    cerr << "Usage: calculatePValues MUT_FILES STRUCT_FILES ";
    cerr << "WEIGHT_FILES NUM_ITER OUTFILE_NAME\n";
    exit(1);
  }

  // store the mutation file names
  vector<string> mutFileNames = storeFileNames(argv[1]);

  // store the mutations
  vector<Mutations> muts = storeMutations(mutFileNames);

  // store the struct file names
  vector<string> structFileNames = storeFileNames(argv[2]);

  // store the struct data
  vector<Structs> structs = storeStruct(structFileNames);

  // store the weight file names
  vector<string> weightFileNames = storeFileNames(argv[3]);

  // store the weights
  vector<Weights> weights = storeWeights(weightFileNames);

  // calculate the p-values
  unsigned int numIter = atoi(argv[4]);
  calculatePVals(muts, structs, weights, numIter, argv[5]);

  return 0;
}
