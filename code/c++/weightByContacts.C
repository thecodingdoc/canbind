//////////////////////////////////////////////////////////////////////
// weightByContacts.C  -- part of the MAPBIND pipeline              //
// Author:  Dario Ghersi                                            //
// Version: 20130902                                                //
// Goal:    The program assigns a weight to each binding residue    //
//          based on the average number of contacts it has in the   //
//          different protein structures                            //
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

#include <iomanip>
#include <iterator>
#include "utils.h"
#include "weightByContacts.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

AvgContacts::AvgContacts(vector<vector<string> > allMappedRes,
			 vector<vector<double> > allWeights)
{
  // constructor for the average contacts object

  // create a vector with the names of the mapped residues
  for (unsigned int i = 0; i < allMappedRes.size(); i++) {
    for (unsigned int j = 0; j < allMappedRes[i].size(); j++) {
      bool found = false;
      for (unsigned int k = 0; k < mappedRes.size(); k++) {
	if (mappedRes[k] == allMappedRes[i][j]) {
	  found = true;
	  break;
	}
      }
      if (!found) {
	mappedRes.push_back(allMappedRes[i][j]);
	avgWeights.push_back(0.0);
      }
    }
  }

  // create a vector with the frequency of occurrence of each residue
  vector<unsigned int> freqRes(mappedRes.size(), 0);

  // calculate the average weight
  for (unsigned int i = 0; i < allWeights.size(); i++) {
    for (unsigned int j = 0; j < allWeights[i].size(); j++) {
      for (unsigned int k = 0; k < mappedRes.size(); k++) {
	if (mappedRes[k] == allMappedRes[i][j]) {
	  freqRes[k]++;
	  avgWeights[k] += allWeights[i][j];
	}
      }
    }
  }
  for (unsigned int i = 0; i < avgWeights.size(); i++) {
    avgWeights[i] /= freqRes[i];
  }
}

//////////////////////////////////////////////////////////////////////

void AvgContacts::sort()
{
  // sort by residue number

  vector<string> sortedNames(mappedRes.size());
  vector<double> sortedWeights(avgWeights.size());
  vector<intInt> nameWeight(avgWeights.size());

  int resNum;
  for (unsigned int i = 0; i < mappedRes.size(); i++) {
    resNum = atoi(mappedRes[i].substr(1,
				      mappedRes[i].length()).c_str());
    nameWeight[i] = make_pair(resNum, i);
  }
  std::sort(nameWeight.begin(), nameWeight.end(), comparator);

  for (unsigned int i = 0; i < nameWeight.size(); i++) {
    sortedNames[i] = mappedRes[nameWeight[i].second];
    sortedWeights[i] = avgWeights[nameWeight[i].second];
  }
  mappedRes = sortedNames;
  avgWeights = sortedWeights;
}

//////////////////////////////////////////////////////////////////////

bool comparator(const intInt &pair1, const intInt &pair2)
{
  // comparison function

  return pair1.first < pair2.first;
}

//////////////////////////////////////////////////////////////////////

bool checkIfClustered(string fileName)
{
  // check if there exists a clusters file for the gene

  // build a path to the putative .clusters file
  string clustersFileName = fileName;
  replace(clustersFileName, ".binding", ".clusters");

  // check if the .clusters file name exists
  bool exist = false;
  ifstream infile(clustersFileName.c_str());
  if (infile) {
    exist = true;
  }

  return exist;
}

//////////////////////////////////////////////////////////////////////

vector<double> getCoords(string line)
{
  // extract the coordinates of an atom

  vector<double> coords(3);
  coords[0] = atof(line.substr(X_POS[0], X_POS[1]).c_str());
  coords[1] = atof(line.substr(Y_POS[0], Y_POS[1]).c_str());
  coords[2] = atof(line.substr(Z_POS[0], Z_POS[1]).c_str());

  return coords;
}

//////////////////////////////////////////////////////////////////////

vector<vector<double> > getLigandCoords(string ligandFileName)
{
  // store the coordinates of the ligand into a vector of vectors

  fstream infile;
  infile.open(ligandFileName.c_str(), fstream::in);
  if (infile.fail()) {
    cerr << "Cannot open " + ligandFileName << endl;
    exit(1);
  }

  string line;
  vector<vector<double> > ligandCoords;
  while (getline(infile, line)) {
    if (line.find("ATOM") != string::npos ||
	line.find("HETATM") != string::npos) {
      vector<double> coords = getCoords(line);
      ligandCoords.push_back(coords);
    }
  }

  infile.close();

  return ligandCoords;
}

//////////////////////////////////////////////////////////////////////

void processBindingSite(string fileName, string pdbDir,
			string ligandDir, double distance,
			string outDir)
{
  // process individual binding sites

  // get the name of the pdb structure and the corresponding
  // ligand ID
  fstream infile;
  infile.open(fileName.c_str(), fstream::in);
  if (infile.fail()) {
    cerr << "Cannot open " << fileName << endl;
    exit(1);
  }
  string line;
  getline(infile, line);
  stringstream linestream(line);
  string pdb, bs, pdbRes, mappedRes, ligand;
  getline(linestream, pdb, '\t');
  getline(linestream, bs, '\t');
  getline(linestream, pdbRes, '\t');
  getline(linestream, mappedRes, '\t');
  getline(linestream, ligand, '\t');

  // close the input file
  infile.close();

  // get the ligand coordinates
  string ligandFileName = ligandDir + "/" + pdb + "_" + bs + "_" +
    ligand + ".pdb";
  vector<vector<double> > ligandCoords =
    getLigandCoords(ligandFileName);

  // build the residues object
  string pdbFileName = pdbDir + "/" + pdb + ".pdb";
  Residues residues(pdbFileName, pdbRes);

  // calculate the number of contacts per residue
  residues.calculateNumContacts(ligandCoords, distance * distance);

  // open the output file
  string temp = basename(fileName);
  replace(temp, ".binding", ".weighted");
  string outFileName = outDir + "/" + temp;
  fstream outFile;
  outFile.open(outFileName.c_str(), fstream::out);

  // print the results
  unsigned int count = 0;
  stringstream linestream2(mappedRes);
  bool beginning = true;
  string outline;

  while (getline(linestream2, temp, ' ')) {
    stringstream weight(stringstream::in | stringstream::out);
    weight << std::fixed << std::setprecision(4) <<
      residues.numContacts[count++];
    if (beginning) {
      beginning = false;
      outline = temp + ' ' + weight.str();
    }
    else {
      outline = '\t' + temp + ' ' + weight.str();
    }
    outFile << outline;
  }
  outFile << endl;
  outFile.close();

}

//////////////////////////////////////////////////////////////////////

void processClusters(string fileName, string pdbDir,
		     string ligandDir, double distance,
		     string outDir)
{
  // process clusters of binding sites

  // open the clusters file
  fstream infile;
  string temp = fileName;
  replace(temp, ".binding", ".clusters");
  infile.open(temp.c_str(), fstream::in);
  if (infile.fail()) {
    cerr << "Cannot open " << temp << endl;
    exit(1);
  }

  // open the output file
  temp = basename(fileName);
  replace(temp, ".binding", ".weighted");
  string outFileName = outDir + "/" + temp;
  fstream outFile;
  outFile.open(outFileName.c_str(), fstream::out);

  // process each cluster
  string line;
  string site;
  while (getline(infile, line)) {
    vector<vector<string> > allMappedRes;
    vector<vector<double> > allWeights;
    stringstream linestream(line);

    while (getline(linestream, site, '-')) {

      // remove trailing blank character
      if (site[0] == ' ') {
	site = site.substr(1, site.length());
      }

      // parse the site
      stringstream linestream2(site);
      string pdb, bs, pdbRes, mappedResString, ligand;
      getline(linestream2, pdb, '\t');
      getline(linestream2, bs, '\t');
      getline(linestream2, pdbRes, '\t');
      getline(linestream2, mappedResString, '\t');
      getline(linestream2, ligand, '\t');
      if (ligand[ligand.length() - 1] == ' ') {
	ligand.erase(ligand.length() - 1, 2);
      }

      // store the mapped residues
      stringstream linestream3(mappedResString);
      vector<string> mappedRes;
      while (getline(linestream3, temp, ' ')) {
	mappedRes.push_back(temp);
      }
      allMappedRes.push_back(mappedRes);

      // get the ligand coordinates
      string ligandFileName = ligandDir + "/" + pdb + "_" + bs + "_" +
	ligand + ".pdb";
      vector<vector<double> > ligandCoords =
	getLigandCoords(ligandFileName);

      // build the residues object
      string pdbFileName = pdbDir + "/" + pdb + ".pdb";
      Residues residues(pdbFileName, pdbRes);

      // calculate the number of contacts per residue
      residues.calculateNumContacts(ligandCoords,
				    distance * distance);

      // store the residues' weights
      allWeights.push_back(residues.numContacts);
    }

    // get the average weight
    AvgContacts avgContacts(allMappedRes, allWeights);
    avgContacts.sort();

    // print the results
    for (unsigned int i = 0; i < avgContacts.avgWeights.size();
	 i++) {
      stringstream weight(stringstream::in | stringstream::out);
      weight << std::fixed << std::setprecision(4) <<
	avgContacts.avgWeights[i];
      if (i == 0) {
	outFile << avgContacts.mappedRes[i] << " " << weight.str();
      }
      else {
	outFile << "\t" << avgContacts.mappedRes[i] << " " <<
	  weight.str();
      }
    }
    outFile << endl;
  }
  outFile.close();
}

//////////////////////////////////////////////////////////////////////

void Residues::calculateNumContacts(vector<vector<double> >ligandCoords,
				    double sqDist)
{
  // calculate the number of contacts per residue
  double xDiff, yDiff, zDiff;
  for (unsigned int resIndex = 0; resIndex < coords.size();
       resIndex++) {
    unsigned contacts = 0;
    for (unsigned int resAtomIndex = 0;
	 resAtomIndex < coords[resIndex].size(); resAtomIndex++) {
      for (unsigned int ligandIndex = 0;
	   ligandIndex < ligandCoords.size(); ligandIndex++) {
	// calculate the squared euclidean norm between residue atom
	// coordinates and ligand atom coordinates
	xDiff = coords[resIndex][resAtomIndex][0] -
	  ligandCoords[ligandIndex][0];
	xDiff *= xDiff;

	yDiff = coords[resIndex][resAtomIndex][1] -
	  ligandCoords[ligandIndex][1];
	yDiff *= yDiff;

	zDiff = coords[resIndex][resAtomIndex][2] -
	  ligandCoords[ligandIndex][2];
	zDiff *= zDiff;

	if ((xDiff + yDiff + zDiff) < sqDist) {
	  contacts++;
	  break; // move to the next atom
	}
      }
    }
    numContacts.push_back((double) contacts / coords[resIndex].size());
  }
}

//////////////////////////////////////////////////////////////////////

Residues::Residues(string pdbFileName, string residues)
{
  // constructor for Residues objects

  // get the numbers of the binding residues
  vector<int> numbers;
  string temp;
  stringstream linestream(residues);
  while (getline(linestream, temp, ' ')) {
    temp = temp.substr(1, temp.length() - 1);
    numbers.push_back(atoi(temp.c_str()));
  }

  // open the pdb file
  fstream pdbFile;
  pdbFile.open(pdbFileName.c_str(), fstream::in);
  if (pdbFile.fail()) {
    cerr << "Cannot open " + pdbFileName << endl;
    exit(1);
  }

  // process the pdb file
  int resnum;
  while (getline(pdbFile, temp)) {
    if (temp.find("ATOM") != string::npos ||
	temp.find("HETATM") != string::npos) {
      resnum = atoi(temp.substr(RESNUM_POS[0], RESNUM_POS[1]).c_str());

      // check if the residue is in the binding list, then store
      bool store = false;
      for (vector<int>::iterator number = numbers.begin();
	   number != numbers.end(); number++) {

	if ((*number) == resnum) {
	  store = true;
	  break;
	}
      }

      // store the coordinates
      if (store) {

	// check if the residue has been stored already
	bool isPresent = false;
	unsigned int i;
	for (i = 0; i < resNums.size(); i++) {
	  if (resNums[i] == resnum) {
	    isPresent = true;
	    break;
	  }
	}
	if (!isPresent) {
	  resNums.push_back(resnum);
	}

	// get the coordinates
	vector<double> atomCoords = getCoords(temp);
	
	if (isPresent) {
	  coords[i].push_back(atomCoords);
	}
	else {
	  vector<vector<double> > toStore;
	  toStore.push_back(atomCoords);
	  coords.push_back(toStore);
	}
      }
    }
  }
  pdbFile.close();
}
//////////////////////////////////////////////////////////////////////

static inline string &rtrim(std::string &s) {
  // trim from end
  s.erase(std::find_if(s.rbegin(), s.rend(),
    std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

//////////////////////////////////////////////////////////////////////
// MAIN                                                             //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // parse the command-line arguments
  if (argc != 6) {
    cerr << "Usage: weightByContacts BINDING_BIOLIP PDB_DIRS ...\n";
    cout << "                        LIGAND_DIR DISTANCE OUT_DIR\n";
    exit(1);
  }
  char *listOfFiles = argv[1];
  string pdbDir = argv[2];
  string ligandDir = argv[3];
  double distanceThresh = atof(argv[4]);
  string outDir = argv[5];

  // store the list of file names to process
  vector<string> fileNames = getFileToProcess(listOfFiles);

  // process each file
  for (vector<string>::iterator fileName = fileNames.begin();
       fileName != fileNames.end(); fileName++) {

    // determine whether the binding site has been clustered or
    // it is just one site
    bool isClustered = checkIfClustered(*fileName);
    if (isClustered) {
      processClusters(*fileName, pdbDir, ligandDir,
      		      distanceThresh, outDir);
    }
    else {
      processBindingSite(*fileName, pdbDir, ligandDir,
       			 distanceThresh, outDir);
    }
  }

  return 0;
}
