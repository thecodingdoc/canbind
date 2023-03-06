//////////////////////////////////////////////////////////////////////
// weightByContacts.h  -- part of the MAPBIND pipeline              //
// Author:  Dario Ghersi                                            //
// Version: 20130830                                                //
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

//////////////////////////////////////////////////////////////////////
// CONSTANTS AND VARIABLES                                          //
//////////////////////////////////////////////////////////////////////

const unsigned int PDB_POS = 0, SITE_POS = 1, PDB_RESIDUES = 2;
const unsigned int MAPPED_RESIDUES = 3, LIG_POS = 4;

const unsigned int X_POS[] = {30, 8};
const unsigned int Y_POS[] = {38, 8};
const unsigned int Z_POS[] = {46, 8};
const unsigned int RESNUM_POS[] = {22, 4};

typedef std::pair<int,int> intInt;

//////////////////////////////////////////////////////////////////////
// CLASSES                                                         //
//////////////////////////////////////////////////////////////////////

class AvgContacts {

 public:
  vector<string> mappedRes;
  vector<double> avgWeights;

  AvgContacts(vector<vector<string> >, vector<vector<double > >);
  void sort();
};

//////////////////////////////////////////////////////////////////////

class Residues {

 public:
  vector<int> resNums;
  vector<vector<vector<double> > > coords;
  vector<double> numContacts;

  Residues(string, string);
  void calculateNumContacts(vector<vector<double> >, double);
};

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

bool comparator(const intInt &, const intInt &);
bool checkIfClustered(string);
vector<double> getCoords(string);
void processBindingSite(string, string, string, double, string);
void processClusters(string, string, string, double, string);
static inline string &rtrim(string &);

