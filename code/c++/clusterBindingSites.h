//////////////////////////////////////////////////////////////////////
// clusterBindingSites.h                                            //
// Author:  Dario Ghersi                                            //
// Version: 20130829                                                //
// Goal:    part of the MapBind suite. The program calculate        //
//          pairwise Jaccard distance between binding sites.        //
//          Then, it builds a graph based on a Jaccard threshold    //
//          and identifies the connected components                 //
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

const unsigned int PDB_ID = 0;
const unsigned int BINDING_SITE_ID = 1;
const unsigned int BINDING_SITE_RES = 3;

typedef boost::adjacency_list<boost::vecS, boost::vecS,\
  boost::undirectedS> MyGraph;
vector<unsigned int> visitedVertices; // global variable, to
                                     // store BFS results

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

void clusterSites(string, double);
set<string> getSite(string);
vector<vector<int> > findClusters(MyGraph);
