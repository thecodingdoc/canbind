//////////////////////////////////////////////////////////////////////
// clusterBindingSites.C                                            //
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

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <iterator>
#include <set>
#include "utils.h"
#include "clusterBindingSites.h"

using namespace boost;
using namespace std;

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

void clusterSites(string infileName, double jaccardThresh)
{
  // cluster the sites in the file based on the jaccard threshold

  // open the input file
  fstream infile;
  infile.open(infileName.c_str(), fstream::in);
  if (infile.fail()) {
    cerr << "Cannot open " << infileName << endl;
    return;
  }

  // store the sites
  vector <string> sites;
  string line;
  while (getline(infile, line)) {
    sites.push_back(line);
  }

  // if there is only one site, do nothing
  if (sites.size() == 1) {
    return;
  }

  // build a graph based on the jaccard similarity
  MyGraph g(sites.size());
  for (unsigned int i = 0; i < sites.size() - 1; i++) {

    // get the residues in the site
    set<string> siteA = getSite(sites[i]);

    for (unsigned int j = i + 1; j < sites.size(); j++) {
      // get the residues in the site
      set<string> siteB = getSite(sites[j]);

      // calculate the Jaccard coefficient
      vector<string> intersect;
      set_intersection(siteA.begin(), siteA.end(), siteB.begin(),
		       siteB.end(), std::back_inserter(intersect));

      vector<string> sunion;
      set_union(siteA.begin(), siteA.end(), siteB.begin(),
		siteB.end(), std::back_inserter(sunion));

      double jaccard = float(intersect.size()) / sunion.size();

      // put an edge if jaccard > threshold
      if (jaccard > jaccardThresh) {
	add_edge(i, j, g);
      }
    }
  }

  infile.close();

  // find the connected components
  vector<vector<int> > clusters = findClusters(g);

  // print the results
  string outfileName = infileName;
  replace(outfileName, ".binding", ".clusters");
  fstream outfile;
  outfile.open(outfileName.c_str(), fstream::out);

  for (unsigned int i = 0; i < clusters.size(); i++) {
    string toPrint;
    bool isFirst = true;
    for (unsigned int j = 0; j < clusters[i].size(); j++) {
      if (isFirst) {
	isFirst = false;
	toPrint += sites[clusters[i][j]];
      }
      else {
	toPrint += " - " + sites[clusters[i][j]];
      }
    }
    outfile << toPrint << endl;
  }
  outfile.close();

}

//////////////////////////////////////////////////////////////////////

vector<vector<int> > findClusters(MyGraph g)
{
  // find the connected components in the graph

  typedef adjacency_list<vecS, vecS, undirectedS> Graph;
  vector<int> component(num_vertices(g));
  
  int num = connected_components(g, &component[0]);

  // cluster the nodes
  vector<vector<int> > clusters(num);

  for (unsigned int i = 0; i < component.size(); i++) {
    clusters[component[i]].push_back(i);
  }

  return clusters;
}

//////////////////////////////////////////////////////////////////////

set<string> getSite(string line)
{
  // get the binding residues and store them in a set

  set<string> site;

  string field;
  stringstream linestream(line);

  // get the site
  unsigned int count = 0;
  while (count <= BINDING_SITE_RES) {
    getline(linestream, field, '\t');
    count++;
  }

  // get the individual residues
  stringstream linestream2(field);
  string residue;
  while (getline(linestream2, residue, ' ')) {
    site.insert(residue);
  }

  return site;
}

//////////////////////////////////////////////////////////////////////
// MAIN                                                             //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // parse the command-line arguments
  if (argc != 3) {
    cerr << "Usage: clusterBindingSites LIST_OF_FILES ...\n";
    cerr << "                           JACCARD_THRESHOLD\n";
    exit(1);
  }
  char *listOfFiles = argv[1];
  double jaccThreshold = atoi(argv[2]);

  // store the list of file names to process
  vector<string> fileNames = getFileToProcess(listOfFiles);

  // process each file
  for (vector<string>::iterator fileName = fileNames.begin();
       fileName != fileNames.end(); fileName++) {
    clusterSites(*fileName, jaccThreshold);
  }

  return 0;
}
