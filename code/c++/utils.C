//////////////////////////////////////////////////////////////////////
// utils.C                                                          //
// Author:  Dario Ghersi                                            //
// Version: 20130902                                                //
// Goal:    part of the MapBind suite. The file contains some       //
//          utility functions used by the other programs            //
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

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

vector<string> getFileToProcess(char *listOfFiles)
{
  // store the names of the files to process

  vector<string> fileNames;

  // open the file with the list of files to process
  fstream infile;
  infile.open(listOfFiles, fstream::in);

  // store the names of the files to process
  if (infile.fail()) {
    cerr << "Cannot open " << listOfFiles << endl;
    exit(1);
  }

  string line;
  while (getline(infile, line)) {
    fileNames.push_back(line);
  }

  infile.close();

  return fileNames;
}

//////////////////////////////////////////////////////////////////////

bool replace(string& str, const string& from, const string& to) {
  // replace a substring with another

  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;

  str.replace(start_pos, from.length(), to);
  return true;
}

//////////////////////////////////////////////////////////////////////

string basename(string s)
{
  // return the directory portion of a path name

  size_t pos = 0;
  string token;
  const string delimiter = "/";

  while ((pos = s.find(delimiter)) != string::npos) {
    token = s.substr(0, pos);
    s.erase(0, pos + delimiter.length());
  }

  return s;
}

//////////////////////////////////////////////////////////////////////

string dirname(string s)
{
  // return the file name portion of a path name
  
  size_t pos = 0;
  string token;
  const string delimiter = "/";

  while ((pos = s.find(delimiter)) != string::npos) {
    token = s.substr(0, pos);
    s.erase(pos, pos + delimiter.length());
  }

  return s;
}
