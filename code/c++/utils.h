//////////////////////////////////////////////////////////////////////
// utils.h                                                          //
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

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

string basename(string);
string dirname(string);
vector<string> getFileToProcess(char *);
bool replace(string &, const string &, const string &);
