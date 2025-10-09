#include "utils.h"
#include <sstream>

vector<string> split(const string& str)
{
 /* C++ implementation of perl split(' ', str) */

 vector<string> cols;

 istringstream in(str);
 
 while( !in.eof() )
 {
  string temp;
  in >> temp;
  if ( temp.empty() ) continue;
  cols.push_back(temp);
 }
 
 return cols;
}

int str2int(const string& str)
{
 int val;
 
 istringstream in(str);
 
 in >> val;
 
 return val;
}

double str2double(const string& str)
{
 double val;
 
 istringstream in(str);
 
 in >> val;
 
 return val;
}
