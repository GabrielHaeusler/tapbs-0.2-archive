#include <iostream>
#include <string>
#include <vector>

using namespace std;

/* natural constants & conversion factors */

const double LN10          = 2.3025850929940459;
const double NA            = 6.02214179E23;
const double R             = 0.008314510; /* kJ / (mol K) */
const double PK_TO_KJMOL   = R * LN10;    /* kJ / (mol K) */
const double KJMOL_TO_PK   = 1/PK_TO_KJMOL;
const double MEADUNITS     = 1.3806505E-23 * 4.336663E17; /* e^2/(A * K * mol) */

/* helpful methods */

vector<string> split(const string& str);
int str2int(const string& str);
double str2double(const string& str);

template<class T> int min_index( T a[], int len )
{
 T min = a[0];
 int i_min = 0;
 
 for(int i=0; i<len; i++)
 {
  if ( a[i]<min ) 
  {
   min = a[i];
	i_min = i;
  }
 }
 
 return i_min;
}

template<class T> int max_index( T a[], int len )
{
 T max = a[0];
 int i_max = 0;
 
 for(int i=0; i<len; i++)
 {
  if ( a[i]>max ) 
  {
   max = a[i];
	i_max = i;
  }
 }
 
 return i_max;
}
