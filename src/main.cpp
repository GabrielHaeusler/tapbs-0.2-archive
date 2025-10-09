/***************************************************************************
 *   Copyright (C) 2007 by Gernot Kieseritzky                              *
 *   gernotf@chemie.fu-berlin.de                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**
 *  @mainpage TAPBS
 * TAPBS calculates electrostatic energies needed for continuum 
 * electrostatics based pKa and redox calculations. It prepares 
 * the input files by Karlsberg 2.0 to predict protonation and redox 
 * patterns in a protein. Essentially, TAPBS is a front end to APBS 4.0, 
 * a stand-alone program developed by Nathan Baker at the 
 * University of Washington to solve the linearized Poisson-Boltzmann 
 * equation (LPBE) numerically, which tries to efficiently and conveniently 
 * automate the large number of LPBE calculations needed to compute Born 
 * solvation, background and pair wise energy terms that determine pKa values. 
 * TAPBS supports local conformational changes, multiple charged states, protein 
 * regions of different dielectric constant and implicit membrane slabs.
 *
 *  @author Gernot Kieseritzky
 */

/**
 * @file main.cpp This small file contains nothing but the ::main function.
 */

#include "titration.h"
#include <iostream>
#include <istream>

using namespace std;

int main(int argc, char *argv[])
{
 cout << "This is TAPBS release 0.2." << endl;
 cout << "Karlsberg was written 2006-2007 by Gernot Kieseritzky" << endl;
 cout << "(gernotf@chemie.fu-berlin.de)." << endl;
 cout << "and modified by Ilkay Sakalli in 2011." << endl;
 cout << "(sakalli@chemie.fu-berlin.de)" << endl;
 cout << "TAPBS comes with ABSOLUTELY NO WARRANTY. It is free software, and you" << endl;
 cout << "are welcome to redistribute it under certain conditions. See file COPYING" << endl;
 cout << "being part of this program package for details." << endl << endl;

 try
 {
  // Construct parser and parse from standard input
  InputDescription inputDescription;
  InputLine("legacy", "", false).addToInputDescription(inputDescription);
  InputLine("gunner", "", false).addToInputDescription(inputDescription);
  InputLine("test", "", false).addToInputDescription(inputDescription);
  APBSparms::getInputDescription(inputDescription);
  Titration::getInputDescription(inputDescription);

//  ifstream fin(argv[1]);
//  InputParser parser(inputDescription,  fin);
  InputParser parser(inputDescription,  cin);

  // run Titration
  if ( parser.hasInputLine("legacy") )
  {
   cout << "Starting 'Legacy' titration ... " << endl;
   Legacy titration(parser);
	titration.run();
  }
  else if ( parser.hasInputLine("gunner") )
  {
   cout << "Starting 'MultiConformerGunner' titration ... " << endl;
   MultiConformerGunner titration(parser);
	titration.run();
  }
  else if ( parser.hasInputLine("test") )
  {
   cout << "Starting 'MultiConformerTest' titration ... " << endl;
   MultiConformerTest titration(parser);
	titration.run();
  }  
  else
  {
   cout << "Starting 'MultiConformer' titration ... " << endl;
   MultiConformer titration(parser);
   titration.run();
  }
 }
 catch (GeneralException& e)
 {
  cout << "Abnormal termination due to fatal error:" << endl;
  cout << e << endl;
  return e.getErrorCode();
 }
 catch (...)
 {
  cout << "Abnormal termination due to unexpected error." << endl;
  return 255;
 }

 return 1;
}
