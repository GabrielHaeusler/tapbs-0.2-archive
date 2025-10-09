#include "titration.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <math.h>
#include <stdexcept>

void Titration::deleteTitratableResidues()
{
 for(iter_titratables iter=residues.begin(); iter!=residues.end(); iter++)
  delete (*iter);
 residues.clear();
}

void Titration::parseSITES(const string& filename)
{
 ifstream in(filename.c_str());
 string line, name, segname, stfile;
 int resid, lineNr=0;
 Titratable *residue;

 cout << "Parsing sites file " << filename << " ... " << endl;
 
 if (!in)
 {
  in.close();
  throw FileIOException(filename, "Titration::parseSITES");
 }
  
 while ( !in.eof() )
 {
  getline(in, line);
  vector<string> split_line = split(line);
  lineNr++;

  // ignore empty lines
  if(split_line.empty()) continue;

  if(split_line.size()==4)
  {
   /* 0 -> segname
    * 1 -> resid
    * 2 -> name
    * 3 -> ST file name
    */

   // parse segment name
   segname = split_line.at(0);   

   // parse resid
   if ( sscanf( split_line.at(1).c_str(), "%d", &resid )!=1 )
   {
    in.close();
    throw ParseException(lineNr, filename, "Tried to parse residue number but did not find integer value.");
   }

   // parse residue name
   name = split_line.at(2);

   // parse ST file name
   stfile = split_line.at(3);

   // create Titratable residue
   residue = new McTitratable(resid, name, segname);
   residue->parseST(stfile);
   residue->fillMissing(pqr);
  }
  else if(split_line.size()==3)
  {
   /* 0 -> resid
    * 1 -> name
    * 2 -> ST file name
    */

   // parse resid
   if ( sscanf( split_line.at(0).c_str(), "%d", &resid )!=1 )
   {
    in.close();
    throw ParseException(lineNr, filename, "Tried to parse residue number but did not find integer value.");
   }

   // parse residue name
   name = split_line.at(1);

   // parse ST file name
   stfile = split_line.at(2);

   // create Titratable residue
   residue = new Titratable(resid, name);
   residue->parseST(stfile);
   residue->fillMissing(pqr);
  }
  else
  {
   throw ParseException(lineNr, filename, "Each line has to comply to the structure 'SEGNAME RESID NAME ST-File' or 'RESID NAME ST-File'.");
  }  

  // add Titratable residue
  residues.push_back(residue);
 }
  
 cout << "... finished." << endl;
 
 in.close();
}

void Titration::setGridCenters(APBSparms *apbsparms, double *center[3], Structure *p, Structure *m)
{
 double px, py, pz;
 double mx, my, mz;
 
 p->get_center(&px, &py, &pz);
 m->get_center(&mx, &my, &mz);
 
 for(int i=0; i<apbsparms->get_nfsteps(); i++)
 {
  center[0][i] = (apbsparms->get_centering(i)==ON_PROTEIN)?px:mx;
  center[1][i] = (apbsparms->get_centering(i)==ON_PROTEIN)?py:my;
  center[2][i] = (apbsparms->get_centering(i)==ON_PROTEIN)?pz:mz;
  cout << "Center of grid in focusing step " << i+1 << " (x,y,z): " 
       << center[0][i] << ", " << center[1][i] << ", " << center[2][i] 
       << endl;
 }
}

void Titration::setSoluteDielectricConstants(InputParser& parser)
{
 double pdie = apbsparms->get_pdie();
 
 // first set to uniform dielectric
 while( pqr->hasNext() )
  pqr->getNext().set_soluteDiel(pdie);

 // now overwrite with user-specified atom-wise inhomogenities
 parser.setKeyword("diel_resid");
 while( parser.hasNextInputLine() )
 {
  InputLine inputLine = parser.getNextInputLine();
  int resid    = (int)inputLine.getArgument(1).getLongValue();
  double adiel = inputLine.getArgument(2).getDoubleValue();
  pqr->select_resid(resid);
  pqr->set_soluteDiel_in_selection(adiel);
  cout << "Set solute dielectric constant of atom " << endl;
  pqr->print_pqr_selection(cout);
  cout << "to a value of " << adiel << "." << endl;
  pqr->select_clear();
 }

 parser.setKeyword("diel_segres");
 while( parser.hasNextInputLine() )
 {
  InputLine inputLine = parser.getNextInputLine();
  string segid = inputLine.getArgument(1).getStringValue();
  int resid    = (int)inputLine.getArgument(2).getLongValue();
  double adiel = inputLine.getArgument(3).getDoubleValue();
  pqr->select_resid(segid, resid);
  pqr->set_soluteDiel_in_selection(adiel);
  cout << "Set solute dielectric constant of atom " << endl;
  pqr->print_pqr_selection(cout);
  cout << "to a value of " << adiel << "." << endl;
  pqr->select_clear();
 }
 
 parser.setKeyword("diel_atomno");
 while( parser.hasNextInputLine() )
 {
  InputLine inputLine = parser.getNextInputLine();
  int atomno   = (int)inputLine.getArgument(1).getLongValue();
  double adiel = inputLine.getArgument(2).getDoubleValue();
  pqr->select_atomno(atomno);
  pqr->set_soluteDiel_in_selection(adiel);
  cout << "Set solute dielectric constant of atom " << endl;
  pqr->print_pqr_selection(cout);
  cout << "to a value of " << adiel << "." << endl;
  pqr->select_clear();
 }
}

bool Titration::setupMembrane(InputParser& parser)
{
 if( parser.hasInputLine("membrane") )
 {
  double xs   = parser.getInputLine("membrane").getArgument(1).getDoubleValue();
  double ys   = parser.getInputLine("membrane").getArgument(2).getDoubleValue();
  double zs   = parser.getInputLine("membrane").getArgument(3).getDoubleValue();
  double xm   = parser.getInputLine("membrane").getArgument(4).getDoubleValue();
  double ym   = parser.getInputLine("membrane").getArgument(5).getDoubleValue();
  double zm   = parser.getInputLine("membrane").getArgument(6).getDoubleValue();
  Vector3D vs(xs, ys, zs);    /* support vector of membrane plane */
  Vector3D vm(xm, ym, zm);    /* vector describing orientation and width of membrane slab */
  Vector3D za(0.0, 0.0, 1.0); /* z-axis */
  
  /* rotate solute and membrane to co-align with z-axis */
  double deg  = vm.angle(za);
  Vector3D ax = vm.cross(za); /* axis of rotation */
  ax = 1/(ax.abs())*ax;       /* normalization */
  vs.rotate(ax, deg);
  vm.rotate(ax, deg);
  pqr->rotate(ax, deg);
  
  /* get membrane slab dimensions */
  zmin = vs.get_z_value();
  zmax = (vs+vm).get_z_value();
  cout << "Membrane dimensions: zmin = " << zmin << ", zmax = " << zmax << endl;
  
  /* membrane dielectric constant */
  epsm = parser.getInputLine("membrane").getArgument(7).getDoubleValue();
  cout << "Membrane dielectric constant = " << epsm << endl;

  /* okay! */
  return true;
 }

 /* no membrane slab specified */
 zmin = zmax = epsm = 0.0;
 return false;
}

/* public methods */

Titration::Titration(InputParser& parser)
{
 /* get current system time */
 startTimer();
 
 /* init parameters */
 apbsparms    = new APBSparms(parser);
 pqrfile      = parser.getInputLine("pqr").getArgument(1).getStringValue();
 sitesfile    = parser.getInputLine("sites").getArgument(1).getStringValue();
 jobname      = parser.getInputLine("output").getArgument(1).getStringValue();
 dummy        = parser.hasInputLine("dummy");
 symcheck     = parser.hasInputLine("symmetry_check");
 writemap     = parser.hasInputLine("writemap");

 stateFileExtension = ".st";
 potenFileExtension = ".potat";
 pkintFileExtension = ".pkint";
 wmatrFileExtension = ".g";

 /* init Structure */
 pqr = new Structure();
 pqr->fromPQR(pqrfile);
 setSoluteDielectricConstants(parser);

 /* deal with membrane slab */
 membrane = setupMembrane(parser);

 /* parse titrating sites */
 parseSITES(sitesfile);
}

/* ---------------------------------- Legacy -------------------------------------- */

void Legacy::calc_pka(const string& filestem)
{
 iter_titratables iter;
 Titratable *residue;
 int resid;
 Structure *mback, *mref, *mstat;
 
 double pka, gborn, gback;
 
 cout << "Calculating intrinsic pKAs ... ";
 cout.flush();
 
 // setup output file
 string filename(jobname+pkintFileExtension);
 ofstream out( filename.c_str() );
 out.setf(ios_base::scientific);
 out.precision(6);

 // setup log
 ostringstream log;
 log.setf(ios_base::scientific);
 log.precision(6);
 if (!out)
 {
  out.close();
  throw FileIOException(filename, "Legacy::calc_pka");
 }

 /* all residues should be in reference state in pqr structure by now */

 // run over all titratable residues
 for(iter=residues.begin(); iter!=residues.end(); iter++) // for each residue
 {
  residue = *iter;
  resid   = residue->getResID();
  
  // get potentials
  PotentialsOfResidue p_container;
  p_container.readFromFile( jobname+"."+residue->getName()+filestem+potenFileExtension );

  // get model reference charges
  residue->select(pqr);
  mref   = pqr->createModelFromSelection();
  mref->clear_all_charges();
  residue->switchToReferenceState(mref); /* non-titrating charges zero */
  
  // get charges of non-reference state
  mstat  = new Structure(*mref);
  mstat->clear_all_charges();
  residue->switchToState(1, mstat);

  // background charges
  residue->clearChargesOfState(0, pqr);
  residue->select(pqr);
  mback = pqr->createModelFromSelection(); /* N- and C-terminus charges included! */
  
  // MEAD legacy energy scale
  int sign = residue->isAnionic(1)?-1:1;

  // born
  gborn = sign * p_container.calcBorn(1, mstat, mref) / LN10;
  
  //back
  gback = sign * p_container.calcBack(1, pqr, mback) / LN10;
  
  //intrinsic pKa
  pka = residue->getEnergy(1) * KJMOL_TO_PK / T - gborn + gback;

  log << residue->getName() << ": R->" << (residue->isAnionic(1)?"A":"C") << endl;
  log << "Gborn: " << gborn << endl;
  log << "Gback: " << gback << endl;
  log << "intrinsic pka = " 
       << residue->getEnergy(1) * KJMOL_TO_PK / T
	    << " - " << gborn
	    << " - " << gback
		 << " = " << pka 
		 << endl;
 
  // write pKa and transition identifier to pkint file
  // identifier is A for anionic and C for cationic
  out  << pka << " " << (residue->isAnionic(1)?"A":"C") << " ";
 
  // write residue name to pkint file
  out << residue->getName() << endl;

  // add reference charges of residue to protein background
  residue->setChargesOfState(0, pqr);
  
  // clean up
  delete mback;
  delete mref;
  delete mstat;
 }

 out.close();

 cout << "finished." << endl;
 cout << log.str();
}

void Legacy::calc_wmatrix(const string& filestem)
{
 Structure *pref, *pstat;

 double G;

 cout << "Calculating W matrix ... ";
 cout.flush();
 
 try
 {
 
 // conversion unit
 const double CONVERT = T * MEADUNITS;

 // create wmatrix
 for(unsigned int i = 0; i<residues.size(); i++)
 {
  wmatrix.push_back(vector<double>());
  for(unsigned int j = 0; j<residues.size(); j++) 
   wmatrix.at(i).push_back(0.0);
 }
 
 // calc matrix
 for(unsigned int i = 0; i<residues.size(); i++)
 {
  // read potentials of residue 1
  PotentialsOfResidue p_container;
  p_container.readFromFile( jobname+"."+residues[i]->getName()+filestem+potenFileExtension );

  for(unsigned int j = 0; j<residues.size(); j++)
  {
	// diagonal elements to be zero
	if ( residues[i]->getName()==residues[j]->getName() )
	{
	 wmatrix[i][j] = 0.0;
	 continue;
	}

   /* 
    legacy mead energy scaling factor:
	 
    energy = scale * [Q(p)-Q(d)] * [P(p)-P(d)]
    anionic state 1:  reference is protonated form
    cationic state 1: reference is deprotonated form
   */
   int scale1 = residues[i]->isAnionic(1)?-1:1;
   int scale2 = residues[j]->isAnionic(1)?-1:1;
   int scale  = scale1*scale2;

	// charges of residue 2
	residues[j]->select(pqr);
	pref  = pqr->fromSelection();
   pstat = pqr->fromSelection();
	pqr->select_clear();

	// reference state
	pref->clear_all_charges();
	residues[j]->switchToReferenceState(pref);

   // charged state
   pstat->clear_all_charges();
   residues[j]->switchToState(1, pstat);

	// Wmv
	G = p_container.calcWmv(1, scale, pstat, pref);
	G = G * CONVERT; // kt->e^2/A
   wmatrix[i][j] = G;
   
	// clean up 
	delete pref;
   delete pstat;
  }
 }
 
 }
 catch(out_of_range& e)
 {
  throw OutOfRangeException("Legacy::calc_wmatrix");
 }
 
 cout << "finished." << endl;
}

void Legacy::symmetrize_wmatrix()
{
 double dev=0.0, w1=0.0, w2=0.0, ave=0.0;
 
 cout << "Symmetrizing W matrix ... ";
 cout.flush();

 // setup output stream for symmetry check
 ostringstream log;
 log.setf(ios_base::scientific);
 log.precision(6);
 
 for(unsigned int i=0; i<residues.size(); i++)
 {
  for(unsigned int j=i+1; (j>i) && (j<residues.size()); j++)
  {
   w1    = wmatrix[i][j];
	w2    = wmatrix[j][i];
	dev   = (w1-w2);
	ave   = 0.5*(w1+w2);
	wmatrix[i][j] = ave;
	wmatrix[j][i] = ave;
   // write symmetry check to log file
   log << residues[i]->getName() << "/";
   log << residues[j]->getName() << ": ";
   log << "(" << w1 << "," << w2 << "), diff = " << dev;
   log << ", ave = " << ave << endl;
  }
 }

 cout << "finished." << endl;

 if (symcheck) 
  cout << log.str();
}

void Legacy::write_wmatrix()
{
 string filename(jobname+wmatrFileExtension);
 ofstream out( filename.c_str() );
 out.setf(ios_base::scientific);
 out.precision(6);

 cout << "Writing W matrix into '" << filename << "' ... ";
 cout.flush();

 if (!out)
 {
  out.close();
  throw FileIOException(filename, "Legacy::write_wmatrix");
 }
 
 // run over all matrix elements
 for(unsigned int i=0; i<residues.size(); i++)
 {
  for(unsigned int j=0; j<residues.size(); j++)
  {
   out.width(4);
	out << i+1 << " ";
	out.width(4);
   out << j+1 << "    ";
   out << wmatrix[i][j] << endl;
  }
 }
 
 out.close();
 
 cout << "finished." << endl;
}

void Legacy::calc_protein_potentials(Structure *proteinVolume, const string& filestem)
{
 Titratable *residue;
 MapController *mapcontroller;
 Structure *site;
 double *center[3]; /* grid centers */
 int firstResidue = 1; /* flag */

 // init arrays that store the grid center in each focusing step
 center[0] = new double[apbsparms->get_nfsteps()];
 center[1] = new double[apbsparms->get_nfsteps()];
 center[2] = new double[apbsparms->get_nfsteps()];
 
 // initialize protein mapcontroller
 if(membrane)
  mapcontroller = new MembraneMapController(apbsparms, epsm, zmin, zmax);
 else
  mapcontroller = new MapController(apbsparms);
 
 for(iter_titratables iter=residues.begin(); iter!=residues.end(); iter++)
 {
  // useful local variables
  residue = *iter;
  PotentialsOfResidue potentials;
  string potentialFileName = jobname+"."+residue->getName()+filestem+potenFileExtension;
  
  // mark residue name in log file
  cout << residue->getName() << endl;
 
  // if present read existing potential file
  potentials.readFromFile(potentialFileName);

  // skip calculation if state is present in potential file
  if ( potentials.nrOfProteinStates()==residue->nrOfStates() )
  {
   cout << "Protein potentials already calculated for this residue." << endl;
   continue;
  }

  // center grid in each focusing step on protein proteinVolume center or site
  residue->select(proteinVolume);
  site = proteinVolume->fromSelection();
  proteinVolume->select_clear();
  setGridCenters( apbsparms, center, proteinVolume, site );

  // compute all possible states: 0 = reference
  for(int state=0; state<residue->nrOfStates(); state++)
  {
   // add partial charges to protein volume
	residue->setChargesOfState(state, proteinVolume);

   cout << "Starting protein calculation for state " << state+1 << "." << endl;
	residue->select(proteinVolume);
	proteinVolume->print_pqr_selection(cout);
   proteinVolume->select_clear();

   // debug
   if(writemap) 
	{
	 ostringstream ost;
	 ost << residue->getName() << ".state" << state+1 << ".protein";
	 mapcontroller->writeMaps(ost.str());
	}

   // call APBS
   APBS_runner apbs(1, proteinVolume, apbsparms);
	if (dummy)
	 apbs.setupDummy(center[0], center[1], center[2], mapcontroller);
   else
    apbs.setup(center[0], center[1], center[2], mapcontroller);
   apbs.run();

	// get electrostatic potential
	potentials.addProteinStatePotential( apbs.get_potential() );

	// remove charges from protein volume
	residue->clearChargesOfState(state, proteinVolume);
  }

  // store potentials in file for later use
  potentials.writeToFile(potentialFileName);

  // clean up
  delete site;
 }

 // clean up
 delete[] center[0];
 delete[] center[1];
 delete[] center[2];
 delete mapcontroller;
}

void Legacy::calc_model_potentials(Structure *proteinVolume, const string& filestem)
{
 Titratable *residue;
 MapController *mapcontroller;
 Structure *model, *site;
 APBSparms_model *apbsparms_model;
 double *center[3]; /* grid centers */

 // get model-specific APBS parameters
 apbsparms_model = new APBSparms_model(*apbsparms);

 // init arrays that store the grid center for each focusing step
 center[0] = new double[apbsparms_model->get_nfsteps()];
 center[1] = new double[apbsparms_model->get_nfsteps()];
 center[2] = new double[apbsparms_model->get_nfsteps()];

 // prepare accessibility maps
 mapcontroller = new MapController(apbsparms_model);

 for(iter_titratables iter=residues.begin(); iter!=residues.end(); iter++)
 {
  // useful local variables
  residue = *iter;
  PotentialsOfResidue potentials;
  string potentialFileName = jobname+"."+residue->getName()+filestem+potenFileExtension;

  // mark residue name in log file
  cout << residue->getName() << endl;
 
  // if present read existing potential file
  potentials.readFromFile(potentialFileName);

  // skip calculation if state is present in potential file
  if ( potentials.nrOfModelStates()==residue->nrOfStates() )
  {
   cout << "Model potentials already calculated for this residue." << endl;
   continue;
  }

  // generate model compound
  residue->select(pqr);
  model = pqr->createModelFromSelection();
  model->clear_all_charges();  

  // center grid in each focusing step on protein volume center or site
  residue->select(proteinVolume);
  site = proteinVolume->fromSelection();
  proteinVolume->select_clear();
  setGridCenters(apbsparms_model, center, proteinVolume, site);

  // compute all possible states: 0 = reference, 1 = ionized state
  for(int state=0; state<residue->nrOfStates(); state++)
  {
   // setup partial charges
	residue->switchToState(state, model);

   cout << "Starting model calculation for state " << state+1 << "." << endl;
   model->print_pqr(cout);

   // debug
   if(writemap) 
	{
	 ostringstream ost;
	 ost << residue->getName() << ".state" << state+1 << ".model";
	 mapcontroller->writeMaps(ost.str());
	}

   APBS_runner apbs(1, model, apbsparms_model);

	if (dummy)
    apbs.setupDummy(center[0], center[1], center[2], mapcontroller);
   else
    apbs.setup(center[0], center[1], center[2], mapcontroller);

	apbs.run();

   // save elecrostatic potential
	potentials.addModelStatePotential( apbs.get_potential() );

	// clean up	
	residue->switchToReferenceState(model);
	residue->clearChargesOfState(0, model);
  }

  // store potentials in file for later use
  potentials.writeToFile(potentialFileName);

  // clean up
  delete model;
  delete site;
 }
 
 // clean up
 delete apbsparms_model;
 delete[] center[0];
 delete[] center[1];
 delete[] center[2];
 delete mapcontroller;
}

Structure* Legacy::createProteinDielectricVolume()
{
 Structure *volume;
 
 // switch residues in pqr structure into reference state
 for(iter_titratables iter=residues.begin(); iter!=residues.end(); iter++)
  (**iter).switchToReferenceState(pqr);
 
 pqr->select_all();
 volume = pqr->fromSelection();
 volume->clear_all_charges();
 pqr->select_clear();

 pqr->writeToFile_pqr(jobname+".reference.pqr");
 
 return volume;
}

/* ------------------------------------- MultiConformer -------------------------------------- */

void MultiConformer::calc_pka(const string& filestem)
{
 iter_titratables iter;
 Titratable *residue;
 int resid;
 Structure *mback, *mref, *mstat;

 double pka, gborn, gback;
 
 cout << "Calculating intrinsic pKAs ... ";
 cout.flush();

 // setup output file
 string filename(jobname+pkintFileExtension);
 ofstream out( filename.c_str() );
 out.setf(ios_base::scientific);
 out.precision(6);

 // setup log
 ostringstream log;
 log.setf(ios_base::scientific);
 log.precision(6);
 if (!out)
 {
  out.close();
  throw FileIOException(filename, "MultiConformer::calc_pka");
 }

 /* all residues should be in its reference state after calling of createProteinDielectric */
 
 // run over all titratable residues
 for(iter=residues.begin(); iter!=residues.end(); iter++) // for each residue
 {
  residue = *iter;
  resid   = residue->getResID();
  
  // get potentials
  PotentialsOfResidue p_container;
  p_container.readFromFile( jobname+"."+residue->getName()+filestem+potenFileExtension );

  // get model reference charges
  residue->select(pqr);
  mref = pqr->createModelFromSelection();
  mref->clear_all_charges();
  residue->switchToReferenceState(mref); /* non-titrating charges zero */

  // prepare model state charges
  mstat  = new Structure(*mref);
  mstat->clear_all_charges();

  // background charges
  residue->clearChargesOfState(0, pqr);
  residue->select(pqr);
  mback = pqr->createModelFromSelection(); /* N- and C-terminus charges included! */

  // write reference energy [kJ/mol] to pkint file
  out  << 0.0 << " " << residue->getID(0) << " ";

  // run over all states
  for(int state=1; state<residue->nrOfStates(); state++)
  {
   // get charges of current state
   residue->switchToState(state, mstat);

   // born
   gborn = p_container.calcBorn(state, mstat, mref) * R * T;
  
   // back
   gback = p_container.calcBack(state, pqr, mback) * R * T;

   // intrinsic pKa
   pka = residue->getEnergy(state) * T + gborn + gback; /* kJ/mol */

   // write energies to log file for control
   log << residue->getName() << ": R->" << residue->getID(state) << endl;
   log << "Gborn: " << gborn << endl;
   log << "Gback: " << gback << endl;
   log << "intrinsic pka = " 
        << residue->getEnergy(state) * T
        << " + " << gborn
        << " + " << gback
        << " = " << pka << endl;

   // write pkint [kJ/mol] and transition identifier to pkint file
   out << pka << " " << residue->getID(state) << " ";
	
	// clean up
   residue->switchToReferenceState(mstat);
	residue->clearChargesOfState(0, mstat);
  }
  
  // write residue name to pkint file
  out << residue->getName() << endl;

  // add reference charges of residue to protein background
  residue->setChargesOfState(0, pqr);
  
  // clean up
  delete mback;
  delete mref;
  delete mstat;
 }

 out.close();

 cout << "finished." << endl;
 cout << log.str();
}

void MultiConformer::calc_wmatrix(const string& filestem)
{
 Structure *pref, *pstat;
 
 double G;

 cout << "Calculating W matrix ... ";
 cout.flush();
 
 try
 {
 
 // calc conversion unit
 const double CONVERT = T * MEADUNITS;

 // create wmatrix
 for(unsigned int i = 0; i<residues.size(); i++)
 {
  wmatrix.push_back( vector< vector< vector<double> > >() );
  for(unsigned int j = 0; j<residues.size(); j++)
  {
   wmatrix.at(i).push_back( vector< vector<double> >() );
	/* there are no matrix entries for reference states
	 that's why we take only the number of states-1 */
	for(unsigned int state1 = 0; state1<residues[i]->nrOfStates()-1; state1++)
	{
	 wmatrix.at(i).at(j).push_back( vector<double>() );
	 for(unsigned int state2 = 0; state2<residues[j]->nrOfStates()-1; state2++)
	 {
	  wmatrix.at(i).at(j).at(state1).push_back(0.0);
	 }
	}
  }
 }

 // calc matrix
 for(unsigned int i = 0; i<residues.size(); i++)
 {
  // read potentials of residue 1
  PotentialsOfResidue p_container;
  p_container.readFromFile( jobname+"."+residues[i]->getName()+filestem+potenFileExtension );
  
  for(unsigned int j = 0; j<residues.size(); j++)
  {
	// diagonal elements to be zero
	if ( residues[i]->getName()==residues[j]->getName() ) 
	{
    for(unsigned int state1 = 1; state1<residues[i]->nrOfStates(); state1++)
    {
     for(unsigned int state2 = 1; state2<residues[j]->nrOfStates(); state2++)
     {
      wmatrix[i][j][state1-1][state2-1] = 0.0;
	  }
	 }
	 continue;
	}
	
	// charges of residue 2
	residues[j]->select(pqr);
	pref  = pqr->fromSelection();
   pstat = pqr->fromSelection();
	pqr->select_clear();

	// reference state
	pref->clear_all_charges();
	residues[j]->switchToReferenceState(pref);

   // charged state
   pstat->clear_all_charges();

   // run over all non-reference states
   for(unsigned int state1 = 1; state1<residues[i]->nrOfStates(); state1++)
	{
    for(unsigned int state2 = 1; state2<residues[j]->nrOfStates(); state2++)
 	 {
     // get charge of current state in residue 2
	  residues[j]->switchToState(state2, pstat);

	  // Wmv
	  G = p_container.calcWmv(state1, 1, pstat, pref);
	  G = G * CONVERT; // kt->e^2/A

     // update matrix
	  wmatrix[i][j][state1-1][state2-1] = G;
	  
	  // clean up
     residues[j]->switchToReferenceState(pstat);
	  residues[j]->clearChargesOfState(0, pstat);
	 }
	}
	
	// clean up
	delete pref;
   delete pstat;
  }
 }

 }
 catch(out_of_range& e)
 {
  throw OutOfRangeException("Titration::calc_wmatrix");
 }

 cout << "finished." << endl;
}

void MultiConformer::symmetrize_wmatrix()
{
 double dev=0.0, w1=0.0, w2=0.0, ave=0.0;
 
 cout << "Symmetrizing W matrix ... ";
 cout.flush();

 // setup output stream for symmetry check
 ostringstream log;
 log.setf(ios_base::scientific);
 log.precision(6);

 try
 {

 for(unsigned int i=0; i<residues.size(); i++)
 {
  for(unsigned int j=i+1; (j>i) && (j<residues.size()); j++)
  {
   for(unsigned int state1 = 1; state1<residues[i]->nrOfStates(); state1++)
	{
	 for(unsigned int state2 = 1; state2<residues[j]->nrOfStates(); state2++)
	 {
	  w1  = wmatrix[i][j][state1-1][state2-1];
	  w2  = wmatrix[j][i][state2-1][state1-1];
     dev = (w1-w2);
     ave = 0.5*(w1+w2);
     wmatrix[i][j][state1-1][state2-1] = ave;
     wmatrix[j][i][state2-1][state1-1] = ave;
	  // write symmetry check to log file
	  log << residues[i]->getName() << "[" << state1+1 << "]/";
	  log << residues[j]->getName() << "[" << state2+1 << "]: ";
	  log << "(" << w1 << "," << w2 << "), diff = " << dev;
	  log << ", ave = " << ave << endl;
	 }
	}
  }
 }

 }
 catch(out_of_range& e)
 {
  throw OutOfRangeException("Titration::symmetrize_wmatrix");
 }
 
 cout << "finished." << endl;

 if(symcheck)
  cout << log.str();
}

void MultiConformer::write_wmatrix()
{
 // setup output file
 string filename(jobname+wmatrFileExtension);
 ofstream out( filename.c_str() );
 out.setf(ios_base::scientific);
 out.precision(6);

 cout << "Writing W matrix into '" << filename << "' ... ";
 cout.flush();

 if (!out)
 {
  out.close();
  throw FileIOException(filename, "Titration::write_wmatrix");
 }

 try
 {
 
 // run over all matrix elements
 for(unsigned int i=0; i<residues.size(); i++)
 {
  for(unsigned int j=0; j<residues.size(); j++)
  {
   for(unsigned int state1 = 1; state1<residues[i]->nrOfStates(); state1++)
	{
	 for(unsigned int state2 = 1; state2<residues[j]->nrOfStates(); state2++)
	 {
     out.width(4);
     out << i+1 << " " << state1+1 << " ";
     out.width(4);
     out << j+1 << " " << state2+1 << "    ";
	  out << wmatrix[i][j][state1-1][state2-1] << endl;
	 }
	}
  }
 }

 }
 catch(out_of_range& e)
 {
  out.close();
  throw OutOfRangeException("Titration::write_wmatrix");
 }

 out.close();

 cout << "finished." << endl;
}

Structure* MultiConformer::createProteinDielectricVolume()
{
 int resid;
 Titratable *residue;
 Structure *volume, *local, buffer;
 
 cout << "Creating combined dielectric volume ... ";
 cout.flush();

 // first set pqr structure into reference state
 for(iter_titratables iter=residues.begin(); iter!=residues.end(); iter++)
  (**iter).switchToReferenceState(pqr);
 
 // create a copy of PQR structure
 pqr->select_all();
 volume = pqr->fromSelection();
 pqr->select_clear();
 
 // now collect all conformers into buffer
 for(iter_titratables iter=residues.begin(); iter!=residues.end(); iter++)
 {
  residue = *iter;
  resid = residue->getResID();
  
  // local copy
  residue->select(pqr);
  local = pqr->fromSelection();
  pqr->select_clear();

  // loop over all states>1, since volume already contains reference state
  for(int state=1; state<residue->nrOfStates(); state++)
  {
   residue->switchToState(state, local);
	buffer.append(*local);
	residue->switchToReferenceState(local);
	/*residue->clearChargesOfState(0, local);*/
  }

  // clean up  
  delete local;
 }

 // append to protein dielectric volume
 volume->appendAtUnknownCoordinates(buffer);
 volume->update();

 // remove charges from volume 
 volume->clear_all_charges();

 volume->writeToFile_pqr(jobname+".volume.pqr");
 pqr->writeToFile_pqr(jobname+".reference.pqr");

 cout << "finished." << endl;

 return volume;
}

/* ----------------------------------- MultiConformerGunner -------------------------------------- */

void MultiConformerGunner::calc_intrinsic_protein_potentials(Structure *proteinVolume, const string& filestem)
{
 Titratable *residue;
 MapController *mapcontroller;
 Structure *site, localCopy;
 double *center[3];        /* grid centers */

 // init arrays that store the grid center for each focusing step
 center[0] = new double[apbsparms->get_nfsteps()];
 center[1] = new double[apbsparms->get_nfsteps()];
 center[2] = new double[apbsparms->get_nfsteps()];

 // local copy of PQR, every residue should be in its reference state conformation after createProteinDielectric
 localCopy = *pqr;
 localCopy.clear_all_charges();

 // prepare accessiblity maps
 if(membrane)
  mapcontroller = new MembraneMapController(apbsparms, epsm, zmin, zmax);
 else
  mapcontroller = new MapController(apbsparms);
  
 for(iter_titratables iter=residues.begin(); iter!=residues.end(); iter++)
 {
  // useful local variables
  residue = *iter;
  PotentialsOfResidue potentials;
  string potentialFileName = jobname+"."+residue->getName()+filestem+potenFileExtension;

  // mark residue name in log file
  cout << residue->getName() << endl;
 
  // if present read existing potential file
  potentials.readFromFile(potentialFileName);

  // skip calculation if state is present in potential file
  if ( potentials.nrOfProteinStates()==residue->nrOfStates() )
  {
   cout << "Intrinsic protein potentials already calculated for this residue." << endl;
   continue;
  }

  // center grid in each focusing step on protein volume center or site
  residue->select(proteinVolume);
  site = proteinVolume->fromSelection();
  proteinVolume->select_clear();
  setGridCenters(apbsparms, center, proteinVolume, site);

  // compute all possible states: 0 = reference, 1 = ionized state
  for(int state=0; state<residue->nrOfStates(); state++)
  {
   // setup partial charges
	residue->switchToState(state, &localCopy);

   cout << "Computing intrinsic protein potential for state " << state+1 << "." << endl;
	residue->select(&localCopy);
   localCopy.print_pqr_selection(cout);
	localCopy.select_clear();

   // debug
   if(writemap) 
	{
	 ostringstream ost;
	 ost << residue->getName() << ".state" << state+1 << ".protein_intrinsic";
	 mapcontroller->writeMaps(ost.str());
	}

   APBS_runner apbs(1, &localCopy, apbsparms);

	if (dummy)
	 apbs.setupDummy(center[0], center[1], center[2], mapcontroller);
   else
    apbs.setup(center[0], center[1], center[2], mapcontroller);

   apbs.run();

   // save elecrostatic potential
	potentials.addProteinStatePotential( apbs.get_potential() );
	
	// clean up
   residue->switchToReferenceState(&localCopy);
	residue->clearChargesOfState(0, &localCopy);
  }

  // store potentials in file for later use
  potentials.writeToFile(potentialFileName);

  // clean up
  delete site;
 }
 
 // clean up
 delete[] center[0];
 delete[] center[1];
 delete[] center[2];
 delete mapcontroller;
}

void MultiConformerGunner::calc_pairwise_protein_potentials(Structure *proteinVolume, const string& filestem)
{
 Titratable *residue;
 MapController *mapcontroller;
 Structure *site, localCopy;
 double *center[3];        /* grid centers */

 // init arrays that store the grid center for each focusing step
 center[0] = new double[apbsparms->get_nfsteps()];
 center[1] = new double[apbsparms->get_nfsteps()];
 center[2] = new double[apbsparms->get_nfsteps()];

 // prepare accessibility maps
 if(membrane)
  mapcontroller = new MembraneMapController(apbsparms, epsm, zmin, zmax);
 else
  mapcontroller = new MapController(apbsparms);

 for(iter_titratables iter=residues.begin(); iter!=residues.end(); iter++)
 {
  // useful local variables
  residue = *iter;
  PotentialsOfResidue potentials;
  string potentialFileName = jobname+"."+residue->getName()+filestem+potenFileExtension;

  // mark residue name in log file
  cout << residue->getName() << endl;
 
  // if present read existing potential file
  potentials.readFromFile(potentialFileName);

  // skip calculation if state is present in potential file
  if ( potentials.nrOfProteinStates()==residue->nrOfStates() )
  {
   cout << "Pairwise protein potentials already calculated for this residue." << endl;
   continue;
  }

  // center grid in each focusing step on protein volume center or site
  residue->select(proteinVolume);
  site = proteinVolume->fromSelection();
  proteinVolume->select_clear();
  setGridCenters(apbsparms, center, proteinVolume, site);

  // remove all conformers for this residue from combined volume
  localCopy = *proteinVolume;
  for(int state=0; state<residue->nrOfStates(); state++)
   residue->selectCoordinatesOfState(state, &localCopy);
  localCopy.delete_atoms_in_selection();
  
  // add reference state to combined protein volume
  residue->select(pqr);
  Structure *temp = pqr->fromSelection();
  pqr->select_clear();
  temp->clear_all_charges();
  localCopy.appendAtUnknownCoordinates( *temp );
  localCopy.update();
  delete temp;

  // compute all possible states: 0 = reference, 1 = ionized state
  for(int state=0; state<residue->nrOfStates(); state++)
  {
   // setup partial charges
	residue->switchToState(state, &localCopy);

   // debug
	ostringstream ostr;
	ostr << jobname << "." << residue->getName() << "." << state << ".volume.pqr";
   localCopy.writeToFile_pqr( ostr.str() );

   cout << "Computing pairwise protein potential for state " << state+1 << "." << endl;
   residue->select(&localCopy);
	localCopy.print_pqr_selection(cout);
	localCopy.select_clear();

   // debug
   if(writemap) 
	{
	 ostringstream ost;
	 ost << residue->getName() << ".state" << state+1 << ".protein_pairwise";
	 mapcontroller->writeMaps(ost.str());
	}

   APBS_runner apbs(1, &localCopy, apbsparms);

	if (dummy)
	 apbs.setupDummy(center[0], center[1], center[2], mapcontroller);
   else
    apbs.setup(center[0], center[1], center[2], mapcontroller);

   apbs.run();

   // save elecrostatic potential
	potentials.addProteinStatePotential( apbs.get_potential() );

	// clean up
   residue->switchToReferenceState(&localCopy);
	residue->clearChargesOfState(0, &localCopy);
  }

  // store potentials in file for later use
  potentials.writeToFile(potentialFileName);

  // clean up
  delete site;
 }
 
 // clean up
 delete[] center[0];
 delete[] center[1];
 delete[] center[2];
 delete mapcontroller;
}
