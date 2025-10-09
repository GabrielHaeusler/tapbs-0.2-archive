#include "apbswrapper.h"
#include <stack>
#include <sstream>
#include <stdexcept>

APBSparms::APBSparms(InputParser& parser)
{
 // initialize containers
 dime.push_back( vector< vector<int> >() );       // protein
 dime.push_back( vector< vector<int> >() );       // model
 spacing.push_back( vector< vector<double> >() ); // protein
 spacing.push_back( vector< vector<double> >() ); // model
 centering.push_back( vector<string>() );         // protein
 centering.push_back( vector<string>() );         // model

 // get user defined parameters
 try
 {
  init(parser);
 }
 catch (out_of_range e)
 {
  throw OutOfRangeException("APBSparms::init");
 }
}

APBSparms::APBSparms(const APBSparms& apbsparms)
{
 this->rparm	  = apbsparms.rparm;
 this->iparm	  = apbsparms.iparm;
 this->ion  	  = apbsparms.ion;
 this->dime 	  = apbsparms.dime;
 this->spacing   = apbsparms.spacing;
 this->centering = apbsparms.centering;
}

void APBSparms::getInputDescription(InputDescription& inputDescription)
{

 // grid geometry options
 // protein
 InputLine("dimension_of_protein_grid", "required", true)
 .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 33))
 .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 33))
 .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 33))
 .addToInputDescription(inputDescription);

 InputLine("spacing_of_protein_grid", "required", true)
 .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
 .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
 .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
 .addToInputDescription(inputDescription);

 InputLine("center_of_protein_grid", "required", true)
 .addArgument(StringArgument(true, ""))
 .addToInputDescription(inputDescription);

 // grid geometry options
 // protein
 InputLine("dimension_of_model_grid", "", true)
 .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 33))
 .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 33))
 .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 33))
 .addToInputDescription(inputDescription);
 
 InputLine("spacing_of_model_grid", "", true)
 .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
 .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
 .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
 .addToInputDescription(inputDescription);

 InputLine("center_of_model_grid", "", true)
 .addArgument(StringArgument(true, ""))
 .addToInputDescription(inputDescription);

 // ion species
 InputLine("ion", "", true)
 .addArgument(FloatArgument(true, 0.0))
 .addArgument(FloatArgument(true, 0.0, 0.0))
 .addArgument(FloatArgument(true, 0.0, 0.0))
 .addToInputDescription(inputDescription);

 // PBE options
 InputLine("bcfl", "bcfl sdh", false)
 .addArgument(StringArgument(true, "sdh"))
 .addToInputDescription(inputDescription);

 InputLine("chgm", "chgm spl0", false)
 .addArgument(StringArgument(true, "spl0"))
 .addToInputDescription(inputDescription);

 InputLine("srfm", "srfm mol", false)
 .addArgument(StringArgument(true, "mol"))
 .addToInputDescription(inputDescription);

 InputLine("pdie", "required", false)
 .addArgument(FloatArgument(true, 4.0, 0.0))
 .addToInputDescription(inputDescription);

 InputLine("sdie", "required", false)
 .addArgument(FloatArgument(true, 80.0, 0.0))
 .addToInputDescription(inputDescription);

 InputLine("srad", "required", false)
 .addArgument(FloatArgument(true, 1.4, 0.0))
 .addToInputDescription(inputDescription);

 InputLine("temperature", "temperature 300.0", false)
 .addArgument(FloatArgument(true, 300.0, 0.0))
 .addToInputDescription(inputDescription);

 InputLine("sdens", "sdens 3.0", false)
 .addArgument(FloatArgument(true, 3.0, 0.0))
 .addToInputDescription(inputDescription);

 // MG options
 InputLine("errtol", "errtol 1.0E-6", false)
 .addArgument(FloatArgument(true, 0.01, 0.0))
 .addToInputDescription(inputDescription);

 InputLine("presmooth", "presmooth 2", false)
 .addArgument(IntegerArgument(true, 2, 0))
 .addToInputDescription(inputDescription);

 InputLine("postsmooth", "postsmooth 2", false)
 .addArgument(IntegerArgument(true, 2, 0))
 .addToInputDescription(inputDescription);

 InputLine("iinfo", "iinfo 0", false)
 .addArgument(IntegerArgument(true, 0, 0, 2))
 .addToInputDescription(inputDescription);

 InputLine("mgcoar", "mgcoar 2", false)
 .addArgument(IntegerArgument(true, 2, 0, 2))
 .addToInputDescription(inputDescription);

 InputLine("mgprol", "mgprol 0", false)
 .addArgument(IntegerArgument(true, 0, 0, 2))
 .addToInputDescription(inputDescription);

 InputLine("mgsmoo", "mgsmoo 1", false)
 .addArgument(IntegerArgument(true, 1, 0, 4))
 .addToInputDescription(inputDescription);

 InputLine("mgsolv", "mgsolv 1", false)
 .addArgument(IntegerArgument(true, 1, 0, 1))
 .addToInputDescription(inputDescription);

 InputLine("mgkey", "mgkey 0", false)
 .addArgument(IntegerArgument(true, 0, 0, 1))
 .addToInputDescription(inputDescription);

 InputLine("itmax", "itmax 200", false)
 .addArgument(IntegerArgument(true, 200, 1))
 .addToInputDescription(inputDescription);

 InputLine("ipcon", "ipcon 3", false)
 .addArgument(IntegerArgument(true, 3, 0, 4))
 .addToInputDescription(inputDescription);

 InputLine("iperf", "iperf 0", false)
 .addArgument(IntegerArgument(true, 0, 0, 3))
 .addToInputDescription(inputDescription);
}

void APBSparms::init(InputParser& parser)
{
 string temp;
 int counter;
 
 /**
  * Input from parser is converted into numerical
  * values and plopped into the parameter arrays.
  */
 
 // nlev
 iparm["nlev"] = 4;

 cout << "Parsing and initializeing APBS parameters ... " << endl;

 /*
  At first we parse the parameters defining the grid geometries.
  There can be different geometries for protein and model. There is
  no sanity check regarding the grid artefact canceling.
 */

 // dimension_of_protein_grid;
 parser.setKeyword("dimension_of_protein_grid");
 counter = 0;
 while( parser.hasNextInputLine() )
 {
  InputLine inputline = parser.getNextInputLine();
  vector<int> temp;
  int nx = (int)inputline.getArgument(1).getLongValue();
  int ny = (int)inputline.getArgument(2).getLongValue();
  int nz = (int)inputline.getArgument(3).getLongValue();
  temp.push_back(nx);
  temp.push_back(ny);
  temp.push_back(nz);
  dime.at(0).push_back( temp );
  cout << "Protein grid points at focusing step " << ++counter << " in x-direction set to " << nx << endl;
  cout << "Protein grid points at focusing step " << counter << " in y-direction set to " << ny << endl;
  cout << "Protein grid points at focusing step " << counter << " in z-direction set to " << nz << endl;
 }
 
 iparm["nfsteps_protein"] = counter;

 // dimension_of_model_grid
 parser.setKeyword("dimension_of_model_grid");
 counter = 0;
 while( parser.hasNextInputLine() )
 {
  InputLine inputline = parser.getNextInputLine();
  vector<int> temp;
  int nx = (int)inputline.getArgument(1).getLongValue();
  int ny = (int)inputline.getArgument(2).getLongValue();
  int nz = (int)inputline.getArgument(3).getLongValue();
  temp.push_back(nx);
  temp.push_back(ny);
  temp.push_back(nz);
  dime.at(1).push_back( temp );
  cout << "Model grid points at focusing step " << ++counter << " in x-direction set to " << nx << endl;
  cout << "Model grid points at focusing step " << counter << " in y-direction set to " << ny << endl;
  cout << "Model grid points at focusing step " << counter << " in z-direction set to " << nz << endl;
 }
 iparm["nfsteps_model"] = counter;

 // if no model grid specified
 if (counter==0)
 {
  dime.at(1) = dime.at(0);
  iparm["nfsteps_model"] = iparm["nfsteps_protein"];
  cout << "You did not specify the grid dimension in model calculations. Using same values as in protein calculations." << endl;
 }

 // spacing_of_protein_grid
 parser.setKeyword("spacing_of_protein_grid");
 counter = 0;
 while( parser.hasNextInputLine() )
 {
  InputLine inputline = parser.getNextInputLine();
  vector<double> temp;
  double hx = inputline.getArgument(1).getDoubleValue();
  double hy = inputline.getArgument(2).getDoubleValue();
  double hz = inputline.getArgument(3).getDoubleValue();
  temp.push_back(hx);
  temp.push_back(hy);
  temp.push_back(hz);
  spacing.at(0).push_back( temp );
  cout << "Protein grid resolution in step " << ++counter << " in x-direction set to " << hx << endl;
  cout << "Protein grid resolution in step " << counter << " in y-direction set to " << hy << endl;
  cout << "Protein grid resolution in step " << counter << " in z-direction set to " << hz << endl;
 }

 // sanitize
 if( spacing.at(0).size()!=iparm["nfsteps_protein"] )
 {
  string message = "The number of 'spacing_of_protein_grid' lines must balance the number of 'dimension_of_protein_grid' lines!";
  throw ParameterInconsistencyException(message);
 }
 
 // spacing_of_model_grid
 parser.setKeyword("spacing_of_model_grid");
 counter = 0;
 while( parser.hasNextInputLine() )
 {
  InputLine inputline = parser.getNextInputLine();
  vector<double> temp;
  double hx = inputline.getArgument(1).getDoubleValue();
  double hy = inputline.getArgument(2).getDoubleValue();
  double hz = inputline.getArgument(3).getDoubleValue();
  temp.push_back(hx);
  temp.push_back(hy);
  temp.push_back(hz);
  spacing.at(1).push_back( temp );
  cout << "Model grid resolution in step " << ++counter << " in x-direction set to " << hx << endl;
  cout << "Model grid resolution in step " << counter << " in y-direction set to " << hy << endl;
  cout << "Model grid resolution in step " << counter << " in z-direction set to " << hz << endl;
 }

 // if no model grid spacing specified
 if (counter==0)
 {
  spacing.at(1) = spacing.at(0);
  cout << "You did not specify the grid spacing in model calculations. Using same values as in protein calculations." << endl;
 }
 
 // sanitize
 if( spacing.at(1).size()!=iparm["nfsteps_model"] )
 {
  string message = "The number of 'spacing_of_model_grid' lines must balance the number of 'dimension_of_model_grid' lines!";
  throw ParameterInconsistencyException(message);
 }

 // centering
 parser.setKeyword("center_of_protein_grid");
 counter = 0;
 while( parser.hasNextInputLine() )
 {
  string center = parser.getNextInputLine().getArgument(1).getStringValue();
  if ( center==ON_PROTEIN || center==ON_MODEL )
  {
   centering.at(0).push_back( center );
   cout << "Centering method in step " << ++counter << " in protein calculations: " << center << endl;
  }
  else
  {
   ostringstream ostr;
   ostr << "Unknown centering method '" << center << "'";
	throw ParameterInconsistencyException(ostr.str());
  }
 }
 if( centering.at(0).size()!=iparm["nfsteps_protein"] )
 {
  string message = "The number of 'center_of_protein_grid' lines must balance the number of 'dimension_of_protein_grid' lines!";
  throw ParameterInconsistencyException(message);
 }

 parser.setKeyword("center_of_model_grid");
 counter = 0;
 while( parser.hasNextInputLine() )
 {
  string center = parser.getNextInputLine().getArgument(1).getStringValue();
  if ( center==ON_PROTEIN || center==ON_MODEL )
  {
   centering.at(1).push_back( center );
   cout << "Centering method in step " << ++counter << " in model calculations: " << center << endl;
  }
  else
  {
   ostringstream ostr;
   ostr << "Unknown centering method '" << center << "'";
	throw ParameterInconsistencyException(ostr.str());
  }
 }

 if (counter==0)
 {
  centering.at(1) = centering.at(0);
  cout << "You did not specify the grid centering method in model calculations. Using same values as in protein calculations." << endl;
 }

 if( centering.at(1).size()!=iparm["nfsteps_model"] )
 {
  string message = "The number of 'center_of_model_grid' lines must balance the number of 'dimension_of_model_grid' lines!";
  throw ParameterInconsistencyException(message);
 }

 /** 
  Now the easy part: PBE parameters + surface definition.
 */

 // ion + nions
 parser.setKeyword("ion");
 counter = 0;
 while( parser.hasNextInputLine() )
 {
  InputLine inputline = parser.getNextInputLine();
  vector<double> temp;
  double ionq = inputline.getArgument(1).getDoubleValue();
  double ionc = inputline.getArgument(2).getDoubleValue();
  double ionr = inputline.getArgument(3).getDoubleValue();
  temp.push_back(ionq);
  temp.push_back(ionc);
  temp.push_back(ionr);
  ion.push_back(temp);
  cout << "Ion species " << ++counter << ": "
       << "charge = " << ionq << ", "
       << "conc = " << ionc << " M, "
       << "radius = " << ionr << " A." 
       << endl;
 }
 iparm["nions"] = counter;
 
 // pdie
 rparm["pdie"] = parser.getInputLine("pdie").getArgument(1).getDoubleValue();
 cout << "Protein dielectric constant set to " << rparm["pdie"] << endl;

 // sdie
 rparm["sdie"] = parser.getInputLine("sdie").getArgument(1).getDoubleValue();
 cout << "Solvent dielectric constant set to " << rparm["sdie"] << endl;
 
 // temperature
 rparm["temperature"] = parser.getInputLine("temperature").getArgument(1).getDoubleValue();
 cout << "Temperature set to " << rparm["temperature"] << endl;
 
 // srad
 rparm["srad"] = parser.getInputLine("srad").getArgument(1).getDoubleValue();
 cout << "Solvent probe radius set to " << rparm["srad"] << endl;
 
#ifdef APBS4
 // sdens
 rparm["sdens"] = parser.getInputLine("sdens").getArgument(1).getDoubleValue();
 cout << "sdens set to " << rparm["sdens"] << endl;
#endif

 //bcfl
 temp = parser.getInputLine("bcfl").getArgument(1).getStringValue();
 if (temp == "zero")
  iparm["bcfl"] = 0;
 else if (temp == "mdh")
  iparm["bcfl"] = 2;
 else
  iparm["bcfl"] = 1;
 cout << "Boundary condition set to " << iparm["bcfl"] << endl;

 //chgm
 temp = parser.getInputLine("chgm").getArgument(1).getStringValue();
 if (temp == "spl2")
  iparm["chgm"] = 1;
 else
  iparm["chgm"] = 0;
 cout << "Charge discretization method set to " << iparm["chgm"] << endl;

 //srfm
 temp = parser.getInputLine("srfm").getArgument(1).getStringValue();
 if (temp == "smol")
  iparm["srfm"] = 1;
 else
  iparm["srfm"] = 0;
 cout << "SASA model set to " << iparm["srfm"] << endl;

 // errtol
 rparm["errtol"] = parser.getInputLine("errtol").getArgument(1).getDoubleValue();
 cout << "Error tolerance set to " << rparm["errtol"] << endl;

 // presmooth
 iparm["presmooth"] = (int)parser.getInputLine("presmooth").getArgument(1).getLongValue();
 cout << "presmooth set to " << iparm["presmooth"] << endl;

 // postsmooth
 iparm["postsmooth"] = (int)parser.getInputLine("postsmooth").getArgument(1).getLongValue();
 cout << "postsmooth set to " << iparm["postsmooth"] << endl;

 // iinfo
 iparm["iinfo"] = (int)parser.getInputLine("iinfo").getArgument(1).getLongValue();
 cout << "iinfo set to " << iparm["iinfo"] << endl;

 /* experimental performance tweaks */

 // ipcon
 iparm["ipcon"] = (int)parser.getInputLine("ipcon").getArgument(1).getLongValue();
 cout << "ipcon set to " << iparm["ipcon"] << endl;

 // iperf
 iparm["iperf"] = (int)parser.getInputLine("iperf").getArgument(1).getLongValue();
 cout << "iperf set to " << iparm["iperf"] << endl;

 // mgcoar
 iparm["mgcoar"] = (int)parser.getInputLine("mgcoar").getArgument(1).getLongValue();
 cout << "mgcoar set to " << iparm["mgcoar"] << endl;

 // mgkey
 iparm["mgkey"] = (int)parser.getInputLine("mgkey").getArgument(1).getLongValue();
 cout << "mgkey set to " << iparm["mgkey"] << endl;

 // mgprol
 iparm["mgprol"] = (int)parser.getInputLine("mgprol").getArgument(1).getLongValue();
 cout << "mgprol set to " << iparm["mgprol"] << endl;

 // mgsmoo
 iparm["mgsmoo"] = (int)parser.getInputLine("mgsmoo").getArgument(1).getLongValue();
 cout << "mgsmoo set to " << iparm["mgsmoo"] << endl;

 // mgsolv
 iparm["mgsolv"] = (int)parser.getInputLine("mgsolv").getArgument(1).getLongValue();
 cout << "mgsolv set to " << iparm["mgsolv"] << endl;

 // mgsolv
 iparm["itmax"] = (int)parser.getInputLine("itmax").getArgument(1).getLongValue();
 cout << "itmax set to " << iparm["itmax"] << endl;

 cout << "... finished." << endl;
}

void APBSparms::get_ion(double ionq[], double ionc[], double ionr[])
{
 try
 {
  for (int i=0; i<get_nions(); i++)
  {
   ionq[i] = ion.at(i).at(0);
   ionc[i] = ion.at(i).at(1);
   ionr[i] = ion.at(i).at(2);
  }
 }
 catch (out_of_range e)
 {
  throw OutOfRangeException("APBSparms::get_ion");
 }
}

void APBSparms::get_dime(int fstep, int array[])
{
 try
 {
  array[0] = dime.at(0).at(fstep).at(0);
  array[1] = dime.at(0).at(fstep).at(1);
  array[2] = dime.at(0).at(fstep).at(2);
 }
 catch(out_of_range e)
 {
  throw OutOfRangeException("APBSparms::get_dime");
 }
}

void APBSparms_model::get_dime(int fstep, int array[])
{
 try
 {
  array[0] = dime.at(1).at(fstep).at(0);
  array[1] = dime.at(1).at(fstep).at(1);
  array[2] = dime.at(1).at(fstep).at(2);
 }
 catch(out_of_range e)
 {
  throw OutOfRangeException("APBSparms_model::get_dime");
 }
}

void APBSparms::get_grid(int fstep, double array[])
{
 try
 {
  array[0] = spacing.at(0).at(fstep).at(0);
  array[1] = spacing.at(0).at(fstep).at(1);
  array[2] = spacing.at(0).at(fstep).at(2);
 }
 catch(out_of_range e)
 {
  throw OutOfRangeException("APBSparms::get_grid");
 }
}

void APBSparms_model::get_grid(int fstep, double array[])
{
 try
 {
  array[0] = spacing.at(1).at(fstep).at(0);
  array[1] = spacing.at(1).at(fstep).at(1);
  array[2] = spacing.at(1).at(fstep).at(2);
 }
 catch(out_of_range e)
 {
  throw OutOfRangeException("APBSparms_model::get_grid");
 }
}

string APBSparms::get_centering(int fstep)
{
 string centering_method;

 try
 {
  centering_method = centering.at(0).at(fstep);
 }
 catch(out_of_range e)
 {
  throw OutOfRangeException("APBSparms::get_centering");
 }
 
 return centering_method;
}

string APBSparms_model::get_centering(int fstep)
{
 string centering_method;

 try
 {
  centering_method = centering.at(1).at(fstep);
 }
 catch(out_of_range e)
 {
  throw OutOfRangeException("APBSparms_model::get_centering");
 }
 
 return centering_method;
}

/*---------------------------------- APBS_base ---------------------------------- */

APBS_base::APBS_base(int debuglevel, Structure *pqr, APBSparms *apbsparms)
{
 int rank, size;
 
 this->debug     = debuglevel;
 this->apbsparms = apbsparms;
 ncalc           = apbsparms->get_nfsteps();
 
 /* arrays */
 pbe      = new Vpbe*[ncalc];
 pmg      = new Vpmg*[ncalc];
 pmgp     = new Vpmgp*[ncalc];
 chargeMap= new Vgrid*[ncalc];
 dielX    = new Vgrid*[ncalc];
 dielY    = new Vgrid*[ncalc];
 dielZ    = new Vgrid*[ncalc];
 kappa    = new Vgrid*[ncalc];
 dummyMap = new Vgrid*[1];

 /* initialize */
 nosh     = VNULL;
 mem      = VNULL;
 com      = VNULL;

 /* APBS I/O and memory manager */
 com = Vcom_ctor(1);
 rank = Vcom_rank(com);
 size = Vcom_size(com);
 startVio(); /* FIXME */
 Vnm_setIoTag(rank, size);
 Vnm_tprint( 0, "apbs:  Hello world from PE %d\n", rank);
 mem = Vmem_ctor("MAIN");
 
 for (int i=0; i<ncalc; i++)
 {
  pmg[i]      = VNULL;
  pmgp[i]     = VNULL;
  pbe[i]      = VNULL;
  chargeMap[i]= VNULL;
  dielX[i]    = VNULL;
  dielY[i]    = VNULL;
  dielZ[i]    = VNULL;
  kappa[i]    = VNULL;
 }

 /** 
  NOsh parameters setup -
  the NOsh APBS class is a parser class
  used to parse and store user parameters.
 */

 nosh = NOsh_ctor(rank,size);
 nosh->proc_rank = rank;
 nosh->proc_size = size;
 nosh->ispara = 0;
 nosh->nmol   = 1;  /* we only have a single molecule */
 nosh->ncalc  = 0;  /* how many focusing steps? */
 nosh->nelec  = 1;  /* how many ELEC sections? */
 nosh->nprint = 0;  /* how many PRINT sections? */
 nosh->parsed = 1;  /* finished parsing? */

 /* PQR structure */
 this->pqr = pqr;
 
 // APBS structures have still to be intialized
 setupflag = 0;
}


MGparm* APBS_base::fillMGparm(MGparm_CalcType type)
{
 /** 
  Create and fill MGparm object that stores
  multi grid user parameters.
 */
 
 MGparm *cmgparm = MGparm_ctor(type); /* 0 -> MGMANUAL, 3 -> MGDUMMY */
 
 /* parameters */
 cmgparm->nlev    = apbsparms->get_nlev();
 cmgparm->cmeth   = MCM_POINT;
 cmgparm->chgm    = apbsparms->get_chgm();
 apbsparms->get_dime(0, cmgparm->dime);
 apbsparms->get_grid(0, cmgparm->grid);

 cmgparm->centmol = 1; /* molecule ID we are 
                          centering on 
								  if we are centering 
								  on molecule - that will never happen
							  */

 /* fill flags */

 cmgparm->setdime   = 1;
 cmgparm->setgcent  = 1;
 cmgparm->setcgcent = 1;
 cmgparm->setfgcent = 1;
 cmgparm->setcglen  = 0;
 cmgparm->setfglen  = 0;
 cmgparm->setglen   = 0;
 cmgparm->setgrid   = 1;
 cmgparm->setchgm   = 1;

 /* useless in mgmanual/dummy mode */
 cmgparm->ccentmol  = 1;
 cmgparm->fcentmol  = 1;
 cmgparm->ccmeth    = MCM_POINT;
 cmgparm->fcmeth    = MCM_POINT;
 cmgparm->pdime[0]  = 0;
 cmgparm->pdime[1]  = 0;
 cmgparm->pdime[2]  = 0;
 cmgparm->ofrac     = 0.1;
 /* ------------------------ */

 cmgparm->parsed = 1;
 
 return cmgparm;
}

PBEparm* APBS_base::fillPBEparm()
{
 /**
  Creates and fills PBEparm object with all necessary PBE parameters
 */

 PBEparm *cpbeparm = PBEparm_ctor();

 /* non-user setup */
 cpbeparm->useDielMap  = 0;
 cpbeparm->dielMapID   = 0;
 cpbeparm->useKappaMap = 0;
 cpbeparm->kappaMapID  = 0;
 cpbeparm->useChargeMap= 0;
 cpbeparm->usePotMap   = 0;
 cpbeparm->chargeMapID = 0;
 cpbeparm->writemat    = 0;
 cpbeparm->numwrite    = 0;
 cpbeparm->molid       = 1;
 cpbeparm->pbetype     = PBE_LPBE;
 cpbeparm->calcenergy  = PCE_NO;
 cpbeparm->calcforce   = PCF_NO;
 //cpbeparm->gamma       = 0.105;
 cpbeparm->swin        = 0.3;

 /* user setup */
 cpbeparm->bcfl        = apbsparms->get_bcfl(); /* will be changed to BCFL_FOCUS later */
 cpbeparm->nion        = apbsparms->get_nions();
 cpbeparm->srfm        = apbsparms->get_srfm();
 cpbeparm->pdie        = apbsparms->get_pdie();
 cpbeparm->sdie        = apbsparms->get_sdie();
 cpbeparm->srad        = apbsparms->get_srad();
 cpbeparm->temp        = apbsparms->get_temperature();
#ifdef APBS4
 cpbeparm->sdens       = apbsparms->get_sdens();
#endif

 /* copy ion species */
 apbsparms->get_ion(cpbeparm->ionq, cpbeparm->ionc, cpbeparm->ionr);
 for(int i=0; i<apbsparms->get_nions(); ++i)
  cpbeparm->setion[i] = 1;

 /* set flags */
 cpbeparm->setpbetype    = 1;
 cpbeparm->setbcfl       = 1;
 cpbeparm->setnion       = 1;
 cpbeparm->setpdie       = 1;
 cpbeparm->setsdie       = 1;
 cpbeparm->setsrfm       = 1;
 cpbeparm->setsrad       = 1;
 cpbeparm->setswin       = 1;
 cpbeparm->settemp       = 1;
// cpbeparm->setgamma      = 1;
 cpbeparm->setcalcenergy = 0;
 cpbeparm->setcalcforce  = 0;
#ifdef APBS4
 cpbeparm->setsdens      = 1;
#endif
 cpbeparm->setwritemat   = 1;
 cpbeparm->parsed        = 1;
 cpbeparm->setmolid      = 1;
 
 return cpbeparm;
}

Valist* APBS_base::pqr2alist(Structure *pqr)
{
 int i = 0;
 double pos[3];

 // Create a new Valist-object.
 Valist* alist = Valist_ctor();
         alist->center[0] = 0.;
         alist->center[1] = 0.;
         alist->center[2] = 0.;
         alist->maxcrd[0] = -VLARGE;
         alist->maxcrd[1] = -VLARGE;
         alist->maxcrd[2] = -VLARGE;
         alist->mincrd[0] = VLARGE;
         alist->mincrd[1] = VLARGE;
         alist->mincrd[2] = VLARGE;
         alist->maxrad    = 0.;
         alist->charge    = 0.;
         alist->number    = pqr->nrOfAtoms();
         alist->atoms = (Vatom*)Vmem_malloc(alist->vmem, alist->number, ( sizeof(Vatom) ) );

 // sanitize
 VASSERT(alist->atoms != VNULL);

 // copy data
 pqr->get_center(&(alist->center[0]), &(alist->center[1]), &(alist->center[2]));
 pqr->get_minmax(alist->mincrd, alist->maxcrd);

 while( pqr->hasNext() )
 {
  Satom& atom    = pqr->getNext();
  double radius = atom.get_radius();
  if (radius > alist->maxrad) alist->maxrad = radius;
  pos[0] = atom.coordinates().get_x_value();
  pos[1] = atom.coordinates().get_y_value();
  pos[2] = atom.coordinates().get_z_value();
  Vatom_setPosition(&(alist->atoms)[i], pos);
  Vatom_setCharge  (&(alist->atoms)[i], atom.get_charge());
  Vatom_setRadius  (&(alist->atoms)[i], atom.get_radius());
#ifdef APBS4
  Vatom_setAtomID  (&(alist->atoms)[i], i);
#endif
  i++;
 }

 /* debug
 if (debug > 2) 
 {
  printf("alist->centerx = %f\n", alist->center[0]);
  printf("alist->centery = %f\n", alist->center[1]);
  printf("alist->centerz = %f\n", alist->center[2]);
  for (i=0; i<alist->number; i++)
  {
   printf("apbs: alist->atom: %f %f %f %f %f\n", 
   alist->atoms[i].position[0],
	alist->atoms[i].position[1],
	alist->atoms[i].position[2],
	alist->atoms[i].charge,
	alist->atoms[i].radius);
  }
 }
 */

 return alist;
}

void APBS_base::setupManualFocus(MGparm_CalcType type, double centerx[], double centery[], double centerz[], MapController *mapcontroller)
{
 cout << "Entering APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[], Mapcontroller*)." << endl;

 if (setupflag)
  throw APBSWrapperSetupCalledTwiceException();

 PBEparm **cpbeparm = new PBEparm*[ncalc];
 MGparm  **cmgparm  = new MGparm*[ncalc];

 /* Convert PQR structure into APBS atom list */
 alist = pqr2alist(pqr);

 /* fill PBEparm object */
 cpbeparm[0] = fillPBEparm();

 /* fill MGparm object */
 cmgparm[0]  = fillMGparm(type);

 /* fill external dielectric maps */
 mapcontroller->fillMaps(*pqr, centerx, centery, centerz, this->dielX, this->dielY, this->dielZ, this->kappa);
 
 /** 
  Call NOsh_setup to set up all remaining 
  APBS data structures. Initialize 
  focusing steps.
 */

 for (int i=0; i<ncalc; i++)
 {
  if (i>0)
  {
   /* copy values from intial focusing step */
   cpbeparm[i] = PBEparm_ctor();
	PBEparm_copy(cpbeparm[i], cpbeparm[0]);
	cmgparm[i] = MGparm_ctor(type);
	MGparm_copy(cmgparm[i], cmgparm[0]);
   
   /* Overwrite some existing values */
   cpbeparm[i]->bcfl = BCFL_FOCUS; /* focusing boundary condition */
  }

  /* use external accessibility maps */
  cpbeparm[i]->useDielMap  = 1;
  cpbeparm[i]->useKappaMap = 1;
  cpbeparm[i]->usePotMap   = 0;
  nosh->ndiel++;
  nosh->nkappa++;
  cpbeparm[i]->dielMapID     = i+1;
  cpbeparm[i]->kappaMapID    = i+1;

  /* set grid resolution and center */
  apbsparms->get_dime(i, cmgparm[i]->dime);
  apbsparms->get_grid(i, cmgparm[i]->grid);
  cmgparm[i]->center[0] = centerx[i];
  cmgparm[i]->center[1] = centery[i];
  cmgparm[i]->center[2] = centerz[i];
  
  /* debug */
  if (debug > 1) 
  {
   cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->useDielMap = " << cpbeparm[i]->useDielMap << endl;
   cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->dielMapID = " << cpbeparm[i]->dielMapID << endl;
	cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->useKappaMap = " << cpbeparm[i]->useKappaMap << endl;
	cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->kappaMapID = " << cpbeparm[i]->kappaMapID << endl;
#ifdef APBS4
	cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->sdens = " << cpbeparm[i]->sdens << endl;
#endif
  }
 
  /* check PBE parameters loaded so far */
  if (!PBEparm_check(cpbeparm[i]))
   throw APBSErrorException("PBEparm_check", "APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[], Mapcontroller*)", 0);

  /* check MG parameters loaded so far */
  /* MGparm_check will also reset nlev appropriately! */
  if (!MGparm_check(cmgparm[i]))
   throw APBSErrorException("MGparm_check", "APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[], Mapcontroller*)", 0);

  if (cmgparm[i]->setgrid == 0) {
	  cmgparm[i]->grid[0] = cmgparm[i]->glen[0]/((double)(cmgparm[i]->dime[0]-1));
	  cmgparm[i]->grid[1] = cmgparm[i]->glen[1]/((double)(cmgparm[i]->dime[1]-1));
	  cmgparm[i]->grid[2] = cmgparm[i]->glen[2]/((double)(cmgparm[i]->dime[2]-1));
  }
  if (cmgparm[i]->setglen == 0) {
	  cmgparm[i]->glen[0] = cmgparm[i]->grid[0]*((double)(cmgparm[i]->dime[0]-1));
	  cmgparm[i]->glen[1] = cmgparm[i]->grid[1]*((double)(cmgparm[i]->dime[1]-1));
	  cmgparm[i]->glen[2] = cmgparm[i]->grid[2]*((double)(cmgparm[i]->dime[2]-1));
  }

  //Get the next calculation object and increment the number of calculations
  nosh->calc[nosh->ncalc] = NOsh_calc_ctor(NCT_MG);
  (nosh->ncalc)++;

    PBEparm_copy(nosh->calc[i]->pbeparm, cpbeparm[i]);
  	MGparm_copy(nosh->calc[i]->mgparm, cmgparm[i]);
 }
 
 delete[] cpbeparm;
 delete[] cmgparm;
 
 setupflag = 2;
 
 cout << "Exiting APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[], Mapcontroller*)." << endl;
}

void APBS_base::setupManualFocus(MGparm_CalcType type, double centerx[], double centery[], double centerz[])
{
 cout << "Entering APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[])." << endl;

 if (setupflag)
  throw APBSWrapperSetupCalledTwiceException();
 
 PBEparm **cpbeparm = new PBEparm*[ncalc];
 MGparm  **cmgparm  = new MGparm*[ncalc];

 /* Convert PQR structure into APBS atom list */
 alist = pqr2alist(pqr);

 /* fill PBEparm object */
 cpbeparm[0] = fillPBEparm();

 /* fill MGparm object */
 cmgparm[0]  = fillMGparm(type);

 /** 
  Call NOsh_setup to set up all remaining 
  APBS data structures. Initialize 
  focusing steps.
 */
 
 for (int i=0; i<ncalc; i++)
 {
  if (i>0)
  {
   /* copy values from intial focusing step */
	
   cpbeparm[i] = PBEparm_ctor();
	PBEparm_copy(cpbeparm[i], cpbeparm[0]);
	cmgparm[i] = MGparm_ctor(type);
	MGparm_copy(cmgparm[i], cmgparm[0]);
   
   /* Overwrite some existing values */

   cpbeparm[i]->bcfl = BCFL_FOCUS; /* focusing boundary condition */
  }

  /* set grid length and center */
  
  apbsparms->get_dime(i, cmgparm[i]->dime);
  
  apbsparms->get_grid(i, cmgparm[i]->grid);
  cmgparm[i]->center[0] = centerx[i];
  cmgparm[i]->center[1] = centery[i];
  cmgparm[i]->center[2] = centerz[i];
  
  /* debug */
  
  if (debug > 1) 
  {
   cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->useDielMap = " << cpbeparm[i]->useDielMap << endl;
   cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->dielMapID = " << cpbeparm[i]->dielMapID << endl;
	cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->useKappaMap = " << cpbeparm[i]->useKappaMap << endl;
	cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->kappaMapID = " << cpbeparm[i]->kappaMapID << endl;
#ifdef APBS4
	cout << "APBS_base::setupManualFocus: cpbeparm[" << i << "]->sdens = " << cpbeparm[i]->sdens << endl;
#endif
  }

  /* check PBE parameters loaded so far */
  if (!PBEparm_check(cpbeparm[i]))
   throw APBSErrorException("PBEparm_check", "APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[])", 0); 

  /* check MG parameters loaded so far */
  /* MGparm_check will also reset nlev appropriately! */
  if (!MGparm_check(cmgparm[i])) 
   throw APBSErrorException("MGparm_check", "APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[])", 0);

//  if ( !NOsh_setupMGMANUAL(nosh, cmgparm[i], cpbeparm[i]) )
//     throw APBSErrorException("NOsh_setupMGMANUAL", "APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[])", 0);
 }
 
 delete[] cpbeparm; // PBEparm objects will be destroyed in NOsh_dtor
 delete[] cmgparm;  // MGparm objects will be destroyed in NOsh_dtor
 
 setupflag = 1;

 cout << "Exiting APBS_base::setupManualFocus(MGparm_CalcType double[], double[], double[])." << endl;
}

Vgrid* APBS_base::exportMap(Vpmg *pmg, PBEparm *pbeparm, Vdata_Type maptype)
{
 /**
  This subroutine wraps the 
  data array in pmg specified by maptype
  into a new Vgrid object. It does not
  copy the contents of the data array!
 */
 
 double *data;
  
 int nx = pmg->pmgp->nx;
 int ny = pmg->pmgp->ny;
 int nz = pmg->pmgp->nz;
 double hx = pmg->pmgp->hx;
 double hy = pmg->pmgp->hy;
 double hzed = pmg->pmgp->hzed;
  
 double xmin, ymin, zmin;

 xmin = -0.5*(nx-1)*hx;
 ymin = -0.5*(ny-1)*hy;
 zmin = -0.5*(nz-1)*hzed;
  
 switch(maptype)
 {
  case VDT_DIELX: 
   xmin  += pmg->pmgp->xcent + 0.5*hx;
   ymin  += pmg->pmgp->ycent;
   zmin  += pmg->pmgp->zcent;
	data   = pmg->epsx;
   break;
  case VDT_DIELY: 
   xmin  += pmg->pmgp->xcent; 
   ymin  += pmg->pmgp->ycent + 0.5*hy;
   zmin  += pmg->pmgp->zcent;
   data   = pmg->epsy;
	break;
  case VDT_DIELZ:
   xmin  += pmg->pmgp->xcent; 
   ymin  += pmg->pmgp->ycent;
   zmin  += pmg->pmgp->zcent + 0.5*hzed;
   data   = pmg->epsz;
   break;
  case VDT_KAPPA:
   xmin  += pmg->pmgp->xcent; 
   ymin  += pmg->pmgp->ycent;
   zmin  += pmg->pmgp->zcent;
	data   = pmg->kappa;
   break;
  case VDT_POT:
   xmin  += pmg->pmgp->xcent; 
   ymin  += pmg->pmgp->ycent;
   zmin  += pmg->pmgp->zcent;
	data   = pmg->u;
   break;
  default:
   xmin  += pmg->pmgp->xcent; 
   ymin  += pmg->pmgp->ycent;
   zmin  += pmg->pmgp->zcent;
	data   = VNULL;
	break;
 }

 /* If grid data is ctored, it will not be freed by Vgrid_dtor() */
 
 return Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin, data);
}

APBS_base::~APBS_base()
{
 // accessibility maps
 killDielMaps(nosh, dielX, dielY, dielZ);
 killKappaMaps(nosh, kappa);

 // mg subsystem
 if (setupflag == 1)
 {
  killMG(nosh, pbe, pmgp, pmg);
 }
 else
 {
  Vpbe_dtor(&(pbe[nosh->ncalc-1]));
  Vpmg_run_dtor(&(pmg[nosh->ncalc-1]));
  Vpmgp_dtor(&(pmgp[nosh->ncalc-1]));
 }

 // alists
 killMolecules(nosh, &alist); 

 // nosh parser object
 NOsh_dtor(&nosh);
 
 /* Memory statistics */
 int bytesTotal = Vmem_bytesTotal();
 int highWater = Vmem_highWaterTotal();
 
 if ( debug > 1 )
  cout << "APBS:: Final memory usage: " 
       << (double)(bytesTotal)/(1024.*1024.) << " MB, "
       << (double)(highWater)/(1024.*1024.) << " MB high water." << endl;

 /* Clean up MALOC structures */
 Vcom_dtor(&com);
 Vmem_dtor(&mem);
 Vnm_tstop(26, "APBS WALL CLOCK");
 
 /* kill arrays */
 delete[] pbe;
 delete[] pmg;
 delete[] pmgp;
 delete[] dielX;
 delete[] dielY;
 delete[] dielZ;
 delete[] kappa;
 delete[] chargeMap;
 cout << "Good bye from APBS!" << endl;
}

/* ---------------------------------- APBS_runner ---------------------------------- */

APBS_runner::APBS_runner(int debuglevel, Structure *pqr, APBSparms *apbsparms) : APBS_base(debuglevel, pqr, apbsparms)
{
 cout << "Hello World from APBS_runner!" << endl;
}

int APBS_runner::initAPBScore(int i, MGparm *mgparm, PBEparm *pbeparm, double realCenter[])
{
 /**
  This method replaces initMG for APBS runs using pre-calculated map.
  It is very close to the APBS routine, however, it calls alternative
  constructors for APBS core structures.
 */
 int focusFlag;
 double sparm, iparm;

 /* Set up PBE object */
 Vnm_tprint(0, "Setting up PBE object...\n");
 if (pbeparm->srfm == VSM_SPLINE) 
  sparm = pbeparm->swin;
 else 
  sparm = pbeparm->srad;

 if (pbeparm->nion > 0) 
  iparm = pbeparm->ionr[0];
 else 
  iparm = 0.0;

 // is this a focusing calculation?
 if (pbeparm->bcfl == BCFL_FOCUS) 
 {
  if (i == 0) 
  {
   cerr << "Can't focus first calculation!" << endl;
   return 0;
  }
  focusFlag = 1;
 }
 else 
  focusFlag = 0;

 pbe[i] = Vpbe_run_ctor(alist, pbeparm->nion,
                    pbeparm->ionc, pbeparm->ionr,
                    pbeparm->ionq, pbeparm->temp,
                    pbeparm->pdie, pbeparm->sdie,
                    sparm, focusFlag,
                    pbeparm->sdens, 0, 0, 0, 0);
 /* Set up PDE object */
 Vnm_tprint(0, "Setting up PDE object...\n");
 if (pbeparm->pbetype == PBE_LPBE)
 {

  pmgp[i] = Vpmgp_ctor(mgparm);

  // performance tweaking
  pmgp[i]->errtol = apbsparms->get_errtol();
  pmgp[i]->nu1    = apbsparms->get_presmooth();
  pmgp[i]->nu2    = apbsparms->get_postsmooth();
  pmgp[i]->iinfo  = apbsparms->get_iinfo();
  /* experimental */
  pmgp[i]->mgcoar = apbsparms->get_mgcoar();
  pmgp[i]->mgprol = apbsparms->get_mgprol();
  pmgp[i]->mgsmoo = apbsparms->get_mgsmoo();
  pmgp[i]->mgsolv = apbsparms->get_mgsolv();
  pmgp[i]->mgkey  = apbsparms->get_mgkey();
  pmgp[i]->itmax  = apbsparms->get_itmax();
  pmgp[i]->ipcon  = apbsparms->get_ipcon();
  pmgp[i]->iperf  = apbsparms->get_iperf();
  pmgp[i]->ipkey  = IPKEY_LPBE;
  pmgp[i]->meth   = 2;
  /* experimental */
  if (debug>1)
  {
   cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->errtol = " << pmgp[i]->errtol << endl;
	cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->itmax = " << pmgp[i]->itmax << endl;
   cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->nu1 = " << pmgp[i]->nu1 << endl;
   cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->nu1 = " << pmgp[i]->nu2 << endl;
   cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->iinfo = " << pmgp[i]->iinfo << endl;
	cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->mgcoar = " << pmgp[i]->mgcoar << endl;
	cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->mgprol = " << pmgp[i]->mgprol << endl;
	cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->mgsmoo = " << pmgp[i]->mgsmoo << endl;
	cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->mgsolv = " << pmgp[i]->mgsolv << endl;
	cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->mgkey = " << pmgp[i]->mgkey << endl;
	cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->ipcon = " << pmgp[i]->ipcon << endl;
   cout << "APBS_runner::initAPBScore: pmgp["<<i<<"]->iperf = " << pmgp[i]->iperf << endl;
  }
 }
 else
 {
  Vnm_tprint( 2, "Non-linear PBE solver selected!\n");
  return 0;  
 }
 
 Vnm_tprint(0, "Setting PDE center to local center...\n");
 pmgp[i]->bcfl  = pbeparm->bcfl;
 // since we do not support MGPARA, no center shifting is needed
 realCenter[0]  = mgparm->center[0];
 realCenter[1]  = mgparm->center[1];
 realCenter[2]  = mgparm->center[2];
 pmgp[i]->xcent = realCenter[0];
 pmgp[i]->ycent = realCenter[1];
 pmgp[i]->zcent = realCenter[2];
 
 if (debug>1)
  cout << "APBS_runner::initAPBScore: Grid x/y/z-center: " 
       << pmgp[i]->xcent << "/" 
		 << pmgp[i]->ycent << "/" 
		 << pmgp[i]->zcent << endl;
 
 if (focusFlag)
  pmg[i] = Vpmg_run_ctor(pmgp[i], pbe[i], focusFlag, pmg[i-1], mgparm, pbeparm->calcenergy);
 else
  pmg[i] = Vpmg_run_ctor(pmgp[i], pbe[i], focusFlag, VNULL, mgparm, pbeparm->calcenergy);
	 
 if (i>0) 
 {
  Vpmgp_dtor(&(pmgp[i-1]));
  Vpbe_dtor(&(pbe[i-1]));
 }

 Vpmg_run_fillco(pmg[i], pbeparm->srfm, pbeparm->swin, mgparm->chgm, 
                 pbeparm->useDielMap, dielX[i],
                 pbeparm->useDielMap, dielY[i],
                 pbeparm->useDielMap, dielZ[i],
                 pbeparm->useKappaMap, kappa[i],
                 pbeparm->useChargeMap, VNULL);
 
 return 1;
}

void APBS_runner::setup(double centerx[], double centery[], double centerz[])
{
 setupManualFocus(MCT_MANUAL, centerx, centery, centerz);
}

void APBS_runner::setup(double centerx[], double centery[], double centerz[], MapController *mapcontroller)
{
 setupManualFocus(MCT_MANUAL, centerx, centery, centerz, mapcontroller);
}

void APBS_runner::setupDummy(double centerx[], double centery[], double centerz[])
{
 setupManualFocus(MCT_DUMMY, centerx, centery, centerz);
}

void APBS_runner::setupDummy(double centerx[], double centery[], double centerz[], MapController *mapcontroller)
{
 setupManualFocus(MCT_DUMMY, centerx, centery, centerz, mapcontroller);
}

void APBS_runner::run()
{
 /**
  starts the APBS calculation
  with current parameter setup
 */

 cout << "Entering APBS_runner::run." << endl;
 if (!setupflag)
  throw APBSWrapperSetupNotCalledException("APBS_runner::run");

 double realCenter[3]; 
 MGparm  *mgparm;
 PBEparm *pbeparm;
 
 if (debug > 0)
 {
  Vnm_tprint( 1, "apbs:  Preparing to run %d PBE calculations.\n", nosh->ncalc);
 }
 
 /* for each focusing step call
    initMG, solveMG and save potential
 */
 for (int i=0; i<nosh->ncalc; i++)
 {
  if (debug > 0) Vnm_tprint( 1, "apbs:  CALCULATION #%d: MULTIGRID\n", i+1);

  /* Useful local variables */

  mgparm = nosh->calc[i]->mgparm;
  pbeparm = nosh->calc[i]->pbeparm;

  if (debug > 0) Vnm_tprint( 1, "apbs:    Setting up problem...\n");

  /** 
   initMG() - setups the problem. Calculates 
	or reinits accessibility functions and fills
	coefficient arrays.
  */

  if (setupflag == 1)
  {
   if (!initMG(i, nosh, mgparm, pbeparm, realCenter, pbe, &alist, dielX, dielY, dielZ, kappa, chargeMap, pmgp, pmg, dummyMap))
	 throw APBSErrorException("initMG", "APBS_runner::run", 0);
  }
  else if (setupflag >= 2)
  {
   if ( !initAPBScore(i, mgparm, pbeparm, realCenter) )
    throw APBSErrorException("APBS_runner::initAPBScore", "APBS_runner::run", 0);
  }
  
  /* print problem parameters */
  
  if (debug > 0) 
  {
   cout << "Printing mgparm:" << endl;
   printMGPARM(mgparm, realCenter);
   cout << "Printing pbeparm:" << endl;
	printPBEPARM(pbeparm);
  }

  /**
   solveMG() - solves (L)PBE
  */

  if (solveMG(nosh, pmg[i], mgparm->type) != 1)
   throw APBSErrorException("solveMG", "APBS_runner::run", 0);

  /* Rescue electrostatic potential. */
  
  Vgrid *pot = exportMap(pmg[i], pbeparm, VDT_POT);
  rescue_potential(pot);
  Vgrid_dtor(&pot);
 }
 
 cout << "Exiting APBS_runner::run." << endl;
}

void APBS_runner::rescue_potential(Vgrid *pot)
{
 /**
  Rescue potential calls the APBS method Vgrid_value()
  to calculate the electrostatic potential on all
  (off grid) atom positions. The result is stored in 
  'potential'.
 */

  double position[3];
  double val;
  
  while ( pqr->hasNext() )
  {
   Satom& atom = pqr->getNext();
	atom.coordinates().get_value(position);
	if( Vgrid_value(pot, position, &val) ) // atom is on this grid
	 potential.set_value_at(position[0], position[1], position[2], val);
  }
}

/* --------------------------------- SASA ------------------------- */

Accessible::Accessible(Structure *pqr, double epsw, double epsp, double srad, double maxIrad, double sdens)
{
 double nhash[3];
 int inhash[3];
 double lower_corner[3] = {0.0, 0.0, 0.0};
 double upper_corner[3] = {0.0, 0.0, 0.0};
 double radius = 0.0;

 // copy values
 this->solute  = pqr;
 this->epsp    = epsp;
 this->epsw    = epsw;
 this->srad    = srad;
 this->irad    = maxIrad;
 this->sdens   = sdens;

 // convert PQR into APBS alist
 alist = APBS_base::pqr2alist(pqr);

 /* build Vclist object */

 // determine Vclist hash size
 (pqr->getSoluteLength()*1.5).get_value( nhash );
 for(int i=0; i<3; i++)
 {
  inhash[i] = (int)nhash[i];
  if (inhash[i]<3) inhash[i] = 3;
  if (inhash[i]> MAX_HASH_DIM) inhash[i] = MAX_HASH_DIM;
 }

 // define maximum radius (important for the edges of Vclist)
 if( maxIrad > srad ) 
  radius = maxIrad + 0.5; // 0.5 from vpbe.c: MAX_SPLINE_WINDOW
 else
  radius = srad + 0.5; // 0.5 from vpbe.c: MAX_SPLINE_WINDOW

 // build Vclist object - expensive!
 clist = Vclist_ctor(alist, radius, inhash, CLIST_AUTO_DOMAIN, lower_corner, upper_corner);
 VASSERT(clist != VNULL);

 /* build SASA */
 acc = Vacc_ctor(alist, clist, sdens);
 VASSERT(acc != VNULL);
 
 // I'm an original
 this->copy = false;
}

void Accessible::mark_sphere(double rtot, double tpos[], Vgrid *grid, double markVal)
{
 int nx, ny, nz, i, j, k, imin, imax, jmin, jmax, kmin, kmax;
 double hx, hy, hzed, dx, dx2, dy, dy2, dz, dz2;
 double rtot2, pos[3];
 double *array;

 nx = grid->nx; hx = grid->hx;
 ny = grid->ny; hy = grid->hy;
 nz = grid->nz; hzed = grid->hzed;
 array = grid->data;

 /* Convert to grid reference frame */
 pos[0] = tpos[0] - grid->xmin;
 pos[1] = tpos[1] - grid->ymin;
 pos[2] = tpos[2] - grid->zmin;

 rtot2 = VSQR(rtot);

 dx = rtot + 0.5*hx;
 imin = VMAX2(0,(int)ceil((pos[0] - dx)/hx));
 imax = VMIN2(nx-1,(int)floor((pos[0] + dx)/hx));
 for (i=imin; i<=imax; i++) 
 {
  dx2 = VSQR(pos[0] - hx*i);
  if (rtot2 > dx2) 
  {
   dy = VSQRT(rtot2 - dx2) + 0.5*hy;
  } 
  else 
  {
   dy = 0.5*hy;
  }
  jmin = VMAX2(0,(int)ceil((pos[1] - dy)/hy));
  jmax = VMIN2(ny-1,(int)floor((pos[1] + dy)/hy));
  for (j=jmin; j<=jmax; j++) 
  {
   dy2 = VSQR(pos[1] - hy*j);
   if (rtot2 > (dx2+dy2)) 
	{ 
    dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed; 
   } 
	else 
	{
    dz = 0.5*hzed;
   }
   kmin = VMAX2(0,(int)ceil((pos[2] - dz)/hzed));
   kmax = VMIN2(nz-1,(int)floor((pos[2] + dz)/hzed));
   for (k=kmin; k<=kmax; k++) 
	{
    dz2 = VSQR(k*hzed - pos[2]);
    if ((dz2 + dy2 + dx2) <= rtot2) 
	 {
     array[IJK(i,j,k)] = markVal;
    }
   } /* k loop */
  } /* j loop */
 } /* i loop */
}

void Accessible::calculateDielMap(Vgrid *diel)
{
 int i, nx, ny, nz, iatom, ipt;
 double hx, hy, hzed, xmin, ymin, zmin, xmax, ymax, zmax, arad, position[3];

 /* Mesh info */
 nx   = diel->nx;
 ny   = diel->ny;
 nz   = diel->nz;
 hx   = diel->hx;
 hy   = diel->hy;
 hzed = diel->hzed;

 /* Define the min/max dimensions */
 xmin = diel->xmin;
 ymin = diel->ymin;
 zmin = diel->zmin;
 xmax = diel->xmax;
 ymax = diel->ymax;
 zmax = diel->zmax;

 /* Reset the arrays */
 for (i=0; i<(nx*ny*nz); i++) 
  diel->data[i] = epsw;

 /* loop over all atoms an mark atomic radii */
 while ( solute->hasNext() ) 
 {
  Satom& atom = solute->getNext();
  atom.coordinates().get_value(position);
  arad = atom.get_radius();

  /* Make sure we're on the grid */
  if ((position[0]>xmin) || (position[0]<xmax)  || (position[1]>ymin) || (position[1]<ymax)  || (position[2]>zmin) || (position[2]<zmax))
  {
   if (arad > VSMALL) mark_sphere((arad+srad), position, diel, epsp);
  }
 }

 /* We only need to do the next step for non-zero solvent radii */
 if (srad > VSMALL) 
 {
  /* Now loop over the solvent accessible surface points */
  iatom = 0;
  while ( solute->hasNext() ) 
  {
   solute->getNext();
   Vatom *atom = Valist_getAtom(alist, iatom++);
   VaccSurf *asurf = Vacc_atomSASPoints(acc, srad, atom);
   
   /* Use each point on the SAS to reset the solvent accessibility */
   for (ipt=0; ipt<(asurf->npts); ipt++)
	{
    position[0] = asurf->xpts[ipt];
    position[1] = asurf->ypts[ipt];
    position[2] = asurf->zpts[ipt];
    mark_sphere(srad, position, diel, epsw);
   }
  }
 }

 /* Now we run over all solute atoms and 
    mark grid points with user supplied dielectric constant 
	 values */

 while( solute->hasNext() )
 {
  Satom& atom  = solute->getNext();
  double pdiel = atom.get_soluteDiel();
  if (pdiel != epsp)
  {
   atom.coordinates().get_value(position);
	arad = atom.get_radius();
   mark_sphere(arad, position, diel, pdiel);
  }
 }
}

void Accessible::calculateKappaMap(double ionstr, Vgrid *kappa)
{
 int i, nx, ny, nz;
 double hx, hy, hzed, xmin, ymin, zmin, xmax, ymax, zmax, arad, ionmask, position[3];

 /* Mesh info */
 nx   = kappa->nx;
 ny   = kappa->ny;
 nz   = kappa->nz;
 hx   = kappa->hx;
 hy   = kappa->hy;
 hzed = kappa->hzed;

 /* Define the min/max dimensions */
 xmin = kappa->xmin;
 ymin = kappa->ymin;
 zmin = kappa->zmin;
 xmax = kappa->xmax;
 ymax = kappa->ymax;
 zmax = kappa->zmax;

 if(ionstr > VPMGSMALL)
  ionmask = 1.0; /* ion accessible */
 else
  ionmask = 0.0; /* ion inaccessible */

 /* Reset the arrays */
 for (i=0; i<(nx*ny*nz); i++) 
  kappa->data[i] = ionmask;

 /* exit if ionic strength is zero */
 if (ionstr < VPMGSMALL) return;

 /* loop over all atoms an mark atomic radii */
 while ( solute->hasNext() )
 {
  Satom& atom = solute->getNext();
  atom.coordinates().get_value(position);
  arad = atom.get_radius();

  /* Make sure we're on the grid */
  if ((position[0]>xmin) || (position[0]<xmax)  || \
      (position[1]>ymin) || (position[1]<ymax)  || \
      (position[2]>zmin) || (position[2]<zmax)) 
  {
   if (arad > VSMALL) 
	{
    /* Mark kappa map */
    mark_sphere((arad+irad), position, kappa, 0.0);
   }
  } /* endif (on the mesh) */
 } /* endwhile (over all atoms) */
 //Vgrid_writeDX(kappa, "FILE", "ASC", VNULL, "debug.dx", "", NULL);
}

void Accessible::smoothMaps(Vgrid *dielX, Vgrid *dielY, Vgrid *dielZ)
{
 // perform harmonic averaging smoothing
 double frac, *a1cf, *a2cf, *a3cf;
 int i, j, k, nx, ny, nz, numpts;

 /* Mesh info */
 nx = dielX->nx;
 ny = dielX->ny;
 nz = dielX->nz;

 /* allocate space for work arrays */
 a1cf = new double[nx*ny*nz];
 a2cf = new double[nx*ny*nz];
 a3cf = new double[nx*ny*nz];

 /* Copy the existing diel arrays to work arrays */
 for (i=0; i<(nx*ny*nz); i++) 
 {
  a1cf[i] = dielX->data[i];
  a2cf[i] = dielY->data[i];
  a3cf[i] = dielZ->data[i];
  /* dielX->data[i] = epsw;
  dielY->data[i] = epsw;
  dielZ->data[i] = epsw; */
 }

 /* Smooth the dielectric values */
 for (i=0; i<nx; i++) 
 {
  for (j=0; j<ny; j++) 
  {
   for (k=0; k<nz; k++) 
	{
    /* Get the 8 points that are 1/sqrt(2) grid spacings away */

    /* Points for the X-shifted array */
    frac = 1.0/a1cf[IJK(i,j,k)];
    frac += 1.0/a2cf[IJK(i,j,k)];
    frac += 1.0/a3cf[IJK(i,j,k)];
    numpts = 3;
    if (j > 0) 
	 {
     frac += 1.0/a2cf[IJK(i,j-1,k)];
     numpts += 1;
    } 
    if (k > 0) 
	 {
     frac += 1.0/a3cf[IJK(i,j,k-1)];
     numpts += 1;
    }
    if (i < (nx-1))
	 {
     frac += 1.0/a2cf[IJK(i+1,j,k)];
     frac += 1.0/a3cf[IJK(i+1,j,k)];
     numpts += 2;
     if (j > 0) 
	  {
      frac += 1.0/a2cf[IJK(i+1,j-1,k)];
      numpts += 1;
     }
     if (k > 0) 
	  {
      frac += 1.0/a3cf[IJK(i+1,j,k-1)];
      numpts += 1;
     }
    }
    dielX->data[IJK(i,j,k)] = numpts/frac;
             
    /* Points for the Y-shifted array */
    frac = 1.0/a2cf[IJK(i,j,k)];
    frac += 1.0/a1cf[IJK(i,j,k)];
    frac += 1.0/a3cf[IJK(i,j,k)];
    numpts = 3;

    if (i > 0) 
	 {
     frac += 1.0/a1cf[IJK(i-1,j,k)];
     numpts += 1;
    } 
    if (k > 0) 
	 {
     frac += 1.0/a3cf[IJK(i,j,k-1)];
     numpts += 1;
    }
    if (j < (ny-1))
	 {
     frac += 1.0/a1cf[IJK(i,j+1,k)];
     frac += 1.0/a3cf[IJK(i,j+1,k)];
     numpts += 2;
     if (i > 0) 
	  {
      frac += 1.0/a1cf[IJK(i-1,j+1,k)];
      numpts += 1;
     }
     if (k > 0) 
	  {
      frac += 1.0/a3cf[IJK(i,j+1,k-1)];
      numpts += 1;
     }
    }
    dielY->data[IJK(i,j,k)] = numpts/frac;

    /* Points for the Z-shifted array */
    frac = 1.0/a3cf[IJK(i,j,k)];
    frac += 1.0/a1cf[IJK(i,j,k)];
    frac += 1.0/a2cf[IJK(i,j,k)];
    numpts = 3;

    if (i > 0) 
	 {
     frac += 1.0/a1cf[IJK(i-1,j,k)];
     numpts += 1;
    }
    if (j > 0) 
	 {
     frac += 1.0/a2cf[IJK(i,j-1,k)];
     numpts += 1;
    }
    if (k < (nz-1))
	 {
     frac += 1.0/a1cf[IJK(i,j,k+1)];
     frac += 1.0/a2cf[IJK(i,j,k+1)];
     numpts += 2;
     if (i > 0) 
	  {
      frac += 1.0/a1cf[IJK(i-1,j,k+1)];
      numpts += 1;
     }
     if (j > 0) 
	  {
      frac += 1.0/a2cf[IJK(i,j-1,k+1)];
      numpts += 1;
     }
    }
    dielZ->data[IJK(i,j,k)] = numpts/frac;
   }
  }
 }
 
 /* clean up */
 delete[] a1cf;
 delete[] a2cf;
 delete[] a3cf;
}

/* Membrane slab */

void MembraneAccessible::fillSlab(const Gridpoint& gp, const vector<Vgrid*>& grid, double zmin, double zmax, double fcol, double bcol) 
{
 stack<Gridpoint> s;
 vector<int> kminv;
 vector<int> kmaxv;
 int ig, in, kmin, kmax;
  
 /* transform slab zmin/zmax into kmin/kmax */
 for(ig=0; ig<grid.size(); ig++) 
 {
  kminv.push_back( VMAX2(0, (int)ceil((zmin - grid[ig]->zmin - 0.5*grid[ig]->hzed)/grid[ig]->hzed)) );
  kmaxv.push_back( VMIN2(grid[ig]->nz-1, (int)floor((zmax - grid[ig]->zmin + 0.5*grid[ig]->hzed)/grid[ig]->hzed)) );
 }

 /* apply boundary fill algorithm on all grids */
 s.push(gp); 
 while( !s.empty() ) 
 {
  //cout << "Stack size = " << s.size() << endl;
  Gridpoint pt = s.top(); s.pop();
  if( (pt.getVal()!=bcol) && (pt.getVal()!=fcol) ) 
  {
   for(ig=0; ig<grid.size(); ig++) 
   {
    if( pt.isContained(grid[ig]) ) 
    {
     kmin = kminv[ig];
     kmax = kmaxv[ig];
     vector<Gridpoint> neighbor = pt.transformed(grid[ig]).getNeighbors(0, grid[ig]->nx-1, 0, grid[ig]->ny-1, kmin, kmax);
     for(in=0; in<neighbor.size(); in++)
      s.push(neighbor[in]);
    }
   }
   pt.setVal(fcol);
  }
 }
/* return */
}

MapController::MapController(APBSparms *apbsparms)
{
 int dime[3], n;
 int fsteps   = apbsparms->get_nfsteps();
 int nions    = apbsparms->get_nions();
 
 double *ionq = new double[nions];
 double *ionc = new double[nions];
 double *ionr = new double[nions];

 // init grids
 for(int i=0; i<fsteps; i++)
 {
  gridCenterX.push_back(0.0);
  gridCenterY.push_back(0.0);
  gridCenterZ.push_back(0.0);
  apbsparms->get_dime(i, dime);
  n = dime[0]*dime[1]*dime[2];
  dielx.push_back( new double[n] );
  diely.push_back( new double[n] );
  dielz.push_back( new double[n] );
  kappa.push_back( new double[n] );
 }

 // calculate max ion radius and ionic strength
 ionMaxRad = 0.0;
 ionStrength = 0.0;
 apbsparms->get_ion(ionq, ionc, ionr);
 for(int i=0; i<nions; i++)
 {
  if(ionr[i]>ionMaxRad) ionMaxRad = ionr[i];
  ionStrength += 0.5*ionc[i]*ionq[i]*ionq[i];
 }

 cout << "MapController: Max. ion radius: " << ionMaxRad << endl;
 cout << "MapController: Ionic strength: " << ionStrength << endl;

 // store apbsparms
 this->apbsparms = apbsparms;

 // initialize molecular accessiblity
 smol = apbsparms->get_srfm();
 acc = NULL;

 delete[] ionq;
 delete[] ionc;
 delete[] ionr;
}

void MapController::fillMapsMol(const Structure& pqr, double centerx[], double centery[], double centerz[], Vgrid **dielX, Vgrid **dielY, Vgrid **dielZ, Vgrid **kappa)
{
 bool fillFlag = false;

 // if solute has been changed recalculate accessibility
 if( pqr!=solute )
 {
  cout << "MapController::fillMapsMol: Calculating molecular accessibility ... ";
  cout.flush();
  solute = pqr;
  if (acc!=NULL) delete acc;
  acc = new Accessible(&solute, apbsparms->get_sdie(), apbsparms->get_pdie(), apbsparms->get_srad(), ionMaxRad, apbsparms->get_sdens());
  fillFlag = true; /* all maps will have to be refilled */
  cout << "finished." << endl;
 }

 // loop over all focusing steps and check whether maps have to be refilled
 for(int i=0; i<apbsparms->get_nfsteps(); i++)
 {
  double center[] = {centerx[i], centery[i], centerz[i]};

  // create new Vgrid object
  dielX[i] = getMap(i, center, VDT_DIELX);
  dielY[i] = getMap(i, center, VDT_DIELY);
  dielZ[i] = getMap(i, center, VDT_DIELZ);
  kappa[i] = getMap(i, center, VDT_KAPPA);

  // refill maps if grid center changed or solute has changed
  if( gridCenterHasChanged(i, center) || fillFlag )
  {
   cout << "MapController::fillMapsMol: Refilling accessibility maps in step " << i+1 << " ... ";
	cout.flush();
   acc->calculateDielMap( dielX[i] );
   acc->calculateDielMap( dielY[i] );
   acc->calculateDielMap( dielZ[i] );
	if (smol==VSM_MOLSMOOTH) Accessible::smoothMaps(dielX[i], dielY[i], dielZ[i]);
   acc->calculateKappaMap( ionStrength, kappa[i] );
	cout << "finished." << endl;
  }
  else
  {
   cout << "MapController::fillMapsMol: Nothing to do in step " << i+1 << "." << endl;
	cout.flush();
  }
  
  // write maps to files
  if(dx_filename!="")
  {
   cout << "MapController::fillMapsMol: Writing maps in step " << i+1 << " to file ... ";
   cout.flush();

   ostringstream ost;
   ost << dx_filename << ".dielx." << i+1 << ".dx";
   Vgrid_writeDX(dielX[i], "FILE", "ASC", VNULL, ost.str().c_str(), "", NULL);
  
   ost.str("");
   ost << dx_filename << ".diely." << i+1 << ".dx";
   Vgrid_writeDX(dielY[i], "FILE", "ASC", VNULL, ost.str().c_str(), "", NULL);

   ost.str("");
   ost << dx_filename << ".dielz." << i+1 << ".dx";
   Vgrid_writeDX(dielZ[i], "FILE", "ASC", VNULL, ost.str().c_str(), "", NULL);

   ost.str("");
   ost << dx_filename << ".kappa." << i+1 << ".dx";
   Vgrid_writeDX(kappa[i], "FILE", "ASC", VNULL, ost.str().c_str(), "", NULL);

   cout << "finished." << endl;
  }
 }
}

void MembraneMapController::fillMapsMol(const Structure& pqr, double centerx[], double centery[], double centerz[], Vgrid **dielX, Vgrid **dielY, Vgrid **dielZ, Vgrid **kappa)
{
 bool fillFlag = false;

 // if solute has been changed recalculate accessibility
 if( pqr!=solute )
 {
  cout << "MembraneMapController::fillMapsMol: Calculating molecular accessibility ... ";
  cout.flush();
  solute = pqr;
  if (acc!=NULL) delete acc;
  acc = new MembraneAccessible(&solute, apbsparms->get_sdie(), apbsparms->get_pdie(), epsm, apbsparms->get_srad(), ionMaxRad, apbsparms->get_sdens(), zmin, zmax);
  cout << "finished." << endl;
  fillFlag = true;
 }

 /* Loop over all focusing steps and refill maps.
    In case of membrane slab calculations we have to refill 
	 all maps for each titratable group because 
	 Boundary Fill has to be run for all grids starting with coarsest
	 down to the finest grid. However, we do this only once for reference
	 state charge set.
 */
 
 vector<Vgrid*> vdielx;
 vector<Vgrid*> vdiely;
 vector<Vgrid*> vdielz;
 vector<Vgrid*> vkappa;
 for(int i=0; i<apbsparms->get_nfsteps(); i++)
 {
  double center[] = {centerx[i], centery[i], centerz[i]};
  // create new Vgrid objects
  vdielx.push_back( getMap(i, center, VDT_DIELX) );
  vdiely.push_back( getMap(i, center, VDT_DIELY) );
  vdielz.push_back( getMap(i, center, VDT_DIELZ) );
  vkappa.push_back( getMap(i, center, VDT_KAPPA) );
  dielX[i] = vdielx[i];
  dielY[i] = vdiely[i];
  dielZ[i] = vdielz[i];
  kappa[i] = vkappa[i];
  if( gridCenterHasChanged(i, center) )
   fillFlag = true;
 }

 /* refill all maps if solute/solvent 
    boundary or grid center of
	 any of the grids have changed.
 */
 
 if( fillFlag )
 {
  cout << "MembraneMapController::fillMapsMol: Refilling accessibility maps ... ";
  cout.flush();
  ((MembraneAccessible*)acc)->calculateDielMap( vdielx );
  ((MembraneAccessible*)acc)->calculateDielMap( vdiely );
  ((MembraneAccessible*)acc)->calculateDielMap( vdielz );
  ((MembraneAccessible*)acc)->calculateKappaMap( ionStrength, vkappa );
  for(int i=0; i<apbsparms->get_nfsteps(); i++)
  {
   if (smol==VSM_MOLSMOOTH) 
    Accessible::smoothMaps(dielX[i], dielY[i], dielZ[i]);
  }
  cout << "finished." << endl;
 }
 else
 {
  cout << "MembraneMapController::fillMapsMol: Nothing to do." << endl;
  cout.flush();
 }
 
 // debug
 if( (dx_filename!="") && fillFlag )
 {
  for(int i=0; i<apbsparms->get_nfsteps(); i++)
  {
   cout << "MapController::fillMapsMol: Writing maps in step " << i+1 << " to file ... ";
   cout.flush();

   ostringstream ost;
   ost << dx_filename << ".dielx." << i+1 << ".dx";
   Vgrid_writeDX(dielX[i], "FILE", "ASC", VNULL, ost.str().c_str(), "", NULL);
  
   ost.str("");
   ost << dx_filename << ".diely." << i+1 << ".dx";
   Vgrid_writeDX(dielY[i], "FILE", "ASC", VNULL, ost.str().c_str(), "", NULL);

   ost.str("");
   ost << dx_filename << ".dielz." << i+1 << ".dx";
   Vgrid_writeDX(dielZ[i], "FILE", "ASC", VNULL, ost.str().c_str(), "", NULL);

   ost.str("");
   ost << dx_filename << ".kappa." << i+1 << ".dx";
   Vgrid_writeDX(kappa[i], "FILE", "ASC", VNULL, ost.str().c_str(), "", NULL);

   cout << "finished." << endl;
  }
 }
}
