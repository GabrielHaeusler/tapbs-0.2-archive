#ifndef APBSWRAPPER_H
#define APBSWRAPPER_H 1
#define APBS4

#include "parser.h"
#include "potential.h"

/* APBS stuff */
#include "apbs/apbs.h"
#include "apbs/nosh.h"
#include "apbs/mgparm.h"
#include "apbs/pbeparm.h"
#include "apbs/femparm.h"
#include "apbs/routines.h" 

#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define ON_PROTEIN "oncent"
#define ON_MODEL "onsite"

/**
 * @brief Store user-provided APBS parameters.
 * Stores all APBS relevant user-defined parameters
 * Parses all user-provided APBS parameters from parser
 * objects.
 */

class APBSparms
{
protected:

 /* properties */

 // general LPBE/MG parameters
 map<string, double> rparm;                              /** double-value parameters */
 map<string, int>    iparm;                              /** integer-value parameters */
 vector< vector<double> > ion;                           /** ion species charges, conc, radii */

 // focusing specific parameters
 vector< vector< vector<int> > > dime;                   /** grid points, protein/model->fstep->nx/ny/nz */
 vector< vector< vector<double> > > spacing;             /** grid spacing, protein/model->fstep->hx/hy/hz */
 vector< vector<string> > centering;                     /** grid centering method, protein/model->fstep */

 /* methods */

 /**
  * @brief Gets all user provided parameters from parser.
  * @author Gernot Kieseritzky
  */

 void init(InputParser& parser);

public:

 /* methods */

 /**
  * @brief Defines possible parser key words.
  * @author Gernot Kieseritzky
  */

 static void getInputDescription(InputDescription& inputDescription);
 
 /**
  * @brief Standard constructor.
  * @author Gernot Kieseritzky
  */

 APBSparms(InputParser& parser);

 /**
  * @brief Copy constructor.
  * @author Gernot Kieseritzky
  */

 APBSparms(const APBSparms& apbsparms);

 /**
  * @brief Minimum numer of multigrid levels?
  * @author Gernot Kieseritzky
  * @returns Integer value that represents minimum number of multgrid levels. (4)
  */

 int get_nlev() {return iparm["nlev"];}

 /**
  * @brief Number of ions.
  * @author Gernot Kieseritzky
  * @returns Integer value that represents number of ions.
  */

 int get_nions()   {return iparm["nions"];}

 /**
  * @brief Get default protein dielectric constant. 
  * @author Gernot Kieseritzky
  * @returns Double value that represents protein dielectric constant.
  */

 double get_pdie() {return rparm["pdie"];}

 /**
  * @brief Get default solvent dielectric constant. 
  * @author Gernot Kieseritzky
  * @returns Double value that represents solvent dielectric constant.
  */

 double get_sdie() {return rparm["sdie"];}

 /**
  * @brief Get temperature (K).
  * @author Gernot Kieseritzky
  * @returns Double value that represents temperature.
  */

 double get_temperature() {return rparm["temperature"];}

 /**
  * @brief Get solvent probe radius (A).
  * @author Gernot Kieseritzky
  * @returns Double value that represents the solvent probe radius in A.
  */

 double get_srad() {return rparm["srad"];}

 /**
  * @brief Get accuracy of SASA. (Point densitiy of SASA spheres).
  * @author Gernot Kieseritzky
  * @note Only APBS 0.4.0.
  * @returns Double value that represents the SASA accuracy in A.
  */

 double get_sdens(){return rparm["sdens"];}

 /**
  * @brief Multgrid error tolerance.
  * By default error is defined in terms of the relative residual.
  * APBS sets this to 1E-6 by default. But 5E-2 seems to work as well
  * and is faster by a factor > 2.
  * @author Gernot Kieseritzky
  * @returns Double value that represents the error tolerance.
  */

 double get_errtol() {return rparm["errtol"];}

 /**
  * @brief Get APBS blab level.
  * Level of details in io.mc
  * @author Gernot Kieseritzky
  * @returns Integer value that represents the blab level.
  */

 int get_iinfo() {return iparm["iinfo"];}

 /**
  * @brief Get number of presmoothing steps.
  * @author Gernot Kieseritzky
  * @returns Integer value that represents the number of presmoothing steps.
  */

 int get_presmooth() {return iparm["presmooth"];}

 /**
  * @brief Get number of postsmoothing steps.
  * @author Gernot Kieseritzky
  * @returns Integer value that represents the number of presmoothing steps.
  */

 int get_postsmooth() {return iparm["postsmooth"];}

 /**
  * @brief Get coarsening method
  * @author Gernot Kieseritzky
  * @returns Integer value that represents the coarsening method (See APBS doc).
  */

 int get_mgcoar() {return iparm["mgcoar"];}

 /**
  * @brief Get prolongation method
  * @author Gernot Kieseritzky
  * @returns Integer value that represents the prolongation method (See APBS doc).
  */

 int get_mgprol() {return iparm["mgprol"];}

 /**
  * @brief Get smoothing method
  * @author Gernot Kieseritzky
  * @returns Integer value that represents the smoothing method (See APBS doc).
  */

 int get_mgsmoo() {return iparm["mgsmoo"];}

 /**
  * @brief Get coarse grid solving method
  * @author Gernot Kieseritzky
  * @returns Integer value that represents coarse grid solving method (See APBS doc).
  */

 int get_mgsolv() {return iparm["mgsolv"];}

 /**
  * @brief Get maximum number of iterations
  * @author Gernot Kieseritzky
  * @returns Integer value that represents max. iterations.
  */

 int get_itmax() {return iparm["itmax"];}

 /**
  * @brief Get MG method
  * @author Gernot Kieseritzky
  * @returns Integer value that represents the MG method (0 - Vcycle, 1 - nested iteration).
  */

 int get_mgkey() {return iparm["mgkey"];}

 /**
  * @brief Get preconditioning method
  * @author Gernot Kieseritzky
  * @returns Integer value that represents preconditioning method (see APBS doc).
  */

 int get_ipcon() {return iparm["ipcon"];}

 /**
  * @brief Get type of operator analysis
  * @author Gernot Kieseritzky
  * @returns Integer value that represents operator analysis type (see APBS doc).
  */

 int get_iperf() {return iparm["iperf"];}

 /**
  * @brief Get boundary condition for coarse grids.
  * BCFL_ZERO -> zero potential on boundary grid points
  * BCFL_SDH  -> Single debye-hueckel boundary condition.
  * BCFL_MDH  -> Multiple debye-hueckel boundary condition.
  * @author Gernot Kieseritzky
  * @returns Enumeration that represents the type of boundary condition.
  */

 Vbcfl get_bcfl()
 { 
  switch(iparm["bcfl"])
  {
   case 0: return BCFL_ZERO; break;
   case 1: return BCFL_SDH;  break;
   case 2: return BCFL_MDH;  break;
  }
  return BCFL_SDH;
 }

 /**
  * @brief Get charge discretisation method.
  * VCM_TRIL -> trilinear interpolation.
  * VCM_BSPL2  -> cubic B-spline interpolation.
  * @author Gernot Kieseritzky
  * @returns Enumeration that represents the charge discr. method.
  */

 Vchrg_Meth get_chgm()
 {
  switch(iparm["chgm"])
  {
   case 0: return VCM_TRIL; break;
   case 1: return VCM_BSPL2; break;
  }
  return VCM_TRIL;
 }

 /**
  * @brief Get ion/solvent accessibility coefficient model.
  * VSM_MOL        -> Molecular surface.
  * VSM_MOLSMOOTH  -> Molecular surface with harmonic-average smoothing.
  * VSM_SPLINE     -> Spline based molecular surface.
  * @author Gernot Kieseritzky
  * @returns Enumeration that represents the coefficient model.
  */

 Vsurf_Meth get_srfm()
 {
  switch(iparm["srfm"])
  {
   case 0: return VSM_MOL; break;
   case 1: return VSM_MOLSMOOTH; break;
   case 2: return VSM_SPLINE; break;
  }
  return VSM_MOL;
 }

 /**
  * @brief Copy ion data into arrays.
  * @author Gernot Kieseritzky
  */

 void get_ion(double ionq[], double ionc[], double ionr[]);

 /**
  * @brief Copy number of grid points in focusing step 'fstep' into array.
  * @author Gernot Kieseritzky
  */

 virtual void get_dime(int fstep, int array[]);

 /**
  * @brief Copy grid resolution in focusing step 'fstep' into array.
  * @author Gernot Kieseritzky
  */

 virtual void get_grid(int fstep, double array[]);

 /**
  * @brief Get grid centering method in focusing step 'fstep'.
  * @author Gernot Kieseritzky
  * @return Reference to string that represents the grid centering method. ('oncent' or 'onsite')
  */

 virtual string get_centering(int fstep);

 /**
  * @brief Get number of focusing steps, i.e. the number of calculations.
  * @author Gernot Kieseritzky
  * @return Integer value that represents the number of calculations.
  */

 virtual int get_nfsteps() {return iparm["nfsteps_protein"];}

 /**
  * @brief Standard destructor.
  * @author Gernot Kieseritzky
  */

 virtual ~APBSparms() {cout << "Good bye from APBSparms!" << endl;}
};

/**
 * @brief Specialization of APBSparms for model compound calculations.
 * This specialization of APBSparms is needed for model 
 * calculation runs with a smaller number of grid points. 
 * It simply overrides some interface functions that return
 * varius grid related parameters.
*/

class APBSparms_model : public APBSparms
{
public:
 APBSparms_model(InputParser& parser) : APBSparms(parser)
 {
  // do nothing
 }

 APBSparms_model(const APBSparms& apbsparms) : APBSparms(apbsparms)
 {
  // do nothing
 }

 virtual void get_dime(int fstep, int array[]);

 virtual void get_grid(int fstep, double array[]);

 virtual string get_centering(int fstep);

 virtual int get_nfsteps() {return iparm["nfsteps_model"];}
};

/* ----------------------- SASA calculations using APBS ------------------------------ */

/**
 * @brief Calculates SASA and ion accessibility maps using APBS routines.
 */

class Accessible
{
protected:
 Structure *solute; /** Dielectric volume */

 Valist *alist;     /** APBS Valist object */
 Vclist *clist;     /** APBS Vclist object */
 Vacc *acc;         /** APBS Vacc object -> SASA */

 double epsp;       /** solute dielectric constant */
 double epsw;       /** solvens dielectric constant */
 double srad;       /** solvent probe radius */
 double irad;       /** ion radius */
 double sdens;      /** SASA point density */

 bool copy;         /** flag, whether this instance is a copy */

 /* methods */

 /**
  * @brief Marks spherical region on grid.
  * Shamelessly ripped from APBS source code with minor modifications.
  * @param rtot Total radius of sphere.
  * @param tpos Coordinates of sphere (tpos[0..2]).
  * @param grid instance of Vgrid containing data array and grid geometry
  * @param markVal marking value
  * @author Nathan Baker
  */

 void mark_sphere(double rtot, double tpos[], Vgrid *grid, double markVal);

 void copyMembers(const Accessible& that)
 {
  this->solute = that.solute;
  this->epsp   = that.epsp;
  this->epsw   = that.epsw;
  this->srad   = that.srad;
  this->irad   = that.irad;
  this->sdens  = that.sdens;

  /* flat copying */
  if (!copy) // I'm not a copy
  {
   Vclist_dtor(&clist);
   Vacc_dtor(&acc);
	Valist_dtor(&alist);
  }
  this->clist  = that.clist;
  this->acc    = that.acc;
  this->alist  = that.alist;
  
  /* mark me to be a copy */
  this->copy   = true;
 }

public:

 /**
  * @brief Constructor.
  * @param pqr instance of Structure containing atom specific dielectric constants.
  * @param epsw solvent dielectric constant
  * @param epsp solute dielectric constant
  * @param srad solvent probe radius.
  * @param maxIrad largest ion radius.
  * @param sdens Point density of SASA.
  * @author Gernot Kieseritzky
  */

 Accessible(Structure *pqr, double epsw, double epsp, double srad, double maxIrad, double sdens);

 /**
  * @brief Copy constructor.
  * Will only make a flat copy, i.e. only the pointer to Vclist/Vacc-objects
  * will be copied to save memory.
  * @author Gernot Kieseritzky
  */

 Accessible(const Accessible& that) {copyMembers(that);}

 /**
  * @brief Calculate dielectric constant map.
  * Shamelessly ripped from APBS source code with minor modifications.
  * @param diel Vgrid object containing array to be filled.
  * @author Nathan Baker
  */

 void calculateDielMap(Vgrid *diel);

 /**
  * @brief Calculate kappa map.
  * Shamelessly ripped from APBS source code with minor modifications.
  * @param ionstr bulk ionic strength
  * @param kappa Vgrid object containing array to be filled.
  * @author Nathan Baker
  */

 void calculateKappaMap(double ionstr, Vgrid *kappa);

 /**
  * @brief Perform 9-point harmonic average smoothing.
  * Shamelessly ripped from APBS source code with minor modifications.
  * @param dielX x-shifted dielectric map.
  * @param dielY y-shifted dielectric map.
  * @param dielZ z-shifted dielectric map.
  * @author Nathan baker.
  */

 static void smoothMaps(Vgrid *dielX, Vgrid *dielY, Vgrid *dielZ);

 /**
  * @brief Assignment operator.
  * Will make a flat copy.
  * @param that Original to be copied.
  * @author Gernot Kieseritzky
  */

 Accessible& operator=(const Accessible& that)
 {
  if (this!=&that)
   copyMembers(that);

  return *this;
 }

 /**
  * @brief Destructor.
  * Will only destroy Vclist/Vacc objects if we are not a copy.
  * @author Gernot Kieseritzky
  */

 virtual ~Accessible() 
 {
  if(!copy)
  {
   Vclist_dtor(&clist);
   Vacc_dtor(&acc);
   Valist_dtor(&alist);
  }
 }
};

/**
 * @brief Calculates SASA and ion accessibility maps using APBS routines 
 * and includes an implicit description of a membrane slab.
 */

class MembraneAccessible : public Accessible
{
 protected:
 double epsm, zmin, zmax; /** Membrane dielectric constant and dimensions */

/**
 * @brief Wraps gridpoint information and related methods.
 */

 class Gridpoint 
 {
  protected:
  Vgrid *grid;       /** the grid */
  int i, j, k, n;    /** discrete coordinates, n is the index */

 /**
  * @brief Helper function returning the array index of the grid point at (i,j,k).
  * @param i discrete x coordinate of grid point in grid reference frame.
  * @param j discrete z coordinate of grid point in grid reference frame.
  * @param k discrete z coordinate of grid point in grid reference frame.
  * @return Integer value representing index in array of grid points.
  * @author Gernot Kieseritzky
  */
 
  int ijk(int i, int j, int k) const 
  {
   return k*grid->nx*grid->ny + j*grid->nx + i;
  }

  public:

 /**
  * @brief Constructor 1
  * @param g point to Vgrid object.
  * @param xg x coordinate of grid point in grid reference frame.
  * @param yg y coordinate of grid point in grid reference frame.
  * @param zg z coordinate of grid point in grid reference frame.
  * @author Gernot Kieseritzky
  */

  Gridpoint(Vgrid *g, double xg, double yg, double zg) 
  {
   grid = g;
   i = (int)(xg/g->hx + 0.5*g->hx);
   j = (int)(yg/g->hy + 0.5*g->hy);
   k = (int)(zg/g->hzed + 0.5*g->hzed);
   n = ijk(i,j,k);
  }

 /**
  * @brief Constructor 2
  * @param g point to Vgrid object.
  * @param i discrete x-coordinate of grid point in grid reference frame.
  * @param j discrete y-coordinate of grid point in grid reference frame.
  * @param k discrete z-coordinate of grid point in grid reference frame.
  * @author Gernot Kieseritzky
  */
  
  Gridpoint(Vgrid *g, int i, int j, int k) 
  {
   grid = g;
   this->i = i;
   this->j = j;
   this->k = k;
   n = ijk(i,j,k);
  }

 /**
  * @brief Changes the value of grid point. (e.g. the dielectric constant)
  * @param val value.
  * @author Gernot Kieseritzky
  */
  
  void setVal(double val) 
  {
   grid->data[n] = val;
  }

 /**
  * @brief Returns the value of grid point. (e.g. the dielectric constant)
  * @return Double value.
  * @author Gernot Kieseritzky
  */

  double getVal() const {
   return grid->data[n];
  }

 /**
  * @brief Return the grid neighbors within specified boundaries.
  * @param imin Minimim discrete x-coordinate.
  * @param imax Maximum discrete x-coordinate.
  * @param jmin Minimim discrete y-coordinate.
  * @param jmax Maximum discrete y-coordinate.
  * @param kmin Minimim discrete z-coordinate.
  * @param kmax Maximum discrete z-coordinate.
  * @return Vector container filled with neighboring grid points.
  * @author Gernot Kieseritzky
  */
  
  vector<Gridpoint> getNeighbors(int imin, int imax, int jmin, int jmax, int kmin, int kmax) const {
   int in, jn, kn;
   vector<Gridpoint> neighbors;
   
   /* find neighbors */
   for(int ip=-1; ip<=1; ip++) {
    in = i+ip;
    if( (in<imin) || (in>imax) ) continue;
    for(int jp=-1; jp<=1; jp++) {
     jn = j+jp;
     if( (jn<jmin) || (jn>jmax) ) continue;
     for(int kp=-1; kp<=1; kp++) {
      kn = k+kp;
      if( (kn<kmin) || (kn>kmax) ) continue;
      if( (ip==0) && (jp==0) && (kp==0) ) continue;
      neighbors.push_back( Gridpoint(grid, in, jn, kn) );
     }
    }
   }
   
   return neighbors;
  }

 /**
  * @brief Returns true if grid point is located on the other grid in question.
  * @param g The other grid in question.
  * @return 0 - False, 1 - True
  * @author Gernot Kieseritzky
  */
    
  int isContained(Vgrid *g) const 
  {
   /* g points to the same address as grid */
   if(grid==g)
	 return 1;
   
   /* more work if pointers are not equal */
   /* change reference frame */
   double x = i * grid->hx   + grid->xmin;
   double y = j * grid->hy   + grid->ymin;
   double z = k * grid->hzed + grid->zmin;

   if( (x>=g->xmin) && (x<=g->xmax) && (y>=g->ymin) && (y<=g->ymax) && (z>=g->zmin) && (z<=g->zmax) )
    return 1;
   else
    return 0;
  }

 /**
  * @brief Transform gridpoint coordinates into other's grid reference frame
  * @param g The other grid.
  * @return Gridpoint object with transformed coordinates
  * @author Gernot Kieseritzky
  */

  Gridpoint transformed(Vgrid *g) const 
  {
   /* g points to the same address as grid */
   if(grid==g)
	 return *this;

   /* more work if pointers are not equal */
   /* change reference frame */
   double xp = i * grid->hx   + grid->xmin - g->xmin;
   double yp = j * grid->hy   + grid->ymin - g->ymin;
   double zp = k * grid->hzed + grid->zmin - g->zmin;
   
   return Gridpoint(g, xp, yp, zp);
  }
 /* class Gridpoint */
 };

 public:
 
 /**
 * @brief Constructor.
 * @param g The other grid.
 * @author Gernot Kieseritzky
 */
 
 MembraneAccessible(Structure *pqr, double epsw, double epsp, double epsm, double srad, double maxIrad, double sdens, double zmin, double zmax) : Accessible(pqr, epsw, epsp, srad, maxIrad, sdens) 
 {
  this->epsm = epsm;
  this->zmin = zmin;
  this->zmax = zmax;
 }

 /**
 * @brief Fill dielectric map including membrane slab.
 * @param diel Vector container with all grids to be filled.
 * @author Gernot Kieseritzky
 */

 void calculateDielMap(vector<Vgrid*>& diel)
 {
  /* fill grid */
  for(vector<Vgrid*>::iterator itg=diel.begin(); itg!=diel.end(); itg++)
   Accessible::calculateDielMap(*itg);
  
  /* fill grid with membrane slab */
  Gridpoint pt(diel[0], 0.0, 0.0, zmin-diel[0]->zmin);
  fillSlab(pt, diel, zmin, zmax, epsm, epsp);
  /* return */
 }

 /**
 * @brief Fill kappa map including membrane slab.
 * @param diel Vector container with all grids to be filled.
 * @author Gernot Kieseritzky
 */

 void calculateKappaMap(double ionstr, vector<Vgrid*>& kappa) 
 {
  /* fill grid */
  for(vector<Vgrid*>::iterator itg=kappa.begin(); itg!=kappa.end(); itg++)
   Accessible::calculateKappaMap(ionstr, *itg);

  /* fill grid with membrane slab */
  Gridpoint pt(kappa[0], 0.0, 0.0, zmin-kappa[0]->zmin);
  fillSlab(pt, kappa, zmin, zmax, 0.0, 0.0);
  /* return */
 }

 /**
 * @brief Fill up membrane slab using Boundary Fill.
 * @param gp Grid point to start with.
 * @param grid Vector container with all grids to be filled.
 * @param zmin Minimum z coordinate of membrane slab.
 * @param zmax Maximum z coordinate of membrane slab.
 * @param fcol "Foreground color", i.e. dielectric constant / ion accessiblity of membrane.
 * @param bcol "Background color", i.e. dielectric constant / ion accessiblity of protein.
 * @author Gernot Kieseritzky
 */

 static void fillSlab(const Gridpoint& gp, const vector<Vgrid*>& grid, double zmin, double zmax, double fcol, double bcol);

/* class MembraneAccessible */
};

/**
 * @brief Controls calculation of SASA and ion accessibility maps.
 */

class MapController
{
protected:
 Accessible *acc;      /** accessibility object */
 Structure solute;     /** local copy of solute structure */

 APBSparms *apbsparms;
 vector<double>  gridCenterX; /** x coordinate of grid center */
 vector<double>  gridCenterY; /** y coordinate of grid center */
 vector<double>  gridCenterZ; /** z coordinate of grid center */
 vector<double*> dielx;       /** x-shifted dielectric maps */
 vector<double*> diely;       /** y-shifted dielectric maps */
 vector<double*> dielz;       /** z-shifted dielectric maps */
 vector<double*> kappa;       /** ion accessiblity maps */

 double ionStrength;  /** ionic strength */
 double ionMaxRad;    /** max ion radius */

 Vsurf_Meth smol;     /** type of surface model */

 string dx_filename;  /** filename stem for writing out DX files */

 typedef vector<double*>::iterator iter_map;

 /**
  * @brief Check whether the grid center has changed and update grid centers accordingly.
  * @param fstep Focusing step.
  * @param center coordinates of new grid center
  * @author Gernot Kieseritzky
  */

 bool gridCenterHasChanged(int fstep, double center[])
 {
  if( gridCenterX.at(fstep)!=center[0] || gridCenterY.at(fstep)!=center[1] || gridCenterZ.at(fstep)!=center[2] )
  {
   gridCenterX[fstep] = center[0];
	gridCenterY[fstep] = center[1];
	gridCenterZ[fstep] = center[2];
	return true;
  }

  return false;
 }

 /**
  * @brief Get accessibility map for focusing step.
  * @param fstep Focusing step.
  * @param center coordinates of new grid center
  * @param maptype Type of map being wrapped into Vgrid object
  * @author Gernot Kieseritzky
  */

 Vgrid* getMap(int fstep, double center[], Vdata_Type maptype)
 {
  int dime[3];
  double grid[3];
  double len[3];
  double min[3];
  double *data;
  
  // get grid size and resolution
  apbsparms->get_dime(fstep, dime);
  apbsparms->get_grid(fstep, grid);
  len[0] = grid[0]*(dime[0]-1);
  len[1] = grid[1]*(dime[1]-1);
  len[2] = grid[2]*(dime[2]-1);
  
  // calculate lowest corner of grid
  switch(maptype)
  {
   case VDT_DIELX: 
    min[0]  = center[0] - 0.5*len[0] + 0.5*grid[0];
    min[1]  = center[1] - 0.5*len[1];
    min[2]  = center[2] - 0.5*len[2];
	 data    = dielx.at(fstep);
    break;
   case VDT_DIELY: 
    min[0]  = center[0] - 0.5*len[0];
    min[1]  = center[1] - 0.5*len[1] + 0.5*grid[1];
    min[2]  = center[2] - 0.5*len[2];
	 data    = diely.at(fstep);
 	 break;
   case VDT_DIELZ:
    min[0]  = center[0] - 0.5*len[0];
    min[1]  = center[1] - 0.5*len[1];
    min[2]  = center[2] - 0.5*len[2] + 0.5*grid[2];
	 data    = dielz.at(fstep);
    break;
   case VDT_KAPPA:
    min[0]  = center[0] - 0.5*len[0];
    min[1]  = center[1] - 0.5*len[1];
    min[2]  = center[2] - 0.5*len[2];
	 data    = kappa.at(fstep);
    break;
	default:
	 return VNULL;
  }
  
  return Vgrid_ctor(dime[0], dime[1], dime[2], grid[0], grid[1], grid[2], min[0], min[1], min[2], data);
 }

 /**
  * @brief Refill maps if necessary.
  * @param pqr Solute structure.
  * @param centerx x coordinates of grid centers in each focusing step.
  * @param centery y coordinates of grid centers in each focusing step.
  * @param centerz z coordinates of grid centers in each focusing step.
  * @param dielX x-shifted dielectric maps
  * @param dielY y-shifted dielectric maps
  * @param dielZ z-shifted dielectric maps
  * @param kappa x-shifted dielectric maps
  * @param kappa kappa maps defining ion accessibility
  * @author Gernot Kieseritzky
  */

 virtual void fillMapsMol(const Structure& pqr, double centerx[], double centery[], double centerz[], Vgrid **dielX, Vgrid **dielY, Vgrid **dielZ, Vgrid **kappa);
 
 virtual void fillMapsSpline(const Structure& pqr, double centerx[], double centery[], double centerz[], Vgrid **dielX, Vgrid **dielY, Vgrid **dielZ, Vgrid **kappa)
 {
  fillMapsMol(pqr, centerx, centery, centerz, dielX, dielY, dielZ, kappa);
 }

 /**
  * @brief Deletes all maps in container.
  * @author Gernot Kieseritzky
  */

 void deleteAllMaps()
 {
  for(unsigned int i=0; i<dielx.size(); i++)
  {
   delete[] dielx[i];
   delete[] diely[i];
   delete[] dielz[i];
   delete[] kappa[i];
  }
 }

public:

 /**
  * @brief Constructor.
  * @param apbsparms APBSparms object with all APBS-related parameters
  * @author Gernot Kieseritzky
  */

 MapController(APBSparms *apbsparms);

 /**
  * @brief Get accessibility maps. Will refill maps or SASA if necessary.
  * @param pqr Solute structure.
  * @param centerx x coordinates of grid centers in each focusing step.
  * @param centery y coordinates of grid centers in each focusing step.
  * @param centerz z coordinates of grid centers in each focusing step.
  * @param dielX x-shifted dielectric maps
  * @param dielY y-shifted dielectric maps
  * @param dielZ z-shifted dielectric maps
  * @param kappa kappa maps defining ion accessibility
  * @author Gernot Kieseritzky
  */

 void fillMaps(const Structure& pqr, double centerx[], double centery[], double centerz[], Vgrid **dielX, Vgrid **dielY, Vgrid **dielZ, Vgrid **kappa)
 {
  if (smol==VSM_MOL || smol==VSM_MOLSMOOTH)
   fillMapsMol(pqr, centerx, centery, centerz, dielX, dielY, dielZ, kappa);
  else
   fillMapsSpline(pqr, centerx, centery, centerz, dielX, dielY, dielZ, kappa);
 }

 /**
  * @brief For debugging - write maps into files.
  * Is done in fillMapsMol(), here only the filname will be set.
  * @param filestem Filename stem for DX files.
  * @author Gernot Kieseritzky
  */

 void writeMaps(const string& filestem) {dx_filename = filestem;}

 /**
  * @brief Destructor.
  * @author Gernot Kieseritzky
  */

 ~MapController()
 {
  deleteAllMaps();
  if (acc!=NULL)
   delete acc;
 }
};

/**
 * @brief Alternative map controller including implicit membrane slab.
 */

class MembraneMapController : public MapController
{
 protected:
 double epsm ;      /** Membrane dielectric constant */
 double zmin, zmax; /** Dimensions of membrane slab */

 /**
  * @brief Refill maps if necessary.
  * @param pqr Solute structure.
  * @param centerx x coordinates of grid centers in each focusing step.
  * @param centery y coordinates of grid centers in each focusing step.
  * @param centerz z coordinates of grid centers in each focusing step.
  * @param dielX x-shifted dielectric maps
  * @param dielY y-shifted dielectric maps
  * @param dielZ z-shifted dielectric maps
  * @param kappa x-shifted dielectric maps
  * @param kappa kappa maps defining ion accessibility
  * @author Gernot Kieseritzky
  */

 virtual void fillMapsMol(const Structure& pqr, double centerx[], double centery[], double centerz[], Vgrid **dielX, Vgrid **dielY, Vgrid **dielZ, Vgrid **kappa);

 virtual void fillMapsSpline(const Structure& pqr, double centerx[], double centery[], double centerz[], Vgrid **dielX, Vgrid **dielY, Vgrid **dielZ, Vgrid **kappa)
 {
  fillMapsMol(pqr, centerx, centery, centerz, dielX, dielY, dielZ, kappa);
 }

 public:

 /**
  * @brief Constructor.
  * @author Gernot Kieseritzky
  */

 MembraneMapController(APBSparms *apbsparms, double epsm, double zmin, double zmax) : MapController(apbsparms)
 {
  this->epsm = epsm;
  this->zmin  = zmin;
  this->zmax  = zmax;
 }
};

/* APBS wrapper classes */

/**
 * @brief Base class of the APBS C++ wrapper.
 * APBS base class which does not do anything useful, yet.
 * Declares all necessary APBS structures and provides 
 * methods to initialize 'nosh' object.
 */

class APBS_base
{
protected:
 int debug;                        /** debug level flag */
 int ncalc;                        /** number of calculations to prepare */
 int setupflag;                    /** Has setup() been called already? */
 APBSparms *apbsparms;             /** APBS parameter object */

 Vmem *mem;                        /** APBS memory manager */
 Vcom *com;                        /** APBS I/O object */

 NOsh *nosh;                       /** APBS  parser object */
 
 Structure *pqr;                   /** PQR object */
 Valist *alist;                    /** APBS atom list */
 Vpmg **pmg;                       /** APBS multigrid objects */
 Vpmgp **pmgp;                     /** APBS multigrid parameter objects */
 Vpbe **pbe;                       /** APBS Poisson-Boltzmann objects */

 Vgrid **chargeMap;                /** APBS charge distribution maps (not used) */
 Vgrid **dielX, **dielY, **dielZ;  /** APBS dielectric accessbility functions */
 Vgrid **kappa;                    /** APBS ion accessibility functions */
 Vgrid **dummyMap;

 /* methods */

 /**
  * @brief Constructs APBS PBEparm object that contains 'parsed' PBE-related parameters.
  * Needed for NOsh_setupMG().
  * @author Gernot Kieseritzky
  * @return Pointer to new APBS PBEparm object.
  */

 PBEparm* fillPBEparm();

 /**
  * @brief Constructs APBS MGparm object that contains 'parsed' MG-related parameters.
  * Needed for NOsh_setupMG().
  * @author Gernot Kieseritzky
  * @return Pointer to new APBS MGparm object.
  */

 MGparm*  fillMGparm(MGparm_CalcType type);

 /**
  * @brief Constructs new NOsh object that contains all relevant parameters for a focusing calculation.
  * This is to mimic the parsing of a real APBS input file.
  * Flavor 1: Parameters set up that ion/solvent accessibility maps are going to be calculated.
  * @author Gernot Kieseritzky
  */

 void setupManualFocus(MGparm_CalcType type, double centerx[], double centery[], double centerz[]);

 /**
  * @brief Constructs new NOsh object that contains all relevant parameters for a focusing calculation.
  * This is to mimic the parsing of a real APBS input file.
  * Flavor 2: Parameters set up that ion/solvent accessibility maps are not going to be calculated
  *           but pre-calculated will be used.
  * @author Gernot Kieseritzky
  */

 void setupManualFocus(MGparm_CalcType type, double centerx[], double centery[], double centerz[], MapController *mapcontroller);

 /**
  * @brief Wraps various Vpmg data arrays into Vgrid maps.
  * (dielx, diely, dielz, kappa or u) It does not copy the content of the arrays!
  * @author Gernot Kieseritzky
  * @return Pointer to new APBS Vgrid object.
  */

 Vgrid*   exportMap(Vpmg *pmg, PBEparm *pbeparm, Vdata_Type maptype);

public:

 /**
  * @brief Standard constructor.
  * @author Gernot Kieseritzky
  */

 APBS_base(int debuglevel, Structure *pqr, APBSparms *apbsparms);

 /**
  * @brief Converts PQR into APBS atom list.
  * @author Gernot Kieseritzky
  * @return Pointer to new APBS atom list object.
  */

 static Valist* pqr2alist(Structure *pqr);

 /**
  * @brief Standard destructor.
  * @author Gernot Kieseritzky
  */

 virtual ~APBS_base();
};

/**
 * @brief Multigrid focusing calculation.
 * Performs a complete multigrid focusing calculation using APBS.
 */

class APBS_runner : public APBS_base
{
private:
 Potential potential;       /** Electrostatic potential on atom postions */

 /* methods */
 
 /**
  * @brief Constructs APBS core objects for focusing step i from parameter objects.
  * Essentially a replacement for initMG().
  * @author Gernot Kieseritzky
  * @return 1 - success, 0 - failed.
  */
 
 int initAPBScore(int i, MGparm *mgparm, PBEparm *pbeparm, double realCenter[]);

 /**
  * @brief Rescue electrostatic potential.
  * Potential is interpolated on off-grid positions and stored in 'potential'.
  * @author Gernot Kieseritzky
  */

 void rescue_potential(Vgrid *pot);

public:

 /* methods */

 /**
  * @brief Standard constructor.
  * @author Gernot Kieseritzky
  */

 APBS_runner(int debuglevel, Structure *pqr, APBSparms *apbsparms);

 /**
  * @brief Defines calculation type: Calculate maps and perform real computation.
  * Simply calls setupManualFocus(MCT_MANUAL, ...) and passes grid centers.
  * @author Gernot Kieseritzky
  */

 void setup(double centerx[], double centery[], double centerz[]);

 /**
  * @brief Defines calculation type: perform real computation with pre-calculated maps.
  * Simply calls setupManualFocus(MCT_MANUAL, ...) and passes grid centers and pre-calculated maps.
  * @author Gernot Kieseritzky
  */

 void setup(double centerx[], double centery[], double centerz[], MapController *mapcontroller);

 /**
  * @brief Defines calculation type: Calculate only maps and perform no computation.
  * Simply calls setupManualFocus(MCT_DUM, ...) and passes grid centers.
  * Good for testing.
  * @author Gernot Kieseritzky
  */

 void setupDummy(double centerx[], double centery[], double centerz[]);

 /**
  * @brief Defines calculation type: Perform no computation and use pre-calculated maps.
  * Simply calls setupManualFocus(MCT_DUM, ...) and passes grid centers and maps.
  * Good for testing.
  * @author Gernot Kieseritzky
  */

 void setupDummy(double centerx[], double centery[], double centerz[], MapController *mapcontroller);

 /**
  * @brief Starts APBS multigrid focusing calculation.
  * Can only be called after setup()!
  * @author Gernot Kieseritzky
  */

 void run();

 /**
  * @brief Get eletrostatic potential.
  * @author Gernot Kieseritzky
  * @return reference to electrostatic potential object.
  */

 Potential& get_potential() {return potential;}

 /**
  * @brief Standard destructor.
  * @author Gernot Kieseritzky
  */
  
 virtual ~APBS_runner() {/* do nothing */};
};

#endif /* !APBSWRAPPER_H */
