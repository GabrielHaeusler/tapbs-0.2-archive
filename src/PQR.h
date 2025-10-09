#include "Cartesian.h"

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>

using namespace std;

/** @brief Class describing atoms in PQR format.
 *  @author  Gernot Kieseritzky
 */

class PQR_atom
{
protected:

 int atomno;                /** running atom number */
 string atomid;             /** atom string identifier */
 char altloc;               /** conformer identifier */
 string residue;            /** residue name */
 char chain;                /** chain ID */
 int resid;                 /** residue number */
 Vector3D cartesian;        /** cartesian coordinate */
 double charge;             /** atomic partial charge */
 double radius;             /** atomic VDW radius */
 string segname;            /** segment ID */

public:
 static int instanceCounter;/** counts instantiations */

/** @brief  Standard constructor
 *  @author  Gernot Kieseritzky
 *  @ingroup PQR_atom
 */

 PQR_atom() 
 {
  this->atomno = 1;
  this->atomid = " ";
  this->altloc = ' ';
  this->residue= " ";
  this->chain  = ' ';
  this->resid  = 1;
  this->charge = 0.0;
  this->radius = 0.0;
  this->segname= " ";
  instanceCounter++;
 }

/** @brief  Constructor: creates new PQR_atom from PQR line string.
 *  @author  Gernot Kieseritzky
 *  @ingroup PQR_atom
 */
 
 PQR_atom(const string& pqrline);

/** @brief  Constructor: creates new PQR_atom manually from parameters.
 *  @author  Gernot Kieseritzky
 *  @ingroup PQR_atom
 */
 
 PQR_atom(int atomno, int resid, double charge, double radius, const string& residue, const string& atomid);

/** @brief  Parses PQR line provided by string.
 *  @author  Gernot Kieseritzky
 *  @return 1 - success, 0 - parse error
 *  @ingroup PQR_atom
 */
 
 int fromPQRLine(const string& pqrline);

 /* interface functions to access atom-specific parameters */
 
 void set_atomno(int atomno)             {this->atomno = atomno;}
 void set_atomid(const string& atomid)   {this->atomid = atomid;}
 void set_conf(char altloc)              {this->altloc = altloc;}
 void set_residue(const string& residue) {this->residue = residue;}
 void set_chain(char id)                 {this->chain  = id;}
 void set_resid(int resid)               {this->resid = resid;}
 void set_charge(double charge)          {this->charge = charge;} 
 void set_radius(double radius)          {this->radius = radius;}
 void set_segment(const string& segname) {this->segname = segname;}

 Vector3D& coordinates()     {return this->cartesian;}

 int get_atomno() const            {return this->atomno;}
 const string& get_atomid() const  {return this->atomid;}
 char get_conf() const             {return this->altloc;}
 const string& get_residue() const {return this->residue;}
 char get_chain() const            {return this->chain;}
 int get_resid() const             {return this->resid;}
 double get_charge() const         {return this->charge;}
 double get_radius() const         {return this->radius;}
 const string& get_segment() const {return this->segname;}
 
/** @brief injects PQR_atom into output stream in PQR format
 *  @author  Gernot Kieseritzky
 *  @ingroup PQR_atom
 */
 
 static void PQR(ostream& os, PQR_atom& atom);
};

/** @brief  Overload output stream operator.
 *  @author  Gernot Kieseritzky
 *  @ingroup PQR_atom
 */

ostream& operator<<(ostream& os, PQR_atom& atom);

/** @brief  Satom extends PQR_atom by solute dielectric constant.
 *  Used for calculations with more than 2 dielectric constants.
 *  @author  Gernot Kieseritzky
 */

class Satom : public PQR_atom
{
protected:
 double soluteDiel;
 
public:

/** @brief  Standard constructor
 *  @author  Gernot Kieseritzky
 *  @ingroup Satom
 */

 Satom() : PQR_atom() {soluteDiel = 0.0;}

/** @brief  Manual constructor
 *  @author  Gernot Kieseritzky
 *  @ingroup Satom
 */
 
 Satom(int atomno, int resid, double charge, double radius, double soluteDiel, const string& residue, const string& atomid);

 /* interface functions */
 
 void set_soluteDiel(double soluteDiel) {this->soluteDiel = soluteDiel;}
 double get_soluteDiel() {return soluteDiel;}
};

/** @brief Structure combines single atoms to molecule structures.
 *  Provides basic functionality for structure manipulations.
 *  @author Gernot Kieseritzky
 */

class Structure
{
protected:
 typedef vector<Satom*>                            AtomContainer;

 typedef multimap<int, Satom*>                     AtomsByResid; 
 typedef map<Vector3D, Satom*>                     AtomsByVec;
 typedef map<string, AtomsByResid>                 AtomsBySegname;

 typedef AtomContainer::iterator                   iter_atom;
 typedef AtomContainer::const_iterator             const_iter_atom;
 typedef AtomsByResid::iterator                    iter_resid;
 typedef AtomsByResid::const_iterator              const_iter_resid;
 typedef AtomsByVec::iterator                      iter_coord;
 typedef AtomsByVec::const_iterator                const_iter_coord;
 typedef AtomsBySegname::iterator                  iter_segname;
 typedef AtomsBySegname::const_iterator            const_iter_segname;

 typedef AtomsByVec::value_type                    pvec;
 typedef AtomsByResid::value_type                  presid;
 typedef AtomsBySegname::value_type                psegname;

 /* containers */

 AtomContainer      atoms; 			/** atom storage */
 AtomContainer      atomSelection;	/** selected atoms storage */

 AtomsByVec         atomsByCoord;	/** for atom access by coordinates */
 AtomsBySegname     atomsBySegname;	/** for atom access by segment name */

 AtomContainer::iterator atomIterator; /** used for hasNext/getNext() */
 AtomContainer::iterator seleIterator; /** used for hasNext_in_selection/getNextin_selection() */

 /* properties */

 Vector3D center;                   /** geometric center   */
 Vector3D min;                      /** minimum coordinate */
 Vector3D max;                      /** maximum coordinate */
 Vector3D smin;                     /** minimum solute coordinate */
 Vector3D smax;                     /** maximum solute coordinate */
 int residues;                      /** number of residues **/

 /* methods */

/** @brief  Refill container structures.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 
 
 void refill_containers();

/** @brief   Copy atoms.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */

 void copyMembers(const Structure& that);

public:

/** @brief  Standard constructor
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */

 Structure() {atomIterator=atoms.begin(); seleIterator=atomSelection.begin(); residues=0;}

/** @brief  Copy constructor
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */

 Structure(const Structure& that);

/** @brief  Opens and reads PQR file
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */

 void fromPQR(const string& pqrfile);

/** @brief  Builds new Structure from atom selection
 *  @author  Gernot Kieseritzky
 *  @returns Pointer to new Structure object.
 *  @ingroup Structure
 */

 Structure* fromSelection();

/** @brief  Constructs a MEAD-like model compound from selection.
 *  @author  Gernot Kieseritzky
 *  @returns Pointer to new Structure object.
 *  @ingroup Structure
 */

 Structure* createModelFromSelection();

/** @brief  Number of atoms stored
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 *  @returns Number of atoms
 */
 
 int nrOfAtoms() const {return atoms.size();}
 
/** @brief  Number of residues stored
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 *  @returns Number of residues.
 */ 
 
 int nrOfResids() const 
 {
  return residues;
 }

 /* interface functions for low level atom manipulation */
 
/** @brief  Access to single atom
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 *  @returns Reference to Satom in atom container
 */ 
 
 Satom& get_atom(int atomno) const;

/** @brief  Access to a single atom by coordinates
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 *  @returns Satom
 */ 
 
 Satom& get_atom(const Vector3D& v) const;
 Satom& get_atom(double x, double y, double z) const {Vector3D v(x, y, z); return get_atom(v);}

/** @brief  Adds a single atom
 *  Don't forget to call hash_atoms() & calc_geometry() afterwards!
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void add_atom(Satom* atom);

/** @brief  Delete atoms in selection
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void delete_atoms_in_selection();

/** @brief  Delete all atoms in structure
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void delete_all_atoms();
 
/** @brief  Renumber and sort all atoms.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void renumber();

/** @brief  Calculate center of geometry and min/max coordinates
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void calc_geometry();
 
/** @brief  Update all data structures after manipulations.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void update() 
 {
  refill_containers();
  renumber();
  calc_geometry();
 }

/** @brief  Do coordinates exist?
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 bool hasPoint(const Vector3D& v) {return atomsByCoord.find(v)!=atomsByCoord.end();}

/** @brief  Append another structure.
 *  For performance reasons it does not call update() automatically.
 *  @param other Structure that is to be appended.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void append(const Structure& other);

/** @brief  Like 'append' but takes only unknown atom positions.
 *  For performance reasons it does not call update() automatically.
 *  @param other Structure that is to be appended.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */

 void appendAtUnknownCoordinates(const Structure& other);

/** @brief  Set all charges to zero
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void clear_all_charges();

/** @brief  Set all charges in selection to zero
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void clear_charges_in_selection();

/** @brief  Set solute dielectric constant in selection.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 
 
 void set_soluteDiel_in_selection(double epsp);
 
 /* Selection operators */

/** @brief  Select all atoms
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void select_all();

/** @brief  Selects atoms by residue number.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 *  @returns Number of atoms selected.
 */ 

 int select_resid(int resid);

/** @brief  Selects atoms by segment name and residue number.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 *  @returns Number of atoms selected.
 */ 

 int select_resid(const string& segname, int resid);

/** @brief  Selects atom by its atom number
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void select_atomno(int atomno);

/** @brief  Selects atom by its coordinate
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void select_coordinate(const Vector3D& v);

/** @brief  Clears selection.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void select_clear();

/** @brief Is there an atom we haven't touched yet?
 *  hasNext/getNext is an JAVA-like iterator interface
 *  to atoms in the structure useful for external access.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 bool hasNext();
 bool hasNext_in_selection();

/** @brief Get reference to next atom.
 *  @author  Gernot Kieseritzky
 *  @returns reference to Satom
 *  @ingroup Structure
 */ 
 
 Satom& getNext();
 Satom& getNext_in_selection();

 /* coordinate manipulations */

/** @brief  Returns geometry center
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 const Vector3D& get_center() const {return center;}
 void get_center(double *x, double *y, double *z) const;

/** @brief  Returns minimum coordinate
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 const Vector3D& getMin() const {return min;}

/** @brief  Returns maximum coordinate
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 const Vector3D& getMax() const {return max;}

/** @brief  Returns min and max coordinates.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void get_minmax(double min[3], double max[3]) const;

/** @brief  Returns the solute size.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 Vector3D getSoluteLength() const { return (smax-smin); }

/** @brief  Rotate structure around 'axis' by 'deg' degrees
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void rotate(Vector3D& axis, double deg);

/** @brief  Translate structure along 'axis'
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void translate(Vector3D& axis);

 /* output streams */

/** @brief  Inject structure into output stream.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void print_pqr(ostream& os);

/** @brief  Inject current selection into output stream.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void print_pqr_selection(ostream& os);

/** @brief  Write structure into PQR file.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 void writeToFile_pqr(const string& filename);

/** @brief  Destructor. Frees memory from atoms.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 ~Structure() {delete_all_atoms();}

/** @brief  Assignment operator.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 Structure& operator=(const Structure& that)
 {
  if (this!=&that)
  {
   delete_all_atoms();
	copyMembers(that);
	update();
  }
  return *this;
 }

/** @brief   Equality operator.
 *  Two Structures are equal iff all coordinates & radii are identical
 *  and have the same number of atoms.
 *  We ignore charges because this would not change the SASA. 
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */ 

 bool operator==(const Structure& that) const
 {
  AtomContainer::const_iterator iterThis, iterThat;
 
  /* test number of atoms */
  if(this->nrOfAtoms() != that.nrOfAtoms())
   return false;
  
  /* test coordinates */
  for(iterThis=atoms.begin(), iterThat=that.atoms.begin(); (iterThis!=atoms.end())&&(iterThat!=that.atoms.end()); iterThis++, iterThat++)
  {
   Satom *thisAtom = *iterThis;
	Satom *thatAtom = *iterThat;
   if ( ( thisAtom->coordinates() != thatAtom->coordinates() )
        || ( thisAtom->get_radius() != thatAtom->get_radius() ) )
	 return false;
  }
  
  return true;
 }

/** @brief   Inequality operator.
 *  @author  Gernot Kieseritzky
 *  @ingroup Structure
 */

 bool operator!=(const Structure& that) const
 {
  return !(*this==that);
 }
};
