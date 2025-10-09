#include "PQR.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

using namespace std;

/**
 * @brief A possible state of a titratable residue.
 * That state can either be a set of charges or a conformer.
 * @author Gernot Kieseritzky
 */

class State
{
protected:

 typedef map<string, PQR_atom*> Container;
 typedef Container::iterator iter_atom;
 
 string id;                    /** transition type identifier */
 double energy;                /** non-electrostatic energy shift */

 Container atoms;              /** atomic states */

/** @brief  Deletes all atoms in container.
 *  @author Gernot Kieseritzky
 */

 void delete_atoms();

 void copyMembers(const State& that);

public:

/** @brief   Constructor
 *  @author  Gernot Kieseritzky
 *  @param id Transition type identifier, e.g. 'D'
 *  @param energy Energy relative to reference state.
 */

 State(string id, double energy) {this->id = id; this->energy = energy;}

/** @brief   Copy constructor
 *  @author  Gernot Kieseritzky
 */

 State(const State& that) {copyMembers(that);}

/** @brief   Get transition type identifier ('R', 'D', ...)
 *  @author  Gernot Kieseritzky
 *  @returns String containing identifier.  
 */

 const string& getID() {return id;}

/** @brief   Get non-electrostatic energy.
 *  @author  Gernot Kieseritzky
 *  @returns Double value representing energy.
 */

 double getEnergy() {return energy;}

/** @brief   Set transition type.
 *  @author  Gernot Kieseritzky
 */

 void setID(const string& val) {id = val;}

/** @brief   Set energy.
 *  @author  Gernot Kieseritzky
 */

 void setEnergy(double val) {energy = val;}

/** @brief   Is this state anionic?
 *  @author  Gernot Kieseritzky
 *  @returns True if total charge is negative.
 */

 bool isAnionic();

/** @brief   Add atom to state.
 *  @author  Gernot Kieseritzky
 *  @param atom Atom to be added.
 *  @returns 0 if atom was present already.
 */

 int add_atom(const PQR_atom& atom);

/** @brief   Fill missing information for every atom.
 *  adds atom number, resid, radius, chain, altloc, segname and coordinates if missing.
 *  @author  Gernot Kieseritzky
 *  @param pqr Source structure to complete atoms with missing information.
 */

 void fillMissing(Structure *pqr);

/** @brief Check if atom number exists in PQR structure.
 *  All atoms in state should have a corresponding atom in PQR structure.
 *  Therefore we check whether atom number exists in PQR structure.
 *  During parsing atom number was set to -1 (undef) and should have been
 *  replaced by existing atom number in fillMissing.
 *  @throw AtomNotFoundException
 *  @author  Gernot Kieseritzky
 *  @param pqr Source structure.
 */

 void check(Structure *pqr);

/** @brief   Change charges and coordinates according to state description.
 *  @author  Gernot Kieseritzky
 *  @param pqr Structure to be changed. Should have non-empty atom selection.
 */

 void change(Structure *pqr);

/** @brief   Puts charges on atom positions.
 *  @author  Gernot Kieseritzky
 *  @param pqr Structure to be changed.
 */

 void charge(Structure *pqr);

/** @brief   Removes charges on atom positions.
 *  @author  Gernot Kieseritzky
 *  @param pqr Structure to be changed.
 */

 void uncharge(Structure *pqr);

/** @brief   Select atom positions associated with this state.
 *  @author  Gernot Kieseritzky
 *  @param pqr Structure to be changed.
 */

 void selectCoordinates(Structure *pqr);

/** @brief  Print state.
 *  @param  os Output stream.  
 *  @author Gernot Kieseritzky
 */

 void print(ostream& os);

/** @brief   Standard destructor.
 *  @author  Gernot Kieseritzky
 */
 
 ~State() {delete_atoms();}

/** @brief   Assignment operator.
 *  @author  Gernot Kieseritzky
 */

 State& operator=(const State& that)
 {
  if (this != &that) 
  {
   delete_atoms();
	copyMembers(that);
  }
  return *this;
 }
};

/**
 * @brief Description of a titratable residue.
 * @author Gernot Kieseritzky
 */

class Titratable
{
protected:
 string name;            /** name of titratable residue, e.g. 'HSP-19' */
 int    resid;           /** resid of titratable residue */
 vector<State*> states;  /** states that residue can adopt */

 typedef vector<State*>::iterator iter_states;

 /**
  * @brief Define reference state.
  * Tries to find 'R'-state and moves it to position 0.
  * If no 'R' was defined state 0 is assumed to be reference state.
  * @author Gernot Kieseritzky
  * @returns Unsigned integer that is the index of the reference state
  */

 void defineReferenceState();

 /**
  * @brief Remove all states.
  * @author Gernot Kieseritzky
  */

 void delete_states();

 void copyMembers(const Titratable& that);

public:

 /**
  * @brief Constructor.
  * @param resid Residue number.
  * @param residue Residue name, e.g. 'GLU'
  * @author Gernot Kieseritzky
  */

 Titratable(int resid, const string& residue);

 /**
  * @brief Copy constructor.
  * @param that titratable residue to be copied
  * @author Gernot Kieseritzky
  */

 Titratable(const Titratable& that) {copyMembers(that);}

 /**
  * @brief Get number of states for this residue.
  * @author Gernot Kieseritzky
  * @returns Number of states.
  */

 unsigned int nrOfStates() {return states.size();}

 /**
  * @brief Get name of titratable residue.
  * @author Gernot Kieseritzky
  * @returns String reference containing residue name, e.g. 'GLU-15'
  */

 const string& getName() {return this->name;}

 /**
  * @brief Get residue number of this residue
  * @author Gernot Kieseritzky
  * @returns Integer representing the residue number
  */

 int getResID() {return this->resid;}

 /**
  * @brief Get energy of a state.
  * @author Gernot Kieseritzky
  * @param state Number of state.
  * @returns Double value representing the energy of state.
  */

 double getEnergy(int state);

 /**
  * @brief Get transition type identifier.
  * @author Gernot Kieseritzky
  * @param state Number of state.
  * @returns Constant string reference containing transition type identifier.
  */

 const string& getID(int state);

 /**
  * @brief Parse ST-file.
  * @param filename String containing filename.
  * @author Gernot Kieseritzky
  */

 void parseST(const string& filename);

 /**
  * @brief Fill missing information in states from PQR structure.
  * @param pqr Pointer to pqr structure.
  * @author Gernot Kieseritzky
  */

 void fillMissing(Structure *pqr);

 /**
  * @brief Select me in PQR structure.
  * @author Gernot Kieseritzky
  */

 virtual void select(Structure *pqr)
 {
  //pqr->select_clear();
  pqr->select_resid(resid);
 }

 /**
  * @brief Select coordinates from structure associated with state.
  * @param state State number.
  * @param pqr PQR structure.
  * @throw UnknownStateException
  * @author Gernot Kieseritzky
  */

 void selectCoordinatesOfState(int state, Structure *pqr);

 /**
  * @brief Switch titratable residue to state 'state' in Structure 'pqr'.
  * @param state State number.
  * @param pqr PQR structure.
  * @author Gernot Kieseritzky
  */

 void switchToState(int state, Structure *pqr);

 /**
  * @brief Switch titratable residue to reference state.
  * @param pqr PQR structure.
  * @author Gernot Kieseritzky
  */

 void switchToReferenceState(Structure *pqr) {switchToState(0, pqr);}

 /**
  * @brief Instead of switching it only charges corresponding atom positions.
  * @param state State number.
  * @param pqr PQR structure.
  * @author Gernot Kieseritzky
  */

 void setChargesOfState(int state, Structure *pqr);

 /**
  * @brief Remove charges again.
  * Usually applied after setChargesOfState
  * @param state State number.
  * @param pqr PQR structure.
  * @author Gernot Kieseritzky
  */

 void clearChargesOfState(int state, Structure *pqr);

 /**
  * @brief Print states
  * @param os Output stream.
  * @author Gernot Kieseritzky
  */

 void printStates(ostream& os);

 /**
  * @brief Is state anionic?
  * Used only in legacy format calculations.
  * @author Gernot Kieseritzky
  * @param state Number of state.
  * @returns true if anionic.
  */

 bool isAnionic(int state);

 /**
  * @brief Standard destructor.
  * @author Gernot Kieseritzky
  */

 virtual ~Titratable() {delete_states();}

 Titratable& operator=(const Titratable& that)
 {
  if (this != &that) 
  {
   delete_states();
	copyMembers(that);
  }
  return *this;
 }
};

/**
 * @brief Description of a titratable residue that knows about multi chain proteins.
 * @author Gernot Kieseritzky
 */

class McTitratable : public Titratable
{
protected:
 string segname;

public:

 /**
  * @brief Constructor.
  * @param resid Residue number.
  * @param residue Residue name, e.g. 'GLU'
  * @param segname Segment name, e.g. 'PROTEIN'
  * @author Gernot Kieseritzky
  */

 McTitratable(int resid, const string& residue, const string& segname);

 /**
  * @brief Select me in PQR structure.
  * @author Gernot Kieseritzky
  */

 virtual void select(Structure *pqr)
 {
  //pqr->select_clear();
  pqr->select_resid(segname, resid);
 }

};
