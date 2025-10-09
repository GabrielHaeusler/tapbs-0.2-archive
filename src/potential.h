#include <vector>
#include <map>
#include "titratable.h"

/** 
 * @brief Container for an electrostatic potential function.
 * Class defines an electrostatic potential.
 * It provides a container and interface
 * functions to access or traverse the values.
 */

class Potential
{
protected:
 typedef map<Vector3D, double> Container;
 typedef Container::iterator iter_pot;
 typedef Container::const_iterator const_iter_pot;
 
 Container values;            /** stores potential values */
 iter_pot potIterator;        /** iterator over potential values */

public:

 /** @brief Struct for external access of values.
  *  @author Gernot Kieseritzky
  */

 struct extAccess
 {
  double x;
  double y;
  double z;
  double pot;
 };

 /** @brief Standard constructor.
  *  @author Gernot Kieseritzky
  */

 Potential() {potIterator = values.begin();}

 /** @brief Alternative constructor intitializing from arrays.
  *  @author Gernot Kieseritzky
  */

 Potential(int size, double x[], double y[], double z[], double pot[])
 {
  copy_from_arrays(size, x, y, z, pot);
  potIterator = values.begin();
 }

 /** @brief Copy constructor
  *  @author Gernot Kieseritzky
  */

 Potential(const Potential& pot);

 /** @brief Interface to iterator from outside.
  *  @author Gernot Kieseritzky
  *  @returns True if iterator has not reached the end of container.
  */

 bool hasNext();

 /** @brief Interface to iterator from outside.
  *  @author Gernot Kieseritzky
  *  @returns struct with (x,y,z) and electrostatic potential at this point.
  */

 struct extAccess getNext();

 /** @brief takes a point in space and returns the corresponding el. potential.
  *  If the point is not in the container it will throw an exception.
  *  @author Gernot Kieseritzky
  *  @returns el. potential (double value)
  */

 double get_value_at(const Vector3D& v) const;
 double get_value_at(double x, double y, double z) const;

 /** @brief Puts point & coresp. potential into container
  *  @author Gernot Kieseritzky
  */

 void set_value_at(double x, double y, double z, double val);

 /** @brief Copy coordinates and potentials into arrays.
  *  @author Gernot Kieseritzky
  */

 void copy_to_arrays(double x[], double y[], double z[], double pot[]) const;

 /** @brief Copy coordinates and potentials from arrays.
  *  @author Gernot Kieseritzky
  */

 void copy_from_arrays(int size, double x[], double y[], double z[], double pot[]);

 /** @brief How many points are stored?
  *  @author Gernot Kieseritzky
  *  @returns Size of container.
  */

 int nrOfPoints() const {return values.size();}

};

/**
 * @brief Stores electrostatic potential functions.
 * Container for electrostatic potentials of every charge state 
 * of the protein and model compounds.
 */

class PotentialsOfResidue
{
protected:
 typedef vector<Potential*> Container; /** Container type used for storage */
 typedef Container::iterator iter_pot; /** Container specific iterator */
 
 Container protein;  /** protein potentials for each charge state of residue */
 Container model;    /** model potentials for each charge state of residue */

public:

 /** @brief Standard constructor.
  *  @author Gernot Kieseritzky
  */
  
 PotentialsOfResidue() {/* do nothing */}

 int nrOfProteinStates() {return protein.size();}
 int nrOfModelStates() {return model.size();}

 /** @brief Add protein potential to containers.
  *  @author Gernot Kieseritzky
  */

 void addProteinStatePotential(const Potential& pot) {protein.push_back( new Potential(pot) );}

 /** @brief Add model potential to containers.
  *  @author Gernot Kieseritzky
  */

 void addModelStatePotential(const Potential& pot) {model.push_back( new Potential(pot) );}

 /** @brief get protein potential
  *  @throw UnknownStateException if state does not exist
  *  @author Gernot Kieseritzky
  */

 const Potential& getProteinStatePotential(int state);

 /** @brief get model potential
  *  @throw UnknownStateException if state does not exist
  *  @author Gernot Kieseritzky
  */

 const Potential& getModelStatePotential(int state);

 /** @brief Write potentials into file.
  *  @author Gernot Kieseritzky
  */

 void writeToFile(const string& filename);

 /** @brief Read potentials from file.
  *  @author Gernot Kieseritzky
  *  @returns 1 - if file was found, 0 - if no file was found
  */

 int readFromFile(const string& filename);

 /** @brief Calculate Born energy.
  *  Needs state number and Structures with corresponding charge sets.
  *  Formula (R: reference state, S: state):
  \f[
       \Delta\Delta G_{born} = \sum^{N_{\mu}}_{i=1} Q^{S}_{i,\mu}
               \left[ \phi_{p}\left( r_i; Q^{S}_{i,\mu}\right) - \phi_{m}\left( r_i; Q^{S}_{i,\mu}\right) \right]
          -\sum^{N_{\mu}}_{i=1} Q^{R}_{i,\mu}
               \left[ \phi_{p}\left(r_i; Q^{R}_{i,\mu}\right) - \phi_{m}\left(r_i; Q^{R}_{i,\mu}\right) \right]
  \f]
  *  @author Gernot Kieseritzky
  */

 double calcBorn(int state, Structure *mstat, Structure *mref);

 /** @brief Calculate Back energy.
  *  Needs state number and Structures with background charges (protein/model).
  *  Formula (R: reference state, S: state):
  \f[
  \Delta\Delta G_{back} = \sum^{N_{p}}_{i=1} q_i
                \left[ \phi_{p}\left(r_i; Q^{S}_{i,\mu}\right) - \phi_{p}\left(r_i; Q^{R}_{i,\mu}\right) \right]
                         -\sum^{N_{m}}_{i=1} q_i
                \left[ \phi_{m}\left(r_i; Q^{S}_{i,\mu}\right) - \phi_{m}\left(r_i; Q^{R}_{i,\mu}\right) \right]
  \f]
  *  @author Gernot Kieseritzky
  */

 double calcBack(int state, Structure *pback, Structure *mback);

 /** @brief Calculate Site Site iteraction term for W-Matrix of a state.
  *  Needs state number of residue 1 whose potentials are used 
  *  and structures of residue 2 with reference and state charges.
  * Formula:
  \f[
    W_{\mu\nu} = \sum^{N_{\mu}}_{i=1}  Q^{S}_{i,\mu} \left[ \phi_{p}\left(r_{Q^{S}_{i,\mu}}, Q^{S}_{i,\nu}\right) - \phi_{p}\left(r_{Q^{S}_{i,\mu}}, Q^{R}_{i,\nu}\right) \right]
                                     - Q^{R}_{i,\mu} \left[ \phi_{p}\left(r_{Q^{R}_{i,\mu}}, Q^{S}_{i,\nu}\right) - \phi_{p}\left(r_{Q^{R}_{i,\mu}}, Q^{R}_{i,\nu}\right) \right]
  \f]
  *  @author Gernot Kieseritzky
  *  @returns W matrix entry Wmv (double).
  */

 double calcWmv(int state, int scale, Structure *mstat, Structure *mref);

 /** @brief Frees all memory associated with potentials
  *  @author Gernot Kieseritzky
  */

 void deletePotentials();

 /** @brief Standard destructor.
  *  @author Gernot Kieseritzky
  */

 ~PotentialsOfResidue() { deletePotentials(); }
};
