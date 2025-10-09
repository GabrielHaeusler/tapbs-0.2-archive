#include "apbswrapper.h"

#include <time.h>
#include <sys/times.h>

#define T apbsparms->get_temperature()

using namespace std;

/** @brief Class interface to titration calculations.
 *  @author Gernot Kieseritzky
 */

class Titration
{
protected:

 typedef vector<Titratable*> Titratables;                /** Container type used for titratable residues */
 typedef Titratables::iterator iter_titratables;         /** Iterator type of titr. residue container */
 
 APBSparms *apbsparms;                                   /** all APBS related parameters */
 Structure *pqr;                                         /** PQR containing all partial charges */
 Titratables residues;                                   /** vector of titratable residues */

 string pqrfile;                                         /** path to PQR file */
 string sitesfile;                                       /** path to SITES file */
 string jobname;                                         /** filename stem for output files */

 string stateFileExtension;                              /** '.st'    */
 string potenFileExtension;                              /** '.potat' */
 string pkintFileExtension;                              /** '.pkint' */
 string wmatrFileExtension;                              /** '.g'     */

 bool dummy;                                             /** If set no potentials will be calculated 
                                                             but everything is setup in the usual way.
																				 Good for testing the parameters chosen.
																	        */

 bool symcheck;                                          /** perform symmetry check? */
 bool writemap;                                          /** write map? */

 bool membrane;                                          /** include implicit membrane slab? */
 double zmin;                                            /** membrane slab minimum z-coordinate */
 double zmax;                                            /** membrane slab maximum z-coordinate */
 double epsm;                                            /** membrane dielectric constant */

 clock_t starttime;                                      /** Time when object was instantiated */
 clock_t stoptime;                                       /** Time when object is destroyed     */
 
 void startTimer()
 {
  tms start;
  starttime = times(&start);
 }
 
 string stopTimer()
 {
  ostringstream ostr;
  tms stop;
  stoptime = times(&stop);
  ostr << "Titration completed and used the follwing amount of time:" << endl;
  ostr << "real: " << static_cast<double>(stoptime - starttime)/static_cast<double>(sysconf(_SC_CLK_TCK)) << "s"
       << "  cpu-user: " << static_cast<double>(stop.tms_utime)/static_cast<double>(sysconf(_SC_CLK_TCK)) << "s"
       << "  cpu-sys: " << static_cast<double>(stop.tms_stime)/static_cast<double>(sysconf(_SC_CLK_TCK)) << "s" << endl;

  return ostr.str();
 }

/** @brief Deletes all Titratable residues in 'residues'.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */
 
 void deleteTitratableResidues();

/** @brief Parses sites-File to get a list 
           of titratable residues 
			  (-> sets up 'residues')
 *  @author Gernot Kieseritzky
 *  @param filename Filename of SITES file
 *  @ingroup Titration
 */
 
 void parseSITES(const string& filename);

/** @brief Setup grid centers according to focusing scheme.
 *         Result is stored in 'center'-Array and must have 
 *         the following structure:
 *         center[0][i] = x coordinate of grid center in focusing step i
 *			  center[1][i] = y coordinate of grid center in focusing step i
 *         center[1][i] = z coordinate of grid center in focusing step i
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 void setGridCenters(APBSparms *apbsparms, double *center[3], Structure *p, Structure *m);

/** @brief Set all solute dielectric constants according to user setup.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */
 
 void setSoluteDielectricConstants(InputParser& parser);

/** @brief Rotates molecule to align with membrane axis.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */
 
 bool setupMembrane(InputParser& parser);

/** @brief Create protein dielectric volume.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 *  @returns Pointer to PQR structure representing the combined volume.
 */

 virtual Structure* createProteinDielectricVolume() = 0;

/** @brief Run APBS to calculate protein electrostatic potentials.
 *  @author Gernot Kieseritzky
 *  @param proteinVolume Protein dielectric volume.
 *  @param filestem Additional filename suffix used for potential files, e.g. 'pairs'
 *  @ingroup Titration
 */

 virtual void calc_protein_potentials(Structure *proteinVolume, const string& filestem = "") = 0;

/** @brief Run APBS to calculate model compound electrostatic potentials.
 *  @author Gernot Kieseritzky
 *  @param proteinVolume Protein dielectric volume.
 *  @param filestem Additional filename suffix used for potential files, e.g. 'pairs'
 *  @ingroup Titration
 */

 virtual void calc_model_potentials(Structure *proteinVolume, const string& filestem = "") = 0;

/** @brief Calculate intrinsic pKas.
 *  @author Gernot Kieseritzky
 *  @param filestem Additional filename suffix used for potential files, e.g. 'pairs'
 *  @ingroup Titration
 */

 virtual void calc_pka(const string& filestem = "") = 0;


/** @brief Creates & calculates W matrix.
 *  @author Gernot Kieseritzky
 *  @param filestem Additional filename suffix used for potential files, e.g. 'pairs'
 *  @ingroup Titration
 */
 
 virtual void calc_wmatrix(const string& filestem = "") = 0;
 
/** @brief Symmetrizes the calculated W matrix.
 *  If W matrix is not symmetric it will calculate
 *  the arithmetic mean of both halves and replace
 *  the elements.
 *  Must be called after calc_wmatrix()!
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 virtual void symmetrize_wmatrix() = 0;

/** @brief Write current W matrix into file.
 *  Must be called after calc_wmatrix()!
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 virtual void write_wmatrix() = 0;

public:

 /* methods */

/** @brief Standard constructor, needs Parser object for parameters that control titration.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 Titration(InputParser& parser);

/** @brief Extends parser by keywords needed by class Titration.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 static void getInputDescription(InputDescription& inputDescription)
 {
  InputLine("pqr", "required", false)
  .addArgument(StringArgument(true, ""))
  .addToInputDescription(inputDescription);

  InputLine("sites", "required", false)
  .addArgument(StringArgument(true, ""))
  .addToInputDescription(inputDescription);

  InputLine("output", "required", false)
  .addArgument(StringArgument(true, ""))
  .addToInputDescription(inputDescription);

  InputLine("dummy", "", false)
  .addToInputDescription(inputDescription);

  InputLine("symmetry_check", "", false)
  .addToInputDescription(inputDescription);

  InputLine("writemap", "", false)
  .addToInputDescription(inputDescription);

  InputLine("diel_atomno", "", true)
  .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 0))
  .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
  .addToInputDescription(inputDescription);

  InputLine("diel_resid", "", true)
  .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 0))
  .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
  .addToInputDescription(inputDescription);

  InputLine("diel_segres", "", true)
  .addArgument(StringArgument(true, ""))
  .addArgument(IntegerArgument(true, numeric_limits<int>::min(), 0))
  .addArgument(FloatArgument(true, numeric_limits<int>::min(), 0.0))
  .addToInputDescription(inputDescription);

  /* membrane <vector-spec> <vector-spec> <dielectric constant>
              ^^^^^^^^^^^^^ |||||||||||||
             support vector |||||||||||||
                           membrane vector
  */
  InputLine("membrane", "", false)
  .addArgument(FloatArgument(true, numeric_limits<int>::min()))
  .addArgument(FloatArgument(true, numeric_limits<int>::min()))
  .addArgument(FloatArgument(true, numeric_limits<int>::min()))
  .addArgument(FloatArgument(true, numeric_limits<int>::min()))
  .addArgument(FloatArgument(true, numeric_limits<int>::min()))
  .addArgument(FloatArgument(true, numeric_limits<int>::min()))
  .addArgument(FloatArgument(true, numeric_limits<int>::min()))
  .addToInputDescription(inputDescription);
 }

/** @brief Run complete titration calculation (potentials, pkints, W-matrix)
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 virtual void run() = 0;

/** @brief Destructor
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 virtual ~Titration() 
 {
  deleteTitratableResidues();
  delete apbsparms;
  delete pqr;
  cout << stopTimer() << endl;
 }

};

/**
 * @brief Legacy titration that performs a MEAD-like computation.
 * @author Gernot Kieseritzky
 */

class Legacy : public Titration
{
protected:

 vector< vector<double> > wmatrix; /** W interaction matrix with MEAD-like structure */

 /* methods */

/** @brief Simply removes partial charges from PQR structure.
 *  @author Gernot Kieseritzky
 */

 virtual Structure* createProteinDielectricVolume();

 virtual void calc_protein_potentials(Structure *proteinVolume, const string& filestem="");
 virtual void calc_model_potentials(Structure *proteinVolume, const string& filestem="");

 virtual void calc_potentials()
 {
  Structure *volume = createProteinDielectricVolume();
  calc_model_potentials(volume);
  calc_protein_potentials(volume);
  delete volume;
 }

 virtual void calc_pka(const string& filestem="");
 virtual void calc_wmatrix(const string& filestem="");
 virtual void symmetrize_wmatrix();
 virtual void write_wmatrix();

public:

 Legacy(InputParser& parser) : Titration(parser) {}

 virtual void run()
 {
  calc_potentials();
  calc_pka();
  calc_wmatrix();
  symmetrize_wmatrix();
  write_wmatrix();
 }

/** @brief Destructor.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 virtual ~Legacy() {}
};

/**
 * @brief Titration supporting multiple titratable states and conformers.
 * Intrinsic pKas and W-Matrix potentials are calculated using the same dielectric volumes
 * that combines all possible conformeric states.
 * @author Gernot Kieseritzky
 */

class MultiConformer : public Legacy
{
protected:

 vector< vector< vector< vector<double> > > > wmatrix;  /** W interaction matrix with Karlsberg 2 compatible structure */

 /* methods overwritten */

/** @brief Combines all conformers to a single dielectric volume.
 *  @author Gernot Kieseritzky
 */

 virtual Structure* createProteinDielectricVolume();

 virtual void calc_pka(const string& filestem="");
 virtual void calc_wmatrix(const string& filestem="");
 virtual void symmetrize_wmatrix();
 virtual void write_wmatrix();

public:

 MultiConformer(InputParser& parser) : Legacy(parser) {}

/** @brief Destructor.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 virtual ~MultiConformer() {}
};

/**
 * @brief Titration supporting multiple titratable states and conformers.
 * Intrinsic pKas potentials are calculated using a dielectric volume with the titratable residue
 * in its current conformeric state and all others in its reference state.
 * W Matrix potentials are calculated using a dielectric volume with the titratable residue
 * in its current conformeric state and all others local conformers combined to single dielectric volume.
 * @author Gernot Kieseritzky
 */

class MultiConformerGunner : public MultiConformer
{
protected:

 /* methods overwritten */

 /** @brief Use dielectric boundary of local conformer.
  *  @author Gernot Kieseritzky
  */

 virtual void calc_intrinsic_protein_potentials(Structure *proteinVolume, const string& filestem="");

 /** @brief Use combined dielectric boundary but take only single conformer for each residue.
  *  @author Gernot Kieseritzky
  */

 virtual void calc_pairwise_protein_potentials(Structure *proteinVolume, const string& filestem="");

 virtual void calc_potentials()
 {
  Structure *volume = createProteinDielectricVolume();
  calc_model_potentials(volume, ".intrinsic");
  calc_intrinsic_protein_potentials(volume, ".intrinsic");
  calc_pairwise_protein_potentials(volume, ".pairwise");
  delete volume;
 }

public:

 MultiConformerGunner(InputParser& parser) : MultiConformer(parser) {}

 virtual void run()
 {
  calc_potentials();
  calc_pka(".intrinsic");
  calc_wmatrix(".pairwise");
  symmetrize_wmatrix();
  write_wmatrix();
 }

/** @brief Destructor.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 virtual ~MultiConformerGunner() {}
};

/**
 * @brief Titration supporting multiple titratable states and conformers.
 * Intrinsic pKas potentials are calculated using a dielectric volume with the titratable residue
 * in its current conformeric state and all others in its reference state.
 * W Matrix potentials are calculated using a dielectric volume with the titratable residue conformers
 * combined to a single dielectric volume (like MultiConformer).
 * @author Gernot Kieseritzky
 */

class MultiConformerTest : public MultiConformerGunner
{
protected:

 virtual void calc_potentials()
 {
  Structure *volume = createProteinDielectricVolume();
  calc_model_potentials(volume, ".intrinsic");
  calc_intrinsic_protein_potentials(volume, ".intrinsic");
  calc_protein_potentials(volume, ".pairwise");
  delete volume;
 }

public:

 MultiConformerTest(InputParser& parser) : MultiConformerGunner(parser) {}

/** @brief Destructor.
 *  @author Gernot Kieseritzky
 *  @ingroup Titration
 */

 virtual ~MultiConformerTest() {}
};
