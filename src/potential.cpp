#include "potential.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "exceptions.h"

Potential::Potential(const Potential& pot)
{
 values = pot.values;
 potIterator = values.begin();
}

double Potential::get_value_at(const Vector3D& v) const
{
 const_iter_pot iter = values.find(v);
 
 if ( iter==values.end() )
  throw PointNotFoundException(v.get_x_value(), v.get_y_value(), v.get_z_value(), "Potential::get_value_at");

 return iter->second;
}

double Potential::get_value_at(double x, double y, double z) const
{
 Vector3D v(x, y, z);
  
 return get_value_at(v);
}

void Potential::set_value_at(double x, double y, double z, double pot)
{
 Vector3D v(x, y, z);
 values[v] = pot;
 potIterator = values.begin();
}

void Potential::copy_to_arrays(double x[], double y[], double z[], double pot[]) const
{
 int i = 0;
 
 for(const_iter_pot iter=values.begin(); iter!=values.end(); iter++)
 {
  x[i]  = iter->first.get_x_value();
  y[i]  = iter->first.get_y_value();
  z[i]  = iter->first.get_z_value();
  pot[i]= iter->second;
  i++;
 }
}

void Potential::copy_from_arrays(int size, double x[], double y[], double z[], double pot[])
{
 for(int i=0; i<size; i++)
  set_value_at(x[i], y[i], z[i], pot[i]);
}

bool Potential::hasNext()
{
 bool answer = potIterator!=values.end();
 
 // rewind iterator
 if (answer == false)
  potIterator = values.begin();

 return answer;
}

struct Potential::extAccess Potential::getNext()
{
 struct extAccess ea;
 
 if( potIterator != values.end() )
 {
  ea.x   = potIterator->first.get_x_value();
  ea.y   = potIterator->first.get_y_value();
  ea.z   = potIterator->first.get_z_value();
  ea.pot = potIterator->second;
  potIterator++;
 }
 else
  throw PointNotFoundException(99., 99., 99., "Potential::getNext");

 return ea;
}

/* PotentialsOfResdidue */

double PotentialsOfResidue::calcBorn(int state, Structure *mstat, Structure *mref)
{
 double born0 = 0, born1 = 0;
 Potential *pp_s0, *pp_s1, *mp_s0, *mp_s1;
 Satom charge_s0, charge_s1;

 // get potentials

 try
 {

 pp_s0 = protein.at(0);
 mp_s0 = model.at(0);
 pp_s1 = protein.at(state);
 mp_s1 = model.at(state);
 
 while( mref->hasNext() & mstat->hasNext() ) // & gives same result as && in this case and executes both operands
 {
  // get charges
  charge_s0 = mref->getNext();
  charge_s1 = mstat->getNext();

  // ignore charges that are zero
  if( charge_s1.get_charge()!=0.0 )
  {
   born1 += charge_s1.get_charge() * 
             ( 
                pp_s1->get_value_at( charge_s1.coordinates() ) 
              - mp_s1->get_value_at( charge_s1.coordinates() )
             );
  }

  // ignore charges that are zero
  if( charge_s0.get_charge()!=0.0 )
  {
   born0 += charge_s0.get_charge() * 
             ( 
                pp_s0->get_value_at( charge_s0.coordinates() ) 
              - mp_s0->get_value_at( charge_s0.coordinates() )
             );
  }
 }

 } 
 catch (out_of_range& e)
 {
  throw UnknownStateException(state, "", "PotentialsOfResidue::calcBorn");
 }
 catch (PointNotFoundException& e)
 {
  ostringstream ostr;
  ostr << "PotentialsOfResidue::calcBorn caught a PointNotFoundException:" << endl;
  ostr << "State: " << state << endl;
  ostr << "mstat: " << endl;
  mstat->print_pqr(ostr);
  ostr << "mref: " << endl;
  mref->print_pqr(ostr);
  ostr << e.getMessage();
  throw GeneralException(e.getErrorCode(), ostr.str());
 }

 return 0.5*(born1-born0);
}

double PotentialsOfResidue::calcBack(int state, Structure *pback, Structure *mback)
{
 double back1 = 0, back2 = 0;
 Potential *pp_s0, *pp_s1, *mp_s0, *mp_s1;
 Satom pcharge, mcharge;

 try
 {

 // get potentials
 pp_s0 = protein.at(0);
 mp_s0 = model.at(0);
 pp_s1 = protein.at(state);
 mp_s1 = model.at(state);

 while( pback->hasNext() )
 {
  // get charges
  pcharge = pback->getNext();

  // ignore charges that are zero
  if(pcharge.get_charge()==0.0) continue;
  
  back1 += pcharge.get_charge() * 
            ( 
               pp_s1->get_value_at( pcharge.coordinates() )
             - pp_s0->get_value_at( pcharge.coordinates() )
            );
 }

 while( mback->hasNext() )
 {
  // get charges
  mcharge = mback->getNext();
  
  // ignore charges that are zero
  if(mcharge.get_charge()==0.0) continue;
  
  // back term
  back2 += mcharge.get_charge() * 
            ( 
               mp_s1->get_value_at( mcharge.coordinates() )
             - mp_s0->get_value_at( mcharge.coordinates() )
            );
 }

 } 
 catch (out_of_range& e)
 {
  throw UnknownStateException(state, "", "PotentialsOfResidue::calcBack");
 }
 catch (PointNotFoundException& e)
 {
  ostringstream ostr;
  ostr << "PotentialsOfResidue::calcBack caught a PointNotFoundException:" << endl;
  ostr << "State: " << state << endl;
  ostr << "pback: " << endl;
  pback->print_pqr(ostr);
  ostr << "mback: " << endl;
  mback->print_pqr(ostr);
  ostr << e.getMessage();
  throw GeneralException(e.getErrorCode(), ostr.str());
 }

 return (back1-back2);
}

double PotentialsOfResidue::calcWmv(int state, int scale, Structure *mstat, Structure *mref)
{
 double wmv = 0;
 Potential *pp_s0, *pp_s1;
 Satom charge_s0, charge_s1;

 try
 {

 // get potentials
 pp_s0 = protein.at(0);
 pp_s1 = protein.at(state);

 while( mref->hasNext() & mstat->hasNext() ) // & gives same result as && in this case and executes both operands
 { 
  // get charges
  charge_s0 = mref->getNext();
  charge_s1 = mstat->getNext();
  
  // ignore charges that are zero
  if ( charge_s0.get_charge()==0.0 && charge_s1.get_charge()==0.0 )
   continue;
  
  /*
  
  wmv += scale * ( charge_s1.get_charge() - charge_s0.get_charge() ) *
            (
               pp_s1->get_value_at( charge_s1.coordinates() )
             - pp_s0->get_value_at( charge_s0.coordinates() )
            );
 */

 wmv += scale * ( charge_s1.get_charge() *
                 (  pp_s1->get_value_at(charge_s1.coordinates()) 
                  - pp_s0->get_value_at(charge_s1.coordinates())
                 )
					  - charge_s0.get_charge() * 
                 (  pp_s1->get_value_at(charge_s0.coordinates()) 
                  - pp_s0->get_value_at(charge_s0.coordinates()) 
                 )
					 );

 }

 } 
 catch (out_of_range& e)
 {
  throw UnknownStateException(state, "", "PotentialsOfResidue::calcWmv");
 }
 catch (PointNotFoundException& e)
 {
  ostringstream ostr;
  ostr << "PotentialsOfResidue::calcWmv caught a PointNotFoundException:" << endl;
  ostr << "State: " << state << endl;
  ostr << "mstat: " << endl;
  mstat->print_pqr(ostr);
  ostr << "mref: " << endl;
  mref->print_pqr(ostr);
  ostr << e.getMessage();
  throw GeneralException(e.getErrorCode(), ostr.str());
 }

 return wmv;
}

int PotentialsOfResidue::readFromFile(const string& filename)
{
 iter_pot iter;
 int states; 
 int nrOfValues;
 struct Potential::extAccess ea;
 
 // if there is no file return with 0
 ifstream in( filename.c_str() );
 if (!in)
 {
  in.close();
  return 0;
 }
 
 // delete old potentials
 deletePotentials();
 
 // protein potentials
 in.read( reinterpret_cast<char *>(&states), sizeof(int) );
 for( int state=0; state<states; state++ )
 {
  in.read( reinterpret_cast<char *>(&nrOfValues), sizeof(int) );
  if(nrOfValues>0) protein.push_back( new Potential() );
  for( int i=0; i<nrOfValues; i++ )
  {
   in.read( reinterpret_cast<char *>(&(ea.x)), sizeof(double) );
   in.read( reinterpret_cast<char *>(&(ea.y)), sizeof(double) );
   in.read( reinterpret_cast<char *>(&(ea.z)), sizeof(double) );
   in.read( reinterpret_cast<char *>(&(ea.pot)), sizeof(double) );
	protein[state]->set_value_at( ea.x, ea.y, ea.z, ea.pot );
  }
 }

 // model potentials
 in.read( reinterpret_cast<char *>(&states), sizeof(int) );
 for( int state=0; state<states; state++ )
 {
  in.read( reinterpret_cast<char *>(&nrOfValues), sizeof(int) );
  if(nrOfValues>0) model.push_back( new Potential() );
  for( int i=0; i<nrOfValues; i++ )
  {
   in.read( reinterpret_cast<char *>(&(ea.x)), sizeof(double) );
   in.read( reinterpret_cast<char *>(&(ea.y)), sizeof(double) );
   in.read( reinterpret_cast<char *>(&(ea.z)), sizeof(double) );
   in.read( reinterpret_cast<char *>(&(ea.pot)), sizeof(double) );
	model[state]->set_value_at( ea.x, ea.y, ea.z, ea.pot );
  }
 }

 in.close();
 
 return 1;
}

void PotentialsOfResidue::writeToFile(const string& filename)
{
 iter_pot iter;
 int states; 
 int nrOfValues;
 struct Potential::extAccess ea;

 ofstream out( filename.c_str() );
 if (!out)
 {
  out.close();
  throw FileIOException(filename, "PotentialsOfResidue::writeToFile");
 }

 // protein potentials
 states = protein.size();
 out.write( reinterpret_cast<char *>(&states), sizeof(int) );
 for( iter=protein.begin(); iter!=protein.end(); iter++ )
 {
  nrOfValues = (**iter).nrOfPoints();
  out.write( reinterpret_cast<char *>(&nrOfValues), sizeof(int) );
  while( (**iter).hasNext() )
  {
   ea = (**iter).getNext();
   out.write( reinterpret_cast<char *>(&(ea.x)), sizeof(double) );
	out.write( reinterpret_cast<char *>(&(ea.y)), sizeof(double) );
	out.write( reinterpret_cast<char *>(&(ea.z)), sizeof(double) );
	out.write( reinterpret_cast<char *>(&(ea.pot)), sizeof(double) );
  }
 }

 // model potential 
 states = model.size();
 out.write( reinterpret_cast<char *>(&states), sizeof(int) );
 for( iter=model.begin(); iter!=model.end(); iter++ )
 {
  nrOfValues = (**iter).nrOfPoints();
  out.write( reinterpret_cast<char *>(&nrOfValues), sizeof(int) );
  while( (**iter).hasNext() )
  {
   ea = (**iter).getNext();
   out.write( reinterpret_cast<char *>(&(ea.x)), sizeof(double) );
	out.write( reinterpret_cast<char *>(&(ea.y)), sizeof(double) );
	out.write( reinterpret_cast<char *>(&(ea.z)), sizeof(double) );
	out.write( reinterpret_cast<char *>(&(ea.pot)), sizeof(double) );
  }
 }
 
 out.close();
}

void PotentialsOfResidue::deletePotentials()
{
 iter_pot iter;
 
 // protein
 for( iter=protein.begin(); iter!=protein.end(); iter++ )
  delete *iter;

 // model
 for( iter=model.begin(); iter!=model.end(); iter++ )
  delete *iter;
  
 protein.clear();
 model.clear();
}

const Potential& PotentialsOfResidue::getProteinStatePotential(int state)
{
 try
 {
  return *(protein.at(state));
 }
 catch (out_of_range& e)
 {
  throw UnknownStateException(state, "", "PotentialsOfResidue::getProteinStatePotential");
 }
}

const Potential& PotentialsOfResidue::getModelStatePotential(int state)
{
 try
 {
  return *(model.at(state));
 }
 catch (out_of_range& e)
 {
  throw UnknownStateException(state, "", "PotentialsOfResidue::getModelStatePotential");
 }
}
