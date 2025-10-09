#include "titratable.h"
#include "exceptions.h"
#include "utils.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

/* State */

void State::copyMembers(const State& that)
{
 this->id = that.id;
 this->energy = that.energy;
 
 for(Container::const_iterator iter = that.atoms.begin(); iter != that.atoms.end(); iter++)
  this->atoms.insert( Container::value_type(iter->first, new PQR_atom(*(iter->second))) );
}

int State::add_atom(const PQR_atom& atom)
{
 iter_atom it = atoms.find(atom.get_atomid());
 
 if ( it != atoms.end() )
  return 0;
 else
  atoms[atom.get_atomid()] = new PQR_atom(atom);

 return 1;
}

void State::fillMissing(Structure *pqr)
{
 double   undefined_radius = 99.999;
 double   undefined_charge = 99.999;
 Vector3D undefined_coordinate(9999.999, 9999.999, 9999.999);
 
 while( pqr->hasNext_in_selection() )
 {
  PQR_atom& atom = pqr->getNext_in_selection();
  iter_atom it  = atoms.find( atom.get_atomid() );
  if ( it != atoms.end() )
  {
   it->second->set_atomno( atom.get_atomno() );
   it->second->set_resid( atom.get_resid() );
   it->second->set_chain( atom.get_chain() );
   it->second->set_segment( atom.get_segment() );
	if(it->second->get_radius() >= undefined_radius)
    it->second->set_radius( atom.get_radius() );
	if(it->second->get_charge() >= undefined_charge)
	 it->second->set_charge( atom.get_charge() );
	if(it->second->coordinates() == undefined_coordinate)
	 it->second->coordinates() = atom.coordinates();
  }
 }
}

void State::check(Structure *pqr)
{
 iter_atom it;
 
 try
 {
 /* provoke exception */
 for(it=atoms.begin(); it!=atoms.end(); it++)
  pqr->get_atom( it->second->get_atomno() );
 }
 catch(AtomNotFoundException e)
 {
  throw AtomNotFoundException("atom", it->second->get_atomid(), "State::check");
 }
}

void State::delete_atoms()
{
 while( !atoms.empty() )
 {
  iter_atom it = atoms.begin();
  delete it->second;
  atoms.erase(it->first);
 }
}

bool State::isAnionic()
{
 double charge = 0.0;

 for(iter_atom it=atoms.begin(); it!=atoms.end(); it++)
 {
  charge += it->second->get_charge();
 }
 
 return (charge<0);
}

void State::change(Structure *pqr)
{
 while( pqr->hasNext_in_selection() )
 {
  PQR_atom& atom = pqr->getNext_in_selection();
  iter_atom it  = atoms.find( atom.get_atomid() );
  if ( it != atoms.end() )
  {
   double new_charge   = it->second->get_charge();
	double new_radius   = it->second->get_radius();
   Vector3D new_coords = it->second->coordinates();
	/* change charge */
	pqr->get_atom( atom.get_atomno() ).set_charge(new_charge);
	/* change radius */
	pqr->get_atom( atom.get_atomno() ).set_radius(new_radius);
	/* change coordinates */
	pqr->get_atom( atom.get_atomno() ).coordinates() = new_coords;
  }
 }
}

void State::charge(Structure *pqr)
{
 PQR_atom *atom;

 for(iter_atom itat=atoms.begin(); itat!=atoms.end(); itat++)
 {
  atom = itat->second;
  pqr->get_atom( atom->coordinates() ).set_charge( atom->get_charge() );
 }
}

void State::uncharge(Structure *pqr)
{
 PQR_atom *atom;

 for(iter_atom itat=atoms.begin(); itat!=atoms.end(); itat++)
 {
  atom = itat->second;
  pqr->get_atom( atom->coordinates() ).set_charge(0.0);
 }
}

void State::selectCoordinates(Structure *pqr)
{
 // loop over all atoms of this state
 for(iter_atom itat=atoms.begin(); itat!=atoms.end(); itat++)
 {
  // if coordinates still exist select them
  if( pqr->hasPoint( itat->second->coordinates() ) )
   pqr->select_coordinate( itat->second->coordinates() );
 }
}

void State::print(ostream& os)
{
 os.setf(ios_base::fixed);
 os.precision(2);
 os << energy << " " << id << endl;
 for(iter_atom itat=atoms.begin(); itat!=atoms.end(); itat++)
 {
  os << *(itat->second);
 }
}

/* Titratable */

void Titratable::copyMembers(const Titratable& that)
{
 this->name   = that.name;
 this->resid  = that.resid;
 
 for(vector<State*>::const_iterator iter = that.states.begin(); iter != that.states.end(); iter++)
  this->states.push_back( new State(**iter) );
}

Titratable::Titratable(int resid, const string& residue)
{
 ostringstream out;
 out << residue << "-" << resid;
 this->resid = resid;
 this->name  = out.str();
}

void Titratable::parseST(const string& filename)
{
 int lineNr = 0;
 string inputline;
 vector<string> col;
 double energy;
 State *state;
 PQR_atom atom;
 ostringstream log;

 cout << this->name << ": parsing " << filename << " ... ";

 // open file
 ifstream in( filename.c_str() );
 if(!in)
 {
  in.close();
  throw FileIOException(filename, "Titratable::parseST");
 }

 // parse
 while( !in.eof() )
 {
  lineNr++;
  getline(in, inputline);
  col = split(inputline);
  if(col.empty()) continue; /* ignore empty lines */
  /* parse state identifier */
  if ( sscanf(col[0].c_str(), "%lg", &energy) ) /* energy header */
  {
	/* convert model energy */
	if(col[1]=="pK")
	 energy *= PK_TO_KJMOL;
	else
	 throw ParseException(lineNr, filename, "Unknown energy unit.");
	/* add state */
	state = new State(col[2], energy);
	states.push_back(state);
  }
  else if(col[0]=="ATOM") /* ATOM line */
  {
   if( atom.fromPQRLine(inputline) )
	{
	 atom.set_atomno(-1); /* mark atom not be associated with PQR structure yet */
	 if( state->add_atom(atom)!=1 ) 
	  throw ParseException(lineNr, filename, "Duplicate atom in st-File.");
	}
	else
	{
	 in.close();
	 throw ParseException(lineNr, filename, "Charge and coordinate data must be in PQR format. Also check for and remove tab characters.");
	}
  }
  else /* otherwise */
  {
	in.close();
	throw ParseException(lineNr, filename, "Parse error.");
  }
 }
 /* states.push_back(state); */

 defineReferenceState();
 
 cout << "finished." << endl;
}

void Titratable::defineReferenceState()
{
 State *temp;

 for(unsigned int i=0; i<states.size(); i++)
 {
  if ( states[i]->getID()=="R" )       /* state marked as reference */
  {
   temp = states[0];                   /* rescue first state */
	states[0] = states[i];
	states[i] = temp;
	break;
  }
 }
}

void Titratable::switchToState(int state, Structure *pqr)
{
 /* make sure no atoms selected */
 pqr->select_clear();

 /* select residue */
 select(pqr);

 /* switch residue */
 try
 {
  states.at(state)->change(pqr);
 }
 catch (out_of_range& x)
 {
  throw UnknownStateException(state, name, "Titratable::switchToState");
 }

 /* clear selection */
 pqr->select_clear();
 pqr->update();
}

void Titratable::setChargesOfState(int state, Structure *pqr)
{
 /* make sure no atoms selected */
 pqr->select_clear();

 /* select residue */
 select(pqr);

 /* switch residue */
 try
 {
  states.at(state)->charge(pqr);
 }
 catch (out_of_range& x)
 {
  throw UnknownStateException(state, name, "Titratable::setChargesOfState");
 }

 /* clear selection */
 pqr->select_clear();
}

void Titratable::clearChargesOfState(int state, Structure *pqr)
{
 /* make sure no atoms selected */
 pqr->select_clear();

 /* select residue */
 select(pqr);

 /* switch residue */
 try
 {
  states.at(state)->uncharge(pqr);
 }
 catch (out_of_range& x)
 {
  throw UnknownStateException(state, name, "Titratable::clearChargesOfState");
 }

 /* clear selection */
 pqr->select_clear();
}

void Titratable::fillMissing(Structure *pqr)
{
 unsigned int state = 0;

 /* make sure no atom is selected */
 pqr->select_clear();

 try
 {

 select(pqr); 
 for(state=0; state<states.size(); state++)
 {
  states[state]->fillMissing(pqr);
  states[state]->check(pqr);
 }
 pqr->select_clear();

 } 
 catch(AtomNotFoundException& e)
 {
  pqr->select_clear();
  string stateDecription(getName()+": "+states[state]->getID());
  throw GeneralException(e.getErrorCode(), "Unknown atoms in "+stateDecription+"\n"+e.getMessage());
 }
}

void Titratable::selectCoordinatesOfState(int state, Structure *pqr)
{
 try
 {
  states.at(state)->selectCoordinates(pqr);
 }
 catch (out_of_range& x)
 {
  throw UnknownStateException(state, name, "Titratable::clearChargesOfState");
 }
}

void Titratable::delete_states()
{
 while( !states.empty() )
 {
  iter_states itstat = states.begin();
  delete *itstat;
  states.erase(itstat);
 }
}

void Titratable::printStates(ostream& os)
{
 for(iter_states itst=states.begin(); itst!=states.end(); itst++)
  (*itst)->print(os);
}

double Titratable::getEnergy(int state)
{
 double energy = 0.0;

 try
 {
  energy = states.at(state)->getEnergy();
 }
 catch (out_of_range& x)
 {
  throw UnknownStateException(state, name, "Titratable::getEnergy");
 }
 
 return energy;
}

const string& Titratable::getID(int state)
{
 try
 {
  return states.at(state)->getID();
 }
 catch (out_of_range& x)
 {
  throw UnknownStateException(state, name, "Titratable::getEnergy");
 }
}

bool Titratable::isAnionic(int state)
{
 bool answer = false;

 try
 {
  answer = states.at(state)->isAnionic();
 }
 catch (out_of_range& x)
 {
  throw UnknownStateException(state, name, "Titratable::switchToReferenceAndClearCharges");
 }
 
 return answer;
}

McTitratable::McTitratable(int resid, const string& residue, const string& segname) : Titratable(resid, residue)
{
 ostringstream out;
 out << name << "_" << segname;
 this->name  = out.str();
 this->segname = segname;
}
