/* PQR.cpp */
#include <stdexcept>
#include "PQR.h"
#include "utils.h"
#include "exceptions.h"

/* PQR_atom */

int PQR_atom::instanceCounter = 0;

void PQR_atom::PQR(ostream& os, PQR_atom& atom)
{
 os.setf(ios_base::fixed);
 os.width(6);
 os << left << "ATOM";
 os.width(5);
 os << right << atom.get_atomno();

 os << " "; /* spacer */

 os.width(4);
 os << left << atom.get_atomid();

 os << atom.get_conf(); /* alt indicator */

 os.width(3);
 os << atom.get_residue();

 os << " "; /* spacer */

 os << atom.get_chain(); /* chain id */
 os.width(4);
 os << right << atom.get_resid();
 os.width(1);
 os << " "; /* icode */

 os << "   "; /* spacer */

 os.width(8);
 os.precision(3);
 os << right << atom.coordinates().get_x_value();
 os.width(8);
 os.precision(3);
 os << right << atom.coordinates().get_y_value();
 os.width(8);
 os.precision(3);
 os << right << atom.coordinates().get_z_value();

 os.width(6);
 os.precision(3); // os.precision(2);
 os << right << atom.get_charge();

 os.width(6);
 os.precision(3); //  os.precision(2);
 os << atom.get_radius();

 os << "      "; /* spacer */

 os.width(4);
 os << left << atom.get_segment();

 os << endl;
}

PQR_atom::PQR_atom(const string& pqrline)
{
 fromPQRLine(pqrline);
 instanceCounter++;
}

int PQR_atom::fromPQRLine(const string& pqrline)
{
 string temp;
 double x, y, z;
 
 // check string length
 // should be at least 75 characters long (we're ignoring the element name)
 if ( pqrline.size() < 75 )
  return 0;

 // atom number
 temp = pqrline.substr(6,5);
 if( sscanf(temp.c_str(), "%d", &(this->atomno))!=1 )
  return 0;

 // atom name
 temp = pqrline.substr(12,4);
 this->atomid = split(temp).at(0);

 // conformer id
 temp = pqrline.substr(16,1);
 if( sscanf(temp.c_str(), "%c", &(this->altloc))!=1 )
  this->altloc = ' ';

 // residue name
 this->residue = pqrline.substr(17,3);

 // chain id
 temp = pqrline.substr(21,1);
 if( sscanf(temp.c_str(), "%c", &(this->chain))!=1 )
  this->chain = ' ';

 // residue number
 temp = pqrline.substr(22,4);
 if( sscanf(temp.c_str(), "%i", &(this->resid))!=1 )
  return 0;

 // coordinates
 temp = pqrline.substr(30,8);
 if( sscanf(temp.c_str(), "%lg", &x)!=1 )
  return 0;
 temp = pqrline.substr(38,8);
 if( sscanf(temp.c_str(), "%lg", &y)!=1 )
  return 0;
 temp = pqrline.substr(46,8);
 if( sscanf(temp.c_str(), "%lg", &z)!=1 )
  return 0;
 this->coordinates().set(x, y, z);

 // charge
 temp = pqrline.substr(54,6);
 if( sscanf(temp.c_str(), "%lg", &(this->charge))!=1 )
  return 0;
 
 // radius
 temp = pqrline.substr(60,6);
 if( sscanf(temp.c_str(), "%lg", &(this->radius))!=1 )
  return 0;

 // segment name
 temp = pqrline.substr(72,4);
 try
 {
  this->segname = split(temp).at(0);
 }
 catch(out_of_range& x)
 {
  this->segname = "    ";
 }

 return 1;
}

PQR_atom::PQR_atom(int atomno, int resid, double charge, double radius, const string& residue, const string& atomid)
{
 this->atomno = atomno;
 this->atomid = atomid;
 this->altloc = ' ';
 this->residue= residue;
 this->chain  = ' ';
 this->resid  = resid;
 this->charge = charge;
 this->radius = radius;
 this->segname= " ";
 instanceCounter++;
}

ostream& operator<<(ostream& os, PQR_atom& atom)
{
 PQR_atom::PQR(os, atom);
 return os; 
}

/* Satom */

Satom::Satom(int atomno, int resid, double charge, double radius, double soluteDiel, const string& residue, const string& atomid) : PQR_atom(atomno, resid, charge, radius, residue, atomid)
{
 this->soluteDiel = soluteDiel;
}

/* Structure */

void Structure::copyMembers(const Structure& that)
{
 Satom *satom;
 for(AtomContainer::const_iterator iter=that.atoms.begin(); iter!=that.atoms.end(); iter++)
 {
  satom = *iter;
  add_atom( new Satom(*satom) );
 }
}

void Structure::add_atom(Satom* atom)
{
 atoms.push_back(atom);
 atomIterator=atoms.begin();
}

Structure::Structure(const Structure& that)
{
 copyMembers(that);
 update();
}

Structure* Structure::fromSelection()
{
 Structure* clone = new Structure();
 
 for(iter_atom i=atomSelection.begin(); i!=atomSelection.end(); i++)
  clone->add_atom( new Satom(**i) );

 clone->update();
 
 return clone;
}

void Structure::fromPQR(const string& pqrfile)
{
 string line;
 Satom *satom;
 int lineNr=0;

 cout << "Parsing PQR file " << pqrfile << " ... ";

 ifstream in( pqrfile.c_str() );

 if (!in)
 {
  in.close();
  throw FileIOException(pqrfile, "Structure::fromPQR");
 }
 
 while( !in.eof() )
 {
  getline(in, line);
  lineNr++;
  // ignore lines that do not begin with ATOM
  if ( line.find("ATOM", 0)==string::npos ) continue;
  satom = new Satom();
  if ( satom->fromPQRLine(line)!=1 )
  {
   in.close();
	throw ParseException(lineNr, pqrfile, "This doesn't seem to be a valid PQR file!");
  }
  add_atom(satom);
 }
 
 in.close();

 update();

 cout << "finished." << endl;
}

Structure* Structure::createModelFromSelection()
{
 Structure *model;
 int resid;
 string segment;
 iter_resid itres;
 
 /* copy current selection */
 model = fromSelection();
 
 /* get selected resid and segment name */
 if( hasNext_in_selection() )
 {
  Satom atom = getNext_in_selection();
  resid = atom.get_resid();
  segment = atom.get_segment();
 }
 select_clear();
 
 /* look for previous residue in chain */
 itres = atomsBySegname[segment].find(resid-1);
 if(itres!=atomsBySegname[segment].end())
 {
  select_resid(segment, resid-1);
  while( hasNext_in_selection() )
  {
   Satom atom = getNext_in_selection();
   if ( atom.get_atomid()=="C" || atom.get_atomid()=="O" )
    model->add_atom( new Satom(atom) );
  }
  select_clear();
 }

 /* look for next residue in chain */
 itres = atomsBySegname[segment].find(resid+1);
 if(itres!=atomsBySegname[segment].end())
 {
  select_resid(segment, resid+1);
  while( hasNext_in_selection() )
  {
   Satom atom = getNext_in_selection();
   if ( atom.get_atomid()=="N" || atom.get_atomid()=="HN" || atom.get_atomid()=="CA" )
    model->add_atom( new Satom(atom) );
  }
  select_clear();
 }
 
 model->update();
 
 return model;
}

void Structure::refill_containers()
{
 iter_coord itvec;
 iter_segname itseg;
 string segname;

 atomsByCoord.clear();
 atomsBySegname.clear();

 for(iter_atom i=atoms.begin(); i!=atoms.end(); i++)
 {
  /* We also check whether two atoms have the same coordinates. */
  itvec = atomsByCoord.find( (**i).coordinates() );
  if( itvec!=atomsByCoord.end() )
  {
   this->print_pqr(cout);
   throw DuplicateCoordinatesException((**i).get_atomno(), itvec->second->get_atomno(), "Structure::refill_containers");
  }
  else
  {
   atomsByCoord.insert( pvec((**i).coordinates(), *i) );
  }
  
  /* store segnames */
  segname = (**i).get_segment();
  itseg = atomsBySegname.find(segname);
  pair<iter_segname, bool> p = atomsBySegname.insert( psegname(segname, AtomsByResid()) );
  p.first->second.insert( presid((**i).get_resid(), *i) );
 }
}

void Structure::renumber()
{
 int atomno = 0;
 int resid = -1;
 int oldResid = -1;
 residues = 0;

 atoms.clear();
 atomSelection.clear();

 for(iter_segname seg=atomsBySegname.begin(); seg!=atomsBySegname.end(); seg++) /* for all segments */
 {
  for(iter_resid res=seg->second.begin(); res!=seg->second.end(); res++) /* for all residues in segment */
  {
   resid = res->second->get_resid(); /* res->first */
	if(resid!=oldResid)
	{
	 residues++;
	 oldResid = resid;
	}
   res->second->set_atomno(++atomno);
   atoms.push_back(res->second);
  }
 }

 atomIterator = atoms.begin();
 seleIterator = atomSelection.begin();
}

void Structure::append(const Structure& other)
{
 copyMembers(other);
}

void Structure::appendAtUnknownCoordinates(const Structure& other)
{
 Satom *other_atom;
 iter_coord it_my_atom;
 
 for(AtomContainer::const_iterator iter=other.atoms.begin(); iter!=other.atoms.end(); iter++)
 {
  other_atom = *iter;
  it_my_atom = atomsByCoord.find( other_atom->coordinates() );
  
  if( it_my_atom==atomsByCoord.end() )  /* coordinates unknown */
  {
   other_atom = new Satom(*other_atom); /* make a copy of other atom */
   add_atom( other_atom ); /* add atom to atom container */
	atomsByCoord.insert( pvec(other_atom->coordinates(), other_atom) ); /* add to coordinates map */
  }
  else /* coordinates known */
  {
   double my_radius = it_my_atom->second->get_radius();
	double other_radius = other_atom->get_radius();
	if(other_radius > my_radius) /* if other radius larger, inflate my radius */
	 it_my_atom->second->set_radius(other_radius);
  }
 }
}

/** low level manipulations */

void Structure::delete_atoms_in_selection()
{
 AtomContainer atomsNew;
 set<int> toErase;
 
 // first collect atom numbers to be deleted
 for(iter_atom iter=atomSelection.begin(); iter!=atomSelection.end(); iter++)
 {
  toErase.insert( (**iter).get_atomno() );
 }
 
 // create new atom container
 for(iter_atom iter=atoms.begin(); iter!=atoms.end(); iter++)
 {
  if( toErase.find( (**iter).get_atomno() ) != toErase.end() ) // atomno marked for deletion
   delete *iter;
  else
   atomsNew.push_back(*iter);
 }
 atoms.clear();
 atoms = atomsNew;
 atomSelection.clear();
 update();
}

void Structure::delete_all_atoms()
{
 for(iter_atom iter=atoms.begin(); iter!=atoms.end(); iter++)
  delete *iter;
 atoms.clear();
 update();
}

Satom& Structure::get_atom(int atomno) const
{
 Satom *atom;
 
 try
 {
  atom = atoms.at(atomno-1);
 }
 catch(out_of_range& x)
 {
  throw AtomNotFoundException("atom", atomno, "Structure::get_atom");
 }

 return *atom;
}

Satom& Structure::get_atom(const Vector3D& v) const
{
 const_iter_coord itvec = atomsByCoord.find(v);
 
 if ( itvec==atomsByCoord.end() )
  throw PointNotFoundException(v.get_x_value(), v.get_y_value(), v.get_z_value(), "Structure::get_atom");

 return *( itvec->second );
}

/* selection operators */

void Structure::select_atomno(int atomno)
{ 
 try
 {
  atomSelection.push_back( atoms.at(atomno-1) ); // -1 vector is 0-based
 }
 catch(out_of_range& x)
 {
  throw AtomNotFoundException("atom", atomno, "Structure::select_atomno");
 }
 
 // init atom iterator
 seleIterator = atomSelection.begin();
}

void Structure::select_coordinate(const Vector3D& v)
{
 iter_coord itvec = atomsByCoord.find(v);
 
 if ( itvec==atomsByCoord.end() )
  throw PointNotFoundException(v.get_x_value(), v.get_y_value(), v.get_z_value(), "Structure::select_coordinate");
 else
  atomSelection.push_back(itvec->second);
  //atomSelection[itvec->second->get_atomno()] = itvec->second;

 // init atom iterator
 seleIterator = atomSelection.begin();
}

int Structure::select_resid(int resid)
{
 int count  = 0;
 int atomno = 0;

 iter_resid b = atomsBySegname.begin()->second.lower_bound(resid);
 iter_resid e = atomsBySegname.begin()->second.upper_bound(resid);

 if (b==e)
  throw AtomNotFoundException("residue", resid, "Structure::select_resid");

 for(iter_resid iter=b; iter!=e; iter++)
 {
  atomSelection.push_back(iter->second);
  count++;
 }

 // init atom iterator
 seleIterator = atomSelection.begin();
 
 return count;
}

int Structure::select_resid(const string& segname, int resid)
{
 int count  = 0;
 int atomno = 0;
 
 iter_segname seg = atomsBySegname.find(segname);
 
 if(seg==atomsBySegname.end())
  throw AtomNotFoundException("segment name", segname, "Structure::select_resid");
 
 iter_resid b = seg->second.lower_bound(resid);
 iter_resid e = seg->second.upper_bound(resid);

 if(b==e)
  throw AtomNotFoundException("residue", resid, "Structure::select_resid");

 for(iter_resid iter=b; iter!=e; iter++)
 {
  atomSelection.push_back(iter->second);
  count++;
 }

 // init atom iterator
 seleIterator = atomSelection.begin();
 
 return count;
}

void Structure::select_all()
{
 atomSelection = atoms;

 // init atom iterator
 seleIterator = atomSelection.begin();
}

void Structure::select_clear()
{
 atomSelection.clear();

 // init atom iterator
 seleIterator = atomSelection.begin();
}

bool Structure::hasNext()
{
 bool answer = atomIterator != atoms.end();
 
 // rewind iterator
 if (answer==false)
  atomIterator = atoms.begin();

 return answer;
}

bool Structure::hasNext_in_selection()
{
 bool answer = seleIterator != atomSelection.end();
 
 // rewind iterator
 if (answer==false)
  seleIterator = atomSelection.begin();

 return answer;
}

Satom& Structure::getNext()
{
 Satom *atom;
 
 if( atomIterator != atoms.end() )
 {
  atom = *atomIterator;
  atomIterator++;
 }
 
 return (*atom);
}

Satom& Structure::getNext_in_selection()
{
 Satom *atom;
 
 if( seleIterator != atomSelection.end() )
 {
  atom = *seleIterator;
  seleIterator++;
 }
  
 return (*atom);
}

void Structure::clear_charges_in_selection()
{
 for(iter_atom iter=atomSelection.begin(); iter!=atomSelection.end(); iter++)
  (**iter).set_charge(0.0);
}

void Structure::clear_all_charges()
{
 for(iter_atom iter=atoms.begin(); iter!=atoms.end(); iter++)
  (**iter).set_charge(0.0);
}

/** coordinate manipulations */

void Structure::calc_geometry()
{
 double xmin=0., ymin=0., zmin=0.;
 double xmax=0., ymax=0., zmax=0.;
 double sxmin=0., symin=0., szmin=0.;
 double sxmax=0., symax=0., szmax=0.;
 double cx=0., cy=0., cz=0.;
 double x, y, z, rad;
 
 if ( atoms.begin()!=atoms.end() )
 {
  iter_atom i = atoms.begin();
  xmin = (**i).coordinates().get_x_value(); xmax = xmin;
  ymin = (**i).coordinates().get_y_value(); ymax = ymin;
  zmin = (**i).coordinates().get_z_value(); zmax = zmin;

  sxmin = (**i).coordinates().get_x_value() - (**i).get_radius();
  sxmax = (**i).coordinates().get_x_value() + (**i).get_radius();
  symin = (**i).coordinates().get_y_value() - (**i).get_radius();
  symax = (**i).coordinates().get_y_value() + (**i).get_radius();
  szmin = (**i).coordinates().get_z_value() - (**i).get_radius();
  szmax = (**i).coordinates().get_z_value() + (**i).get_radius();
 
  for(;i!=atoms.end(); i++)
  {
   x   = (**i).coordinates().get_x_value();
   y   = (**i).coordinates().get_y_value();
   z   = (**i).coordinates().get_z_value();
	rad = (**i).get_radius();

   /* look for smallest and biggest value of x, y and z */
   if ( x <= xmin ) xmin = x;
   if ( y <= ymin ) ymin = y;
   if ( z <= zmin ) zmin = z;
   if ( x >= xmax ) xmax = x;
   if ( y >= ymax ) ymax = y;
   if ( z >= zmax ) zmax = z;

   /* look for smallest and biggest value of x+/-rad, y+/-rad and z+/-rad */
   if ( x-rad <= sxmin ) sxmin = x-rad;
   if ( y-rad <= symin ) symin = y-rad;
   if ( z-rad <= szmin ) szmin = z-rad;
   if ( x+rad >= sxmax ) sxmax = x+rad;
   if ( y+rad >= symax ) symax = y+rad;
   if ( z+rad >= szmax ) szmax = z+rad;
  }
 
  /* geometric middle of min/max is center */
  cx = 0.5*(xmax + xmin);
  cy = 0.5*(ymax + ymin);
  cz = 0.5*(zmax + zmin);
 }

 /* update properties */
 min    = Vector3D(xmin, ymin, zmin);
 max    = Vector3D(xmax, ymax, zmax);
 smin   = Vector3D(sxmin, symin, szmin);
 smax   = Vector3D(sxmax, symax, szmax);
 center = Vector3D(cx, cy, cz);
}

void Structure::get_center(double *x, double *y, double *z) const
{
 *x = center.get_x_value();
 *y = center.get_y_value();
 *z = center.get_z_value();
}

void Structure::get_minmax(double min[3], double max[3]) const
{
 min[0] = this->min.get_x_value(); min[1] = this->min.get_y_value(); min[2] = this->min.get_z_value();
 max[0] = this->max.get_x_value(); max[1] = this->max.get_y_value(); max[2] = this->max.get_z_value();
}

void Structure::rotate(Vector3D& axis, double deg)
{
 for(iter_atom i=atoms.begin(); i!=atoms.end(); i++)
 {
  (**i).coordinates().rotate(axis, deg);
 }
 
 calc_geometry();
}

void Structure::translate(Vector3D& axis)
{
 for(iter_atom i=atoms.begin(); i!=atoms.end(); i++)
 {
  (**i).coordinates().translate(axis);
 }

 calc_geometry();
}

void Structure::set_soluteDiel_in_selection(double epsp)
{
 for(iter_atom i=atomSelection.begin(); i!=atomSelection.end(); i++)
  (**i).set_soluteDiel(epsp);
}

/* output */

void Structure::print_pqr(ostream& os)
{
 for(iter_atom i=atoms.begin(); i!=atoms.end(); i++)
 {
  os << **i;
 }
 os << "END" << endl;
}

void Structure::print_pqr_selection(ostream& os)
{
 for(iter_atom i=atomSelection.begin(); i!=atomSelection.end(); i++)
 {
  os << **i;
 }
 os << "END" << endl;
}

void Structure::writeToFile_pqr(const string& filename)
{
 // setup output file
 ofstream out( filename.c_str() );
 if (!out)
 {
  out.close();
  throw FileIOException(filename, "Structure::writeToFile_pqr");
 }

 // print pqr into output stream
 print_pqr(out); 
 out.close();  
}
