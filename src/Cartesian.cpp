#include "Cartesian.h"

/* -----------------------------Vector3D---------------------------- */

Vector3D::Vector3D()
{
 this->x = 0e0;
 this->y = 0e0;
 this->z = 0e0;
}

Vector3D::Vector3D(double x, double y, double z)
{
 this->x = x;
 this->y = y;
 this->z = z;
}

void Vector3D::set(double x, double y, double z)
{
 this->x = x;
 this->y = y;
 this->z = z;
}

double Vector3D::abs()
{
 return ( sqrt(x*x+y*y+z*z) );
}

double Vector3D::angle(Vector3D& v)
{
 double scalprod = (*this)*v;
 double abs_sqr  = (this->abs()) * (v.abs());
 double cos = scalprod/abs_sqr;
 double rad = acos(cos);
 return rad2deg(rad);
}

void Vector3D::rotate(Vector3D& axis, double deg)
{
 Matrix *m = Matrix::rotation_matrix( axis, deg2rad(deg) );
  
 /* calculate m*v */
 *this = (*m)*(*this);
 
 /* clean up */
 delete m;
}

void Vector3D::translate(Vector3D& axis)
{
 *this = (*this)+axis;
}

double Vector3D::deg2rad(double deg)
{
 return(deg/180*PI);
}

double Vector3D::rad2deg(double rad)
{
 return(rad/PI*180);
}

/* operators */

bool Vector3D::operator==(const Vector3D& v)
{
 return (x == v.get_x_value() && y == v.get_y_value() && z == v.get_z_value() );
}

bool Vector3D::operator<(const Vector3D& v) const
{
 /* 
  An arbitrary definition since
  there is no natural definition.
  It's defined such that points in 
  space will be sorted by x, y, z
  in increasing order:
  
  x1, y1, z1
  x2, y1, z1
  ...
  xN, y1, z1
  x1, y2, z1
  ...
  xN, yN, z1
  x1, y1, z2
  ...
  xN, yN, zN
 */
 
 bool answer = false;
 
 if ( x < v.get_x_value() )
  answer = true;
 else if ( x==v.get_x_value() && y < v.get_y_value() )
  answer = true;
 else if ( x==v.get_x_value() && y==v.get_y_value() && z < v.get_z_value() )
  answer = true;
 else
  answer = false;
 
 return answer;
}

Vector3D Vector3D::operator+ (const Vector3D& b) const
{
 return Vector3D( x + b.get_x_value(), y + b.get_y_value(), z + b.get_z_value() );
}

Vector3D Vector3D::operator- (const Vector3D& b) const
{
 return Vector3D( x - b.get_x_value(), y - b.get_y_value(), z - b.get_z_value() );
}

double Vector3D::operator* (const Vector3D& b) const
{
 return ( x*b.get_x_value() + y*b.get_y_value() + z*b.get_z_value() );
}

Vector3D Vector3D::cross(const Vector3D& b) const
{
 double cx = y*b.get_z_value() - z*b.get_y_value();
 double cy = z*b.get_x_value() - x*b.get_z_value();
 double cz = x*b.get_y_value() - y*b.get_x_value();
 return Vector3D(cx, cy, cz);
}

Vector3D operator* (double val, const Vector3D& v)
{
 return Vector3D(v.get_x_value()*val, v.get_y_value()*val, v.get_z_value()*val);
}

Vector3D operator* (const Vector3D& v, double val)
{
 return Vector3D(v.get_x_value()*val, v.get_y_value()*val, v.get_z_value()*val);
}

ostream& operator<<(ostream& os, const Vector3D& atom)
{
 os << "(" << atom.get_x_value() << "," << atom.get_y_value() << "," << atom.get_z_value() << ")";
 return os;
}

/* -----------------------------Matrix---------------------------- */
Matrix::Matrix(int rows, Vector3D *row)
{
 this->row = new Vector3D[rows];
 for(int i=0; i<rows; i++)
 {
  this->row[i] = row[i];
 }
 this->rows = rows;
}

Matrix::~Matrix()
{
 delete[] row;
}

Matrix* Matrix::rotation_matrix(Vector3D& axis, double rad)
{
 /* axis components */
 double ax = axis.get_x_value();
 double ay = axis.get_y_value();
 double az = axis.get_z_value();
 Vector3D *v = new Vector3D[3];
 
 /* build rotation matrix */
 v[0] = Vector3D( cos(rad)+(1-cos(rad))*ax*ax, (1-cos(rad))*ax*ay-sin(rad)*az, (1-cos(rad))*ax*az+sin(rad)*ay );
 v[1] = Vector3D( (1-cos(rad))*ax*ay+sin(rad)*az, cos(rad)+(1-cos(rad))*ay*ay, (1-cos(rad))*ay*az-sin(rad)*ax );
 v[2] = Vector3D( (1-cos(rad))*az*ax-sin(rad)*ay, (1-cos(rad))*az*ay+sin(rad)*ax, cos(rad)+(1-cos(rad))*az*az );
 Matrix *m = new Matrix(3, v);

 /* clean up */
 delete[] v;
 
 return m;
}

/* operators */
Vector3D& Matrix::operator[] (unsigned int row)
{
 return this->row[row];
}
 
Vector3D Matrix::operator* (Vector3D& v)
{
 double comp[3];
 for(int i=0; i<rows; i++)
 {
  comp[i] = row[i]*v;
 }
 return Vector3D(comp[0],comp[1],comp[2]);
}
