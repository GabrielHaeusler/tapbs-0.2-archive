#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

#define PI 3.1415926535897932384626433832795028841971

using namespace std;

/**
 * @brief Vector in 3D cartesian space.
 * Provides elementary vector operations
 * and implements an arbitrary less-than operator.
 */

class Vector3D
{
private:
 double x;                            /** x component */
 double y;                            /** y component */
 double z;                            /** z component */

public:

 /**
  * @brief Standard constructor
  * @author Gernot Kieseritzky
  */

 Vector3D();

 /**
  * @brief Constructor with intial values.
  * @author Gernot Kieseritzky
  */

 Vector3D(double x, double y, double z);

 /**
  * @brief Get x-component.
  * @author Gernot Kieseritzky
  * @return Double value that represents the x-component.
  */

 double get_x_value() const {return this->x;}

 /**
  * @brief Get y-component.
  * @author Gernot Kieseritzky
  * @return Double value that represents the y-component.
  */

 double get_y_value() const {return this->y;}

 /**
  * @brief Get z-component.
  * @author Gernot Kieseritzky
  * @return Double value that represents the z-component.
  */

 double get_z_value() const {return this->z;}

 /**
  * @brief Get all components at once.
  * @author Gernot Kieseritzky
  */

 void get_value(double array[]) const
 { 
  array[0] = this->get_x_value();
  array[1] = this->get_y_value();
  array[2] = this->get_z_value();
 }

 /**
  * @brief Calculate absolute value.
  * @author Gernot Kieseritzky
  * @return Double value that represents the absolute value.
  */

 double abs();

 /**
  * @brief Set all three components.
  * @author Gernot Kieseritzky
  */

 void set(double x, double y, double z);

 /**
  * @brief Rotate vector around axis by 'deg' degree.
  * @author Gernot Kieseritzky
  */

 void rotate(Vector3D& axis, double deg);

 /**
  * @brief Translate vector along axis.
  * @author Gernot Kieseritzky
  */

 void translate(Vector3D& axis);

 /**
  * @brief Calculate angle between to values.
  * @author Gernot Kieseritzky
  * @return Double values that represents the angle in degrees.
  */

 double angle(Vector3D& v);

 /**
  * @brief Convert degrees into radians.
  * @author Gernot Kieseritzky
  * @return Double values that represents the angle in radians.
  */

 static double deg2rad(double deg);

 /**
  * @brief Convert radians into degrees.
  * @author Gernot Kieseritzky
  * @return Double values that represents the angle in degrees.
  */

 static double rad2deg(double rad);

 /* operators */

 /**
  * @brief Equals operator
  * Two vectors are equal iff all components are equal.
  * @author Gernot Kieseritzky
  * @return Boolean.
  */

 bool operator==(const Vector3D& v);

 /**
  * @brief Not-Equals operator
  * Two vectors are not equal if one of the components differ
  * @author Gernot Kieseritzky
  * @return Boolean.
  */

 bool operator!=(const Vector3D& v) {return !(*this==v);}

 /**
  * @brief Less than
  * Vector a is smaller than vector b if
  * x-component or, if x>=, y-component or, if x>= and y >=, z-component is smaller.
  * @author Gernot Kieseritzky
  * @return Boolean.
  */

 bool operator<(const Vector3D& v) const;

 /**
  * @brief Add two vectors.
  * Simply adds components.
  * @author Gernot Kieseritzky
  * @return Result vector.
  */

 Vector3D operator+ (const Vector3D& b) const;

 /**
  * @brief Subtract two vectors.
  * @author Gernot Kieseritzky
  * @return Result vector.
  */

 Vector3D operator- (const Vector3D& b) const;

 /**
  * @brief Scalar product.
  * @author Gernot Kieseritzky
  * @return Result vector.
  */

 double operator* (const Vector3D& b) const;

 /**
  * @brief Vector product.
  * @author Gernot Kieseritzky
  * @return Result vector.
  */

 Vector3D cross(const Vector3D& b) const;
};

/**
 * @brief Product of scalar and vector. (a*V)
 * @author Gernot Kieseritzky
 * @return Result vector.
 */

Vector3D operator* (double val, const Vector3D& v);

/**
 * @brief Product of scalar and vector. (V*a)
 * Should be commutative.
 * @author Gernot Kieseritzky
 * @return Result vector.
 */

Vector3D operator* (const Vector3D& v, double val);

/**
 * @brief Overload output stream operator.
 * (x, y, z)
 * @author Gernot Kieseritzky
 * @return Modified output stream.
 */

ostream& operator<<(ostream& os, const Vector3D& atom);

/**
 * @brief Basic matrix.
 * Matrix is constructed from row vectors.
 * @author Gernot Kieseritzky
 */

class Matrix
{
private:
 Vector3D *row;                   /** Row vectors */
 int rows;                        /** number of rows */
public:

/**
 * @brief Standard constructor.
 * @author Gernot Kieseritzky
 */

 Matrix(int rows, Vector3D *row);

/**
 * @brief Calculate a rotation matrix.
 * @author Gernot Kieseritzky
 * @return Pointer to new rotation 3x3 matrix.
 */

 static Matrix* rotation_matrix(Vector3D& axis, double rad);

 /* operators */

/**
 * @brief Index operator.
 * @author Gernot Kieseritzky
 * @return Reference to row vector.
 */

 Vector3D& operator[] (unsigned int row);

/**
 * @brief Matrix vector product.
 * @author Gernot Kieseritzky
 * @return Result vector.
 */

 Vector3D operator* (Vector3D& v);

/**
 * @brief Standard destructor.
 * @author Gernot Kieseritzky
 */

 ~Matrix();
};
