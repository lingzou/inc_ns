//
//  Point.h
//
//
//  Created by Ling Zou on 5/27/13.
//
//

#ifndef _Point_h
#define _Point_h

#include <cmath>
#include <stdlib.h>

class Vec3d {
public:
  Vec3d() {}
  ~Vec3d() {}

  Vec3d(double xx, double yy, double zz) :
    _x(xx), _y(yy), _z(zz)
  {}

  void print() { std::cout << "(" << _x << ", " << _y << ", " << _z << ")" << std::endl; }
  const double x() const { return _x; }
  const double y() const { return _y; }
  const double z() const { return _z; }
  double & x() { return _x; }
  double & y() { return _y; }
  double & z() { return _z; }

  inline Vec3d operator+ (const Vec3d & v) const {
    return Vec3d(x() + v.x(), y() + v.y(), z() + v.z()); }

  inline Vec3d operator- (const Vec3d & v) const {
    return Vec3d(x() - v.x(), y() - v.y(), z() - v.z()); }

  inline double operator * (const Vec3d & v) const {
    return x() * v.x() + y() * v.y() + z() * v.z(); }

  inline Vec3d operator / (const double & s) const
  {
    if(s == 0)
    {
      std::cerr << "File: " << __FILE__ << std::endl;
      std::cerr << "Function: " << __FUNCTION__ << std::endl;
      std::cerr << "Line: " << __LINE__ << std::endl;
      std::cerr << "ERROR: Vec3d divided by zero.\n";
      exit(1);
    }
    else
      return Vec3d(x()/s, y()/s, z()/s);
  }

  inline Vec3d cross(const Vec3d & v) {
    return Vec3d(y() * v.z() - z() * v.y(),
                 z() * v.x() - x() * v.z(),
                 x() * v.y() - y() * v.x());
  }

  double norm() {
    return std::sqrt(_x * _x + _y * _y + _z * _z); }

  Vec3d unitVector()
  {
    if(norm() == 0)
      return Vec3d(0, 0, 0);
    else
      return (*this) / norm();
  }

protected:
  double _x, _y, _z;
};

class Point
{
public:
  Point()
  {}

  Point(double xx, double yy, double zz) :
    _x(xx),
    _y(yy),
    _z(zz)
  {
  }

  ~Point()
  {}

  double x() const { return _x; }
  double y() const { return _y; }
  double z() const { return _z; }

  double & x() { return _x; }
  double & y() { return _y; }
  double & z() { return _z; }

  inline
  Point operator+ (const Point & p) const
  {
    return Point(x() + p.x(), y() + p.y(), z() + p.z());
  }

  inline
  Vec3d operator- (const Point & p) const
  {
    return Vec3d(x() - p.x(), y() - p.y(), z() - p.z());
  }

protected:
  double _x, _y, _z;
};

#endif
