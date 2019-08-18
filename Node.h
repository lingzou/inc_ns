//
//  Node.h
//
//
//  Created by Ling Zou on 5/27/13.
//
//

#ifndef _Node_h
#define _Node_h

#include "Point.h"

class Node
{
public:
  // FIXME: consturct without id given?
  Node(Point p) :
    _point(p),
    _id(-1) // an invalid id to initialize
  {
  }

  Node(double x, double y, double z)
  {
    _point.x() = x;
    _point.y() = y;
    _point.z() = z;
  }

  const unsigned int id() const { return _id; }
  unsigned int& id() { return _id; }

  inline const Point & point() const { return _point; }

  const double x() const { return _point.x(); }
  const double y() const { return _point.y(); }
  const double z() const { return _point.z(); }

protected:
  Point _point;
  unsigned int _id;
};


#endif
