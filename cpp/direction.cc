#include "direction.h"

direction::direction(Vector3d dir)
{
  this->dir = dir;
}

direction::direction(double x, double y, double z) {
  dir << x,y,z;
}