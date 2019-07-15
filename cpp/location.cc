#include "location.h"

location::location(Vector3d loc)
{
  this->loc = loc;
}

location::location(double x, double y, double z) {
  loc << x,y,z;
  
}
