#ifndef ____location__
#define ____location__

#include "Eigen/Dense"
using namespace std;
using Eigen::Vector3d;

class location {
  
public:
  
  location(Vector3d loc);
  location(double x,double y,double z);
  
  Vector3d loc;
  
};

#endif /* defined(____location__) */
