#ifndef ____direction__
#define ____direction__

#include <iostream>
#include "Eigen/Dense"
using namespace std;
using Eigen::Vector3d;

class direction {
  
public:
  
  direction(Vector3d dir);
  direction(double x, double y, double z);

  Vector3d dir;

};

#endif /* defined(____direction__) */
