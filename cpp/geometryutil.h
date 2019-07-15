#ifndef __geometryutil__
#define __geometryutil__

#include "location.h"
#include "direction.h"
#include <cmath>
#include "Eigen/Dense"
#include <string>
#include <fstream>
#include <limits>
#define _USE_MATH_DEFINES
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::numeric_limits;

// code snippet from http://stackoverflow.com/questions/18400596/how-can-a-eigen-matrix-be-written-to-file-in-csv-format

inline void save_matrixxd_csv(MatrixXd* m,string filename)
{
  // construct and initiate
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  
  for(int  i_row = 0; i_row < m->rows(); i_row++){
    for(int i_col = 0; i_col < m->cols(); i_col++){
      ofs.precision(numeric_limits<double>::digits10 + 1);
      ofs.setf( std::ios::fixed, std:: ios::floatfield );
      if(i_col+1 == m->cols()){
        ofs << (*m)(i_row,i_col);
      }else{
        ofs << (*m)(i_row,i_col) << ',';
      }
    }
    ofs<<'\n';
  }
  
  ofs.close();
}

inline double norm(VectorXd v) {
   return sqrt(v.array().pow(2).sum());
}

inline double stdev(VectorXd v) {
  double mean = v.array().sum()/v.size();
  double mean_square = v.array().pow(2).sum()/v.size();
  
  return sqrt(mean_square-pow(mean,2));
}

inline double distance(location *loc1P, location *loc2P) {  
  return norm(loc1P->loc.array() - loc2P->loc.array());
}

/*
Given a location, a direction/vector and the locations at which that vector begins, finds the point at which the normal vector would intersect the original vector
 */

inline location location_normal_to_vector(location *ini,location *v_base,direction *v)
{
  return location(v_base->loc+(v->dir.dot(ini->loc-v_base->loc)/v->dir.array().pow(2).sum())*v->dir);
}

/* from the interwebz - acos can run into rounding error issues, resulting in nan being returned with the input is very near +/-1
 http://stackoverflow.com/questions/8489792/is-it-legal-to-take-acos-of-1-0f-or-1-0f
 */

inline double safe_acos (double x)
{
  if (x < -1.0) x = -1.0 ;
  else if (x > 1.0) x = 1.0 ;
  return acos (x) ;
}

/*
 Given a location, a direction/vector and the locations at which that vector begins, determines if the location is in the same "side" of the space as the vector (e.g., if you made a plane perpendicular to the base of the vector, would the location be on the same side of that plane as the vector's direction?
 */
inline bool in_vector_direction(location *ini,location *v_base,direction *v)
{  
  return safe_acos(v->dir.dot(ini->loc-v_base->loc)/(norm(v->dir)*norm(ini->loc-v_base->loc))) < M_PI/2;
}

#endif