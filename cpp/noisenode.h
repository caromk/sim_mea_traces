#ifndef ____noisenode__
#define ____noisenode__

#include <iostream>
#include "location.h"
#include <vector>
#include "neuron.h"
#include "simparams.h"
#include "Eigen/Dense"
#include <random>
#include "randomutil.h"
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class noisenode {
  
public:
  // constructors
  noisenode(location* locP,std::default_random_engine* generatorP,params* sim_paramsP);
  noisenode(location* locP,VectorXd trace,params* sim_paramsP);
  
  void init_noise();
  
  VectorXd noise_trace(location* det_locP);
  
  double distance_snr_multiplier(location* det_locP);
  double distance_to_detector(location* det_locP);
  double distance_to_detector_over_probe(location* det_locP);
  
  // variable access
  VectorXd* return_trace();
  location* return_location();
  
private:
  location* locP;
  VectorXd trace;
  params* sim_paramsP;
  default_random_engine* generatorP;
  
};

#endif /* defined(____noisenode__) */
