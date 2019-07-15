#ifndef ____detector__
#define ____detector__

#include <iostream>
#include "location.h"
#include <vector>
#include "neuron.h"
#include "simparams.h"
#include "Eigen/Dense"
#include <random>
#include "randomutil.h"
#include "noisenode.h"
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class detector {
  
public:
  // constructor
  detector(location* locP,vector<neuron*> *neuronVP,vector<noisenode*> *noisenodeVP,std::default_random_engine* generatorP,params* sim_paramsP);
  
  // calculate all the parts
  void calc_spike_train();
  void calc_johnson_noise();
  void calc_shared_noise();
  void calc_trace();
  
  // variable access
  VectorXd* return_trace();
  location* return_location();
  
private:
  location* locP;
  VectorXd johnson_noise;
  VectorXd shared_noise;
  VectorXd spike_train;
  VectorXd trace;
  vector<neuron*> *neuronVP;
  vector<noisenode*> *noisenodeVP;
  params* sim_paramsP;
  int n_neurons;
  default_random_engine* generatorP;
  
  int n_noisenode;
};

#endif /* defined(____detector__) */
