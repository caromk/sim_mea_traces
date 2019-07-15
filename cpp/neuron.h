#ifndef __neuron__
#define __neuron__

#include <vector>
#include <random>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <iostream>
#include "Eigen/Dense"
#include "location.h"
#include "direction.h"
#include "waveform.h"
#include "geometryutil.h"
#include "simparams.h"
using namespace std;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::MatrixXi;

class neuron {

public:

  /**
   * default constructor
   */

  neuron(double poisson_rate,location *locP, direction *dirP,waveform* waveP,std::default_random_engine* generatorP,params* sim_paramsP);
  
  neuron(vector<int> spike_times,MatrixXi index_spikes_drawn,location *locP, direction *dirP,waveform* waveP,std::default_random_engine* generatorP,params* sim_paramsP);
  
  // setup functions
  void create_spike_times();
  void draw_spike_index();
  
  // called by detector
  double distance_to_detector(location* det_locP);
  double distance_to_detector_over_probe(location* det_locP);
  double normal_distance_to_detector(location* det_locP);
  double distance_snr_multiplier(location* det_locP);
  VectorXd spike_train(location* det_locP);
  
  // variable access
  vector<int>* return_spk_tms();
  location* return_location();
  direction* return_direction();
  MatrixXi* return_index_spikes_drawn();
  vector<MatrixXd>* return_jiggled_turn_coord();

private:
  double poisson_rate;
  int add_n_pts;
  vector<int> spike_times;
  int n_spikes;
  MatrixXi index_spikes_drawn;
  params* sim_paramsP;
  
  location *locP;
  direction *dirP;
  std::default_random_engine* generatorP;
  waveform* waveP;
  
};

#endif
