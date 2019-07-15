#ifndef ____waveform__
#define ____waveform__

#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include "Eigen/Dense"
#include "randomutil.h"
#include "simparams.h"
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;

class waveform {
  
public:
  waveform(MatrixXd template_turn_coord,double ratio_stdev_jiggle_x,double ratio_stdev_jiggle_y,int n_draw,std::default_random_engine* generatorP,params* sim_paramsP);
  waveform(vector<MatrixXd> jiggled_turn_coord,params* sim_paramsP);
  
  // setup functions
  void create_n_draws();
  MatrixXd jiggle_turns();
  void set_time_values();
  void set_x_values();
  
  // called from neuron object
  MatrixXd waves_from_distance(MatrixXi *index_spikes_drawnP,double dist_from_neuron);
  VectorXd calc_wave(VectorXd curr_turns_x,VectorXd curr_turns_y);
  VectorXd fit_cubic_two_turns(double x1,double x2,double y1,double y2);
  MatrixXd fit_cubics_n_turns(VectorXd turn_coord_x,VectorXd turn_coord_y);
  
  // variable access functions
  vector<MatrixXd>* return_jiggled_turn_coord();

private:

  MatrixXd template_turn_coord;
  double ratio_stdev_jiggle_x;
  double ratio_stdev_jiggle_y;
  double stdev_jiggle_x;
  double stdev_jiggle_y;
  int n_draw;
  std::default_random_engine* generatorP;
  
  int n_turns;
  vector<MatrixXd> jiggled_turn_coord;
  
  params* sim_paramsP;
  
  VectorXd x_values;
};

#endif /* defined(____waveform__) */
