#ifndef _simparams_
#define _simparams_

#include <cmath>

using namespace std;

struct params {
  // general
  int sampling_rate = 30000;
  double length_in_sec = 1;
  int n_pts;
  
  // electrode
  int n_electrode_rows = 2; //12;
  int n_electrode_columns = 2;
  double dist_betw_electrodes = 10; // microns
  bool rand_detectors = 0;
  double x_edge = 0;
  double y_edge = 0;
  int n_detectors;
  
  // waveforms
  int n_draws_per_neuron = 50;
  double ratio_stdev_jiggle_x_neuron = 0.096;
  double ratio_stdev_jiggle_y_neuron = 0.25;
  double ratio_stdev_jiggle_x_draw = 0.025;
  double ratio_stdev_jiggle_y_draw = 0.025;
  double dist_to_fin = 60;
  double dist_to_zero = 10000;
  double spk_length = 1.5/1000;
  double refractory_period = 0.005;
  int n_pts_per_spike;
  
  // amplitude falloff
  // y = 1/(ax^2 + bx + c)
  double a_falloff = 4.083e-06;
  double b_falloff = 5.166e-05;
  double c_falloff = 4.000e-03;
//  double max_amp_lin_regime = 200;
//  double min_amp_lin_regime = 150;
//  double radius_lin_regime = 30;
//  double x_value_in_inv_regime = 100;
//  double y_value_in_inv_regime = 20;
  double train_mult_rand_min = 0.7;
  double train_mult_rand_max = 1.3;
  
  // noise
  double rms_gauss_noise = 15;
  double dist_betw_gauss_noise_nodes = 40;
  double scale_johnson_noise = 10;
  
  // space layout
  double density_neurons = 60*pow(10,-6);
  double volume_radius = 350; //microns
  double volume_cylinder_height = 800;
  
  // firing rate
  bool log_normal_firing_distrib = 0;
  double poisson_rate_lognorm_mean = 0;
  double poisson_rate_lognorm_var = 0;
  double poisson_rate_lognorm_mu = 0.92;
  double poisson_rate_lognorm_sigma = 1.18;
  double poisson_rate_max = 20;
  double poisson_rate_min = 0.5;
  double active_neuron_rate = 1;

  // rand (initialized in code)
  int rand_seed = time(NULL);
};

#endif
