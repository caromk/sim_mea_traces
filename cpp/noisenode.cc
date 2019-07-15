#include "noisenode.h"
using namespace std;

noisenode::noisenode(location* locP,std::default_random_engine* generatorP,params* sim_paramsP)
{
  this->locP = locP;
  this->sim_paramsP = sim_paramsP;
  this->generatorP = generatorP;

  // make the noise
  init_noise();
}

noisenode::noisenode(location* locP,VectorXd trace,params* sim_paramsP)
{
  this->locP = locP;
  this->trace = trace;
  this->sim_paramsP = sim_paramsP;
}

void noisenode::init_noise()
{
  trace = rand_matrix_norm_amp_limit(sim_paramsP->n_pts,1,0,1,3.5,generatorP);
}

VectorXd noisenode::noise_trace(location* det_locP)
{
  return distance_snr_multiplier(det_locP)*trace;
}

double noisenode::distance_to_detector(location* det_locP) {
  return distance(locP,det_locP);
}

VectorXd* noisenode::return_trace()
{
  return &trace;
}

// NOTE - this function is duplicated in noisenode.cc
double noisenode::distance_to_detector_over_probe(location* det_locP) {
  
  // new assignment of variables for readability / sanity
  
  // electrode coordinates
  double xd = det_locP->loc(0);
  double yd = det_locP->loc(1);
  
  // probe coordinates
  double xp = -sim_paramsP->x_edge;
  // yp - unknown, minimized distnace over yp to find the equation below
  double zp = det_locP->loc(2); // pad z coordinate is the same as probe
  
  // neuron coordinates
  double xs = locP->loc(0);
  double ys = locP->loc(1);
  double zs = locP->loc(2);
  
  /**
   minimum distance calculated by minimizing over the y crossing point on the probe - did manually to a point, then moved over to mathematica for accuracy, equation copied from mathematica, mathematica file titled probe_min_dist_final
   **/
  
  // calculate for both negative and positive edges of the probe
  double distances_xp_pos = sqrt(pow(xd-xp,2)+pow(pow(xd-xp,2)*(yd-ys)+sqrt(pow(xd-xp,2)*pow(yd-ys,2)*(pow(xp-xs,2)+pow(zp-zs,2))),2)/pow(-(xd-xs)*(xd-2*xp+xs)+pow(zp-zs,2),2))+sqrt(pow(xp-xs,2)+pow((yd-ys)*(pow(xp-xs,2)+pow(zp-zs,2))+sqrt(pow(xd-xp,2)*pow(yd-ys,2)*(pow(xp-xs,2)+pow(zp-zs,2))),2)/pow(-(xd-xs)*(xd-2*xp+xs)+pow(zp-zs,2),2)+pow(zp-zs,2));
  double distances_xp_neg = sqrt(pow(xd+xp,2)+pow(pow(xd+xp,2)*(yd-ys)+sqrt(pow(xd+xp,2)*pow(yd-ys,2)*(pow(-xp-xs,2)+pow(zp-zs,2))),2)/pow(-(xd-xs)*(xd-2*-xp+xs)+pow(zp-zs,2),2))+sqrt(pow(-xp-xs,2)+pow((yd-ys)*(pow(-xp-xs,2)+pow(zp-zs,2))+sqrt(pow(xd+xp,2)*pow(yd-ys,2)*(pow(-xp-xs,2)+pow(zp-zs,2))),2)/pow(-(xd-xs)*(xd-2*-xp+xs)+pow(zp-zs,2),2)+pow(zp-zs,2));
  
  // take the shorter distance
  double distance = min(distances_xp_pos,distances_xp_neg);
  
  return distance;
}


double noisenode::distance_snr_multiplier(location* det_locP) {
  double dist;
  if (locP->loc(2) >= 0) {
    dist = distance_to_detector(det_locP);
  }
  else {
    dist = distance_to_detector_over_probe(det_locP);
  }
  
  double mult = 1/(sim_paramsP->a_falloff*pow(dist,2)+sim_paramsP->b_falloff*dist+sim_paramsP->c_falloff);
  
  return mult;
}

location* noisenode::return_location()
{
  return locP;
}