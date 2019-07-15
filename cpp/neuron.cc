#include "neuron.h"
using namespace std;

neuron::neuron(double poisson_rate, location *locP, direction *dirP,waveform* waveP,std::default_random_engine* generatorP,params* sim_paramsP) {
  this->poisson_rate = poisson_rate;
  this->locP = locP;
  this->dirP = dirP;
  this->generatorP = generatorP;
  this->waveP = waveP;
  this->sim_paramsP = sim_paramsP;
  this->add_n_pts = ceil(sim_paramsP->sampling_rate * 2.5 * (sim_paramsP->spk_length+sim_paramsP->refractory_period));
  
  create_spike_times();
  draw_spike_index();
}

neuron::neuron(vector<int> spike_times,MatrixXi index_spikes_drawn,location *locP, direction *dirP,waveform* waveP,std::default_random_engine* generatorP,params* sim_paramsP) {
  
  this->spike_times = spike_times;
  this->index_spikes_drawn = index_spikes_drawn;
  this->locP = locP;
  this->dirP = dirP;
  this->waveP = waveP;
  this->generatorP = generatorP;
  this->sim_paramsP = sim_paramsP;
  this->add_n_pts = ceil(sim_paramsP->sampling_rate * 2.5 * (sim_paramsP->spk_length+sim_paramsP->refractory_period));
  
  n_spikes = spike_times.size();
}

void neuron::draw_spike_index(){
  index_spikes_drawn = rand_matrix_uniform_int(n_spikes,1,0,sim_paramsP->n_draws_per_neuron-1,generatorP);
}

VectorXd neuron::spike_train(location* det_locP) {
  // initialize output trace
  VectorXd spike_trace(sim_paramsP->n_pts);
  spike_trace = VectorXd::Zero(sim_paramsP->n_pts);
  // calculate distances
  double dist_along_spk_nonlin=0;
  if(in_vector_direction(det_locP,locP,dirP))
    dist_along_spk_nonlin= normal_distance_to_detector(det_locP);
  // make waveforms based on distances
  MatrixXd waves = waveP->waves_from_distance(&index_spikes_drawn,dist_along_spk_nonlin);
  //cout <<  endl << "waves = " << waves << endl;
  // values to calculate amplitude
  double dist_mult = distance_snr_multiplier(det_locP);
  std::uniform_real_distribution<double> distribution(sim_paramsP->train_mult_rand_min,sim_paramsP->train_mult_rand_max);
  double rand_mult = distribution(*generatorP);
  // place waveforms at trace
  for(int i_spk=0;i_spk<n_spikes;i_spk++)
  {
    // starts before begining of trace
    if (spike_times[i_spk]<0) {
      if (-spike_times[i_spk] < sim_paramsP->n_pts_per_spike )
      {
        int curr_n_pts_per_spike = sim_paramsP->n_pts_per_spike+spike_times[i_spk];
        //std::cout << spike_times[i_spk] << " " << curr_n_pts_per_spike << endl;
        spike_trace.segment(0,curr_n_pts_per_spike)= dist_mult*rand_mult*waves.row(i_spk).segment(-spike_times[i_spk],curr_n_pts_per_spike);
      }
    }
    // fully within the trace
    else if (spike_times[i_spk]+sim_paramsP->n_pts_per_spike-1 < sim_paramsP->n_pts) {
      spike_trace.segment(spike_times[i_spk],sim_paramsP->n_pts_per_spike) = dist_mult*rand_mult*waves.row(i_spk);
    }
    // end beyond end of trace
    else {
      int curr_n_pts_per_spike = sim_paramsP->n_pts - spike_times[i_spk];
      spike_trace.segment(spike_times[i_spk],curr_n_pts_per_spike) = dist_mult*rand_mult*waves.row(i_spk).segment(0,curr_n_pts_per_spike);
    }
  }
  return spike_trace;
}

double neuron::distance_to_detector(location* det_locP) {
  return distance(locP,det_locP);
}

double neuron::normal_distance_to_detector(location* det_locP) {
  location normal_cross_loc = location_normal_to_vector(det_locP,locP,dirP);
  return distance(&normal_cross_loc,det_locP);
}

// NOTE - this function is duplicated in noisenode.cc
double neuron::distance_to_detector_over_probe(location* det_locP) {
  
  // new assignment of variables for readability / sanity
    
  // electrode coordinates
  double xd = det_locP->loc(0);
  double yd = det_locP->loc(1);
  
  // probe coordinates
  double xp = -sim_paramsP->x_edge;
  // yp - unknown, minimized distance over yp to find the equation below
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

/**
  distance_snr_multiplier
    for a given location (x,y,z), determine the amplitude multiplier based on the location of this neuron
    - distance: for neurons on the front side of the probe (y>=0), distance between neuron and site is shortest path between two points,
      for neurons on the backside of the probe (y<0), the distance is the shortest from the site over the area of the probe to the neuron
    - multiplier: for distance <= radius_lin_regime, linear fall off, past that 1/d^2 fall off
**/
double neuron::distance_snr_multiplier(location* det_locP) {
  double dist;
  // determine distance between site and probe
  if (locP->loc(2) >= 0) {
    dist = distance_to_detector(det_locP);
  }
  else {
    dist = distance_to_detector_over_probe(det_locP);
  }

  // determine multiplier based on that distance
  double mult = 1/(sim_paramsP->a_falloff*pow(dist,2)+sim_paramsP->b_falloff*dist+sim_paramsP->c_falloff);
  
  return mult;
}

void neuron::create_spike_times() {
  
  // init calcs
  const double lambda = poisson_rate / sim_paramsP->sampling_rate; // for this sampling rate
  const double refractory_period_at_sr = sim_paramsP->refractory_period * sim_paramsP->sampling_rate;
  // sr = sampling rate 
  
  // initialize random number generator for lambda
  std::exponential_distribution<double> distribution(lambda);
  
  double cumm_isi = -add_n_pts;
  int i_isi = 0;
  while (cumm_isi < sim_paramsP->n_pts) {
    double curr_isi = distribution(*generatorP);
    //std::cout << "curr_isi: " << curr_isi << std::endl;
    // must be longer than the refractory period
    if (curr_isi>refractory_period_at_sr)
    {
      // running total of isis
      cumm_isi += curr_isi;
      // see if have run past total time you want
      if(cumm_isi < sim_paramsP->n_pts)
      {
        spike_times.push_back(floor(cumm_isi));
        //std::cout << "cumm_isi: "  << floor(cumm_isi) << std::endl;
        i_isi++;
      }
    }
  }
  
  n_spikes = i_isi;
//  std::cout << spike_times.size() << std::endl;
  
}

vector<int> *neuron::return_spk_tms()
{
  return &spike_times;
}

location* neuron::return_location()
{
  return locP;
}

direction* neuron::return_direction()
{
  return dirP;
}

MatrixXi* neuron::return_index_spikes_drawn()
{
  return &index_spikes_drawn;
}

vector<MatrixXd>* neuron::return_jiggled_turn_coord()
{
  return waveP->return_jiggled_turn_coord();
}