
#include "detector.h"
using namespace std;

detector::detector(location* locP,vector<neuron*> *neuronVP,vector<noisenode*> *noisenodeVP,std::default_random_engine* generatorP,params* sim_paramsP)
{
  this->locP = locP;
  this->neuronVP = neuronVP;
  this->noisenodeVP = noisenodeVP;
  this->sim_paramsP = sim_paramsP;
  this->n_neurons = neuronVP->size();
  this->generatorP = generatorP;
  n_noisenode = noisenodeVP->size();
  
  // get neural portion of signal
  calc_spike_train();
  // todo: get shared noise portion of signal
  calc_shared_noise();
  // todo: get johnson noise portion of signal
  calc_johnson_noise();
  // todo: add it all up
  calc_trace();
  //cout << trace << endl;
}

VectorXd* detector::return_trace()
{
  return &trace;
}

void detector::calc_spike_train()
{
  spike_train = VectorXd::Zero(sim_paramsP->n_pts);
  
  for (int i_neuron=0;i_neuron<this->n_neurons;i_neuron++)
  {
    spike_train += neuronVP->at(i_neuron)->spike_train(locP);
  }
  
}

void detector::calc_johnson_noise()
{
  johnson_noise = rand_matrix_norm_amp_limit(sim_paramsP->n_pts,1,0,sim_paramsP->scale_johnson_noise,3.5,generatorP);

}

void detector::calc_shared_noise()
{
  shared_noise = VectorXd::Zero(sim_paramsP->n_pts);
  
  for (int i_noisenode=0;i_noisenode<n_noisenode;i_noisenode++)
  {
    shared_noise += noisenodeVP->at(i_noisenode)->noise_trace(locP);
  }
  
  shared_noise *= sim_paramsP->rms_gauss_noise / stdev(shared_noise);
}

void detector::calc_trace()
{
  trace = spike_train + johnson_noise + shared_noise;
}

location* detector::return_location()
{
  return locP;
}