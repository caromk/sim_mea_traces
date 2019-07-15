#include <vector>
#include <string>
#include <iostream> 
#include <random>
#include "neuron.h"
#include "location.h"
#include "direction.h"
#include "geometryutil.h"
#include "waveform.h"
#include "detector.h"
#include "Eigen/Dense"
#include "simparams.h"
#include <fstream>
#include <cmath>
#include "time.h"
#include "noisenode.h"
#include <sstream>
#include <stdlib.h> 
#include <unordered_map>
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

void prep_neurons(vector<neuron*> *neuronVP,waveform template_w,std::default_random_engine* generatorP,params* sim_paramsP);
void prep_detectors(vector<detector*> *detectorVP,vector<neuron*> *neuronVP,vector<noisenode*> *noisenodeVP,std::unordered_map<std::string,int> *detector_locP,string out_file_prefix,std::default_random_engine* generatorP,params* sim_paramsP);
void prep_noisenodes(vector<noisenode*> *noisenodeVP,std::default_random_engine* generatorP,params* sim_paramsP);

void save_detector_traces(string file_prefix,vector<detector*> *detectorVP);
void save_spike_times(string file_prefix,vector<neuron*> *neuronVP);
void save_detector_locations(string file_prefix,vector<detector*> *detectorVP);
void save_neuron_locations(string file_prefix,vector<neuron*> *neuronVP);
void save_neuron_directions(string file_prefix,vector<neuron*> *neuronVP);
void save_noisenode_locations(string file_prefix,vector<noisenode*> *noisenodeVP);
void save_noisenode_traces(string file_prefix,vector<noisenode*> *noisenodeVP);
void save_index_spikes_drawn(string file_prefix,vector<neuron*> *neuronVP);
void save_jiggled_turn_coord(string file_prefix,vector<neuron*> *neuronVP);

void waveform_from_file(vector<waveform*> *waveformVP,string file_prefix,params* sim_paramsP);
void neuron_from_file(vector<neuron*> *neuronVP,vector<waveform*> *waveformVP,string file_prefix,std::default_random_engine* generatorP,params* sim_params);
void noisenode_from_file(vector<noisenode*> *noisenodeVP,string file_prefix,params* sim_paramsP);

void read_detector_loc_file(string detector_loc_file,std::unordered_map<std::string,int> *detector_locP);

/**
Usage example, runs in reasonable time:
  ./simtraces -cols 2 -rows 2 -out "test"
**/

int main(int argc, char* argv[])
{  
  // time
  time_t rawtime;
  
  // 1
  time (&rawtime);
  cout << ctime (&rawtime) << "starting" << endl;
  
  // seed random number generator
  std::default_random_engine generator;
  
  // Prepare parameters
  
  // detector location hash table
  std::unordered_map<std::string,int> detector_loc;
  
  // shared parameters
  params sim_params;
  
  // check arguments
  // prefix for saving files
  string in_file_prefix;
  string out_file_prefix;
  string detector_loc_file;
  
  // input argument check loop
  for (int i = 1; i < argc; i+=2) {
    if (i + 1 != argc) {
      if (string(argv[i]) == "-in") {
        in_file_prefix = argv[i + 1];
      } else if (string(argv[i]) == "-out") {
        out_file_prefix = argv[i + 1];
      } else if (string(argv[i]) == "-rows") {
        sim_params.n_electrode_rows = atoi(argv[i + 1]);
      } else if (string(argv[i]) == "-cols") {
        sim_params.n_electrode_columns = atoi(argv[i + 1]);
      } else if (string(argv[i]) == "-pitch") {
        sim_params.dist_betw_electrodes = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-length") {
        sim_params.length_in_sec = atoi(argv[i + 1]);
      } else if (string(argv[i]) == "-detector_loc") {
        detector_loc_file = argv[i + 1];
      } else if (string(argv[i]) == "-rms_gauss_noise") {
        sim_params.rms_gauss_noise = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-scale_johnson_noise") {
        sim_params.scale_johnson_noise = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-density_neurons") {
        sim_params.density_neurons = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-volume_radius") {
        sim_params.volume_radius = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-volume_cylinder_height") {
        sim_params.volume_cylinder_height = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-n_rand_detectors") {
        sim_params.n_detectors = atoi(argv[i + 1]);
        sim_params.rand_detectors = 1;
      } else if (string(argv[i]) == "-x_edge") {
        sim_params.x_edge = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-y_edge") {
        sim_params.y_edge = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-poisson_rate_max") {
        sim_params.poisson_rate_max = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-poisson_rate_min") {
        sim_params.poisson_rate_min = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-volume_radius") {
        sim_params.volume_radius = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-volume_cylinder_height") {
        sim_params.volume_cylinder_height = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-rand_seed") {
        sim_params.rand_seed = atoi(argv[i + 1]);
      } else if (string(argv[i]) == "-log_normal_firing_distrib") {
        sim_params.log_normal_firing_distrib = !!atoi(argv[i + 1]);
      } else if (string(argv[i]) == "-poisson_rate_lognorm_mean") {
        sim_params.poisson_rate_lognorm_mean = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-poisson_rate_lognorm_var") {
        sim_params.poisson_rate_lognorm_var = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-poisson_rate_lognorm_mu") {
        sim_params.poisson_rate_lognorm_mu = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-poisson_rate_lognorm_sigma") {
        sim_params.poisson_rate_lognorm_sigma = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-active_neuron_rate") {
        sim_params.active_neuron_rate = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-train_mult_rand_min") {
        sim_params.train_mult_rand_min = atof(argv[i + 1]);
      } else if (string(argv[i]) == "-train_mult_rand_max") {
        sim_params.train_mult_rand_max = atof(argv[i + 1]);
      } 
    }
  } 

  // generate random seed
  generator.seed (sim_params.rand_seed);

  // calc shared params
  sim_params.n_pts = ceil(sim_params.sampling_rate*sim_params.length_in_sec);
  sim_params.n_pts_per_spike = floor(sim_params.sampling_rate*sim_params.spk_length)+1;
  // need to calc x_edge and y_edge if not specified (only specified in the case of rand dist detectors)
  if (sim_params.x_edge == 0) {
    sim_params.x_edge = (sim_params.n_electrode_columns-1)*sim_params.dist_betw_electrodes;
    sim_params.y_edge = (sim_params.n_electrode_rows-1)*sim_params.dist_betw_electrodes;
    sim_params.n_detectors = sim_params.n_electrode_columns * sim_params.n_electrode_rows;
  }
  
  // Tell user what's running
  if (!sim_params.rand_detectors)
    cout << "Running " << sim_params.n_electrode_rows << "x" << sim_params.n_electrode_columns << " with " << sim_params.dist_betw_electrodes << " pitch for " << sim_params.length_in_sec << " second" << endl;
  else
    cout << "Running " << sim_params.n_detectors << " randomly distributed in a " << 2*sim_params.x_edge << " by " << 2*sim_params.y_edge << " area" << endl;
  cout << "Neural density: " << sim_params.density_neurons << endl;
  cout << "Firing rate between " <<  sim_params.poisson_rate_min << " and " << sim_params.poisson_rate_max << endl;
  cout << "Infile: " << in_file_prefix << endl << "Outfile: " << out_file_prefix << endl;

  // Create location hash table so as not to duplicate any existing locations
  if (!detector_loc_file.empty()) {
    read_detector_loc_file(detector_loc_file,&detector_loc);
  }
  
  // Prepare waveforms, neurons, and noise nodes
  vector<neuron*> neuronV;
  vector<noisenode*> noisenodeV;

  if (!in_file_prefix.empty()) {
    // get waveforms from file
    vector<waveform*> waveformV;
    waveform_from_file(&waveformV,in_file_prefix,&sim_params);

    time (&rawtime);
    cout << ctime (&rawtime) << "waveforms loaded" << endl;
    
    // get neurons from file
    neuron_from_file(&neuronV,&waveformV,in_file_prefix,&generator,&sim_params);
  
    time (&rawtime);
    cout << ctime (&rawtime) << neuronV.size() << " neurons loaded" << endl;
    
    // get noisenodes from file
    noisenode_from_file(&noisenodeV,in_file_prefix,&sim_params);
    
    time (&rawtime);
    cout << ctime (&rawtime) << noisenodeV.size() << " noise nodes loaded" << endl;
  }
  else {
    // master spike waveform template
    // waveform
    MatrixXd template_spk_turns(4,8);
    template_spk_turns <<
    0,2,4,7,12,16,18,20,
    0,0,0,1,-1.0/4.0,0,0,0,
    0,2,4,7,12,14,18,20,
    0,0,0,1,0,1.0/3.0,0,0;
    waveform template_w = waveform(template_spk_turns,sim_params.ratio_stdev_jiggle_x_neuron,sim_params.ratio_stdev_jiggle_y_neuron,0,&generator,&sim_params);
    
    // make binders of neurons
    prep_neurons(&neuronV,template_w,&generator,&sim_params);
  
    time (&rawtime);
    cout << ctime (&rawtime) << neuronV.size() << " neurons prepared" << endl;
  
    // prepare noise nodes
    prep_noisenodes(&noisenodeV,&generator,&sim_params);

    time (&rawtime);
    cout << ctime (&rawtime) << noisenodeV.size() << " noise nodes prepared" << endl;

    // save the base data
    save_spike_times(out_file_prefix,&neuronV);
    save_neuron_locations(out_file_prefix,&neuronV);
    save_neuron_directions(out_file_prefix,&neuronV);
    save_index_spikes_drawn(out_file_prefix,&neuronV);
    save_noisenode_locations(out_file_prefix,&noisenodeV);
    save_noisenode_traces(out_file_prefix,&noisenodeV);
    save_jiggled_turn_coord(out_file_prefix,&neuronV);
    time (&rawtime);
    cout << ctime (&rawtime) << "base data saved" << endl;
  }
  
  // calculate traces for each detector
  vector<detector*> detectorV;
  prep_detectors(&detectorV,&neuronV,&noisenodeV,&detector_loc,out_file_prefix,&generator,&sim_params);

  time (&rawtime);
  cout << ctime (&rawtime) << "detectors prepared" << endl;
  
  // save all the stuff
  save_detector_traces(out_file_prefix,&detectorV);
  
  time (&rawtime);
  cout << ctime (&rawtime) << "detector traces saved" << endl;
  
  return 0;

}

void save_noisenode_locations(string file_prefix,vector<noisenode*> *noisenodeVP)
{
  int n_noisenodes = noisenodeVP->size();
  string filename = file_prefix + "_noisenode_loc.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_nn = 0; i_nn < n_noisenodes; i_nn++)
  {
    location *locP = noisenodeVP->at(i_nn)->return_location();
    ofs << locP->loc(0) << "," << locP->loc(1) << "," << locP->loc(2) << '\n';
  }
  ofs.close();
}

void save_detector_locations(string file_prefix,vector<detector*> *detectorVP)
{
  int n_detectors = detectorVP->size();
  string filename = file_prefix + "_detector_loc.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_det = 0; i_det < n_detectors; i_det++)
  {
    location *locP = detectorVP->at(i_det)->return_location();
    ofs << locP->loc(0) << "," << locP->loc(1) << "," << locP->loc(2) << '\n';
  }
  ofs.close();
}

void save_detector_locations_from_locV(string file_prefix,vector<location*> *locationVP)
{
  int n_detectors = locationVP->size();
  string filename = file_prefix + "_detector_loc.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_det = 0; i_det < n_detectors; i_det++)
  {
    location *locP = locationVP->at(i_det);
    ofs << locP->loc(0) << "," << locP->loc(1) << "," << locP->loc(2) << '\n';
  }
  ofs.close();
}

void save_neuron_locations(string file_prefix,vector<neuron*> *neuronVP)
{
  int n_neurons = neuronVP->size();
  //cout << n_neurons << endl;
  string filename = file_prefix + "_neuron_loc.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    location *locP = neuronVP->at(i_neuron)->return_location();
    ofs << locP->loc(0) << "," << locP->loc(1) << "," << locP->loc(2) << '\n';
  }
  ofs.close();
}

void save_neuron_directions(string file_prefix,vector<neuron*> *neuronVP)
{
  int n_neurons = neuronVP->size();
  //cout << n_neurons << endl;
  string filename = file_prefix + "_neuron_dir.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    direction *dirP = neuronVP->at(i_neuron)->return_direction();
    ofs << dirP->dir(0) << "," << dirP->dir(1) << "," << dirP->dir(2) << '\n';
  }
  ofs.close();
}

void save_detector_traces(string file_prefix,vector<detector*> *detectorVP)
{
  int n_detectors = detectorVP->size();
  string filename = file_prefix + "_detector.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_det = 0; i_det < n_detectors; i_det++)
  {
    VectorXd *traceP = detectorVP->at(i_det)->return_trace();
    for (int i_tms=0;i_tms < traceP->size();i_tms++) {
      ofs.precision(numeric_limits<double>::digits10 + 1);
      ofs.setf( std::ios::fixed, std:: ios::floatfield );
      if (i_tms == traceP->size()-1)
        ofs << (*traceP)(i_tms);
      else
        ofs << (*traceP)(i_tms) << ",";
    }
    ofs << '\n';
  }
  ofs.close();
}

void save_noisenode_traces(string file_prefix,vector<noisenode*> *noisenodeVP)
{
  int n_noisenodes = noisenodeVP->size();
  string filename = file_prefix + "_noisenode.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_nn = 0; i_nn < n_noisenodes; i_nn++)
  {
    VectorXd *traceP = noisenodeVP->at(i_nn)->return_trace();
    for (int i_tms=0;i_tms < traceP->size();i_tms++) {
      ofs.precision(numeric_limits<double>::digits10 + 1);
      ofs.setf( std::ios::fixed, std:: ios::floatfield );
      if (i_tms == traceP->size()-1)
        ofs << (*traceP)(i_tms);
      else
        ofs << (*traceP)(i_tms) << ",";
    }
    ofs << '\n';
  }
  ofs.close();
}

void save_spike_times(string file_prefix,vector<neuron*> *neuronVP)
{
  int n_neurons = neuronVP->size();
  string filename = file_prefix + "_neuron.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    vector<int> *tmsP = neuronVP->at(i_neuron)->return_spk_tms();
    if (tmsP->size() == 0) {
      ofs << "-10000\n";
    }
    else {
      for (int i_tms=0;i_tms < tmsP->size();i_tms++) {
        if (i_tms == tmsP->size()-1)
          ofs << (*tmsP)[i_tms];
        else
          ofs << (*tmsP)[i_tms] << ",";
      }
      ofs << '\n';
    }
  }
  ofs.close();
}

void save_index_spikes_drawn(string file_prefix,vector<neuron*> *neuronVP)
{
  int n_neurons = neuronVP->size();
  string filename = file_prefix + "_index_spikes_drawn.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    MatrixXi* index_spikes_drawnP = neuronVP->at(i_neuron)->return_index_spikes_drawn();
    for (int i_spk=0; i_spk < index_spikes_drawnP->rows(); i_spk++) {
      if (i_spk == index_spikes_drawnP->rows()-1)
        ofs << (*index_spikes_drawnP)(i_spk);
      else
        ofs << (*index_spikes_drawnP)(i_spk) << ",";
    }
    if (index_spikes_drawnP->rows() == 0)
      ofs << -10000;
    ofs << '\n';
  }
  ofs.close();
}

void save_jiggled_turn_coord(string file_prefix,vector<neuron*> *neuronVP)
{
  int n_neurons = neuronVP->size();
  string filename = file_prefix + "_waveform_jiggled_turn_coord.txt";
  std::ofstream ofs;
  ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    vector<MatrixXd>* jiggled_turn_coordP = neuronVP->at(i_neuron)->return_jiggled_turn_coord();
    for(int i_draw=0;i_draw < jiggled_turn_coordP->size();i_draw++)
    {
      int n_rows = (*jiggled_turn_coordP)[i_draw].rows();
      int n_cols = (*jiggled_turn_coordP)[i_draw].cols();
      ofs << i_draw << ":" << n_rows << "," << n_cols << ":";
      for (int i_row=0; i_row < n_rows; i_row++) {
        for (int i_col=0; i_col < n_cols; i_col++) {
          ofs << (*jiggled_turn_coordP)[i_draw].coeff(i_row,i_col);
          if (i_row != n_rows-1 || i_col != n_cols-1)
            ofs << ",";
        }
      }
      ofs << '\n';
    }
  }
  ofs.close();
}

void waveform_from_file(vector<waveform*> *waveformVP,string file_prefix,params* sim_paramsP)
{
  string filename = file_prefix + "_waveform_jiggled_turn_coord.txt";
  ifstream ifs(filename);
  vector<MatrixXd> jiggled_turn_coord;
  
  int i_neuron = -1;
  
  while (ifs.good()) {
    string val;
    // draw number
    getline(ifs,val,':');
    int i_draw = atoi(val.c_str());
    // check if new set of draws started
    if (i_draw == 0) {
      i_neuron++;
      // if previous set of draws complete, save waveform and clear jiggled_turn_coord to restart
      if (i_neuron > 0) {
        waveform *wP = new waveform(jiggled_turn_coord,sim_paramsP);
        waveformVP->push_back(wP);
        jiggled_turn_coord.clear();
      }
    }
    // n_rows, n_col
    getline(ifs,val,',');
    int n_rows = atoi(val.c_str());
    getline(ifs,val,':');
    int n_cols = atoi(val.c_str());
    // create a draw
    MatrixXd jiggled_turn(n_rows,n_cols);
    for (int i_row=0; i_row < n_rows; i_row++)
    {
      for (int i_col=0; i_col < n_cols; i_col++)
      {
        if (i_row != n_rows-1 || i_col != n_cols-1)
          getline(ifs,val,',');
        else
          getline(ifs,val,'\n');
        jiggled_turn(i_row,i_col) = atof(val.c_str());
      }
    }
    jiggled_turn_coord.push_back(jiggled_turn);
  }
  ifs.close();
}

void neuron_from_file(vector<neuron*> *neuronVP,vector<waveform*> *waveformVP,string file_prefix,std::default_random_engine* generatorP,params* sim_paramsP)
{
  int n_neurons = waveformVP->size();
  
  // get locations
  vector<location*> locV;
  
  string filename_loc = file_prefix +  "_neuron_loc.txt";
  ifstream ifs_loc(filename_loc);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    string val;
    double x,y,z;
    
    getline(ifs_loc,val,',');
    x = atof(val.c_str());
    getline(ifs_loc,val,',');
    y = atof(val.c_str());
    getline(ifs_loc,val,'\n');
    z = atof(val.c_str());
    location *lP = new location(x,y,z);
    locV.push_back(lP);
  }
  //cout << locV[2]->loc(0) << " " << locV[2]->loc(1) << " " << locV[2]->loc(2) << endl;
  ifs_loc.close();
  
  // get directions
  vector<direction*> dirV;
  
  string filename_dir = file_prefix +  "_neuron_dir.txt";
  ifstream ifs_dir(filename_dir);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    string val;
    double x,y,z;
    
    getline(ifs_dir,val,',');
    x = atof(val.c_str());
    getline(ifs_dir,val,',');
    y = atof(val.c_str());
    getline(ifs_dir,val,'\n');
    z = atof(val.c_str());
    direction *dP = new direction(x,y,z);
    dirV.push_back(dP);
  }
  ifs_dir.close();
  
  // get spike times
  vector<vector<int>> spike_timesV;
  
  string filename_spike = file_prefix + "_neuron.txt";
  ifstream ifs_spike(filename_spike);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    string line;
    getline( ifs_spike, line );
    istringstream iss( line );
    string result;
    vector<int> timesV;
    while( std::getline( iss, result, ',' ) )
    {
      int t = atoi(result.c_str());
      if (t!=-10000)
        timesV.push_back(t);
    }
    spike_timesV.push_back(timesV);
  }
  ifs_spike.close();
  
  // get index spikes drawn
  vector<MatrixXi> index_spikes_drawnV;
  
  string filename_spike_index = file_prefix + "_index_spikes_drawn.txt";
  ifstream ifs_spike_index(filename_spike_index);
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    string line;
    getline( ifs_spike_index, line );
    istringstream iss( line );
    string result;
    MatrixXi index_spikes_drawn(spike_timesV[i_neuron].size(),1);
    int i_draw = 0;
    while( std::getline( iss, result, ',' ) )
    {
      int t = atoi(result.c_str());
      if (t!=-10000) {
        index_spikes_drawn(i_draw,0) = t;
      }
      i_draw++;
    }
    index_spikes_drawnV.push_back(index_spikes_drawn);
  }
  
  // create neurons
  for (int i_neuron = 0; i_neuron < n_neurons; i_neuron++)
  {
    neuron *nP = new neuron(spike_timesV[i_neuron],index_spikes_drawnV[i_neuron],locV[i_neuron],dirV[i_neuron],(*waveformVP)[i_neuron],generatorP,sim_paramsP);
    neuronVP->push_back(nP);
  }
  
}

void noisenode_from_file(vector<noisenode*> *noisenodeVP,string file_prefix,params* sim_paramsP)
{
  // get noisenode locations from file
  vector<location*> locV;
  
  string filename_loc = file_prefix + "_noisenode_loc.txt";
  ifstream ifs_loc(filename_loc);
  while (ifs_loc.good())
  {
    string val;
    double x,y,z;
    
    getline(ifs_loc,val,',');
    // check for an empty line
    if (val.length()!=0)
    {
      x = atof(val.c_str());
      getline(ifs_loc,val,',');
      y = atof(val.c_str());
      getline(ifs_loc,val,'\n');
      z = atof(val.c_str());
      location *lP = new location(x,y,z);
      locV.push_back(lP);
    }
  }
  ifs_loc.close();
  int n_noisenodes = locV.size();
  
  // get noisenode traces from file
  vector<VectorXd> traceV;
  
  string filename_trace = file_prefix + "_noisenode.txt";
  ifstream ifs_trace(filename_trace);
  for (int i_nn = 0; i_nn < n_noisenodes; i_nn++)
  {
    VectorXd trace(sim_paramsP->n_pts);
    string line;
    getline( ifs_trace, line );
    istringstream iss( line );
    string result;
    int i_pt = 0;
    while( std::getline( iss, result, ',' ) )
    {
      trace(i_pt) = atof(result.c_str());
      i_pt++;
    }
    traceV.push_back(trace);
  }

  // create noisenodes
  for (int i_nn = 0; i_nn < n_noisenodes; i_nn++) {
    noisenode *nP = new noisenode(locV[i_nn],traceV[i_nn],sim_paramsP);
    noisenodeVP->push_back(nP);
  }
}

void read_detector_loc_file(string detector_loc_file,std::unordered_map<std::string,int> *detector_locP)
{
  ifstream ifs_loc(detector_loc_file);
  while (ifs_loc.good())
  {
    string val;    
    getline(ifs_loc,val);
    // account for possibly having a blank last line
    if (val.length()!=0)
      (*detector_locP)[val.c_str()] = 1;
  }
  
  cout << detector_locP->size() << " location(s) already been calculated" << endl;
  
  ifs_loc.close();
}

void prep_detectors(vector<detector*> *detectorVP,vector<neuron*> *neuronVP,vector<noisenode*> *noisenodeVP,std::unordered_map<std::string,int> *detector_locP,string out_file_prefix,std::default_random_engine* generatorP,params* sim_paramsP)
{
  vector<location*> locationV;
  
  // calc detector locations

  // locations specified case
  if (!sim_paramsP->rand_detectors) {
    int i_row = 0;
    int i_col = 0;
    for (int i_detector=0; i_detector<sim_paramsP->n_detectors; i_detector++) {

      // coord for next location
      double loc_x = -sim_paramsP->dist_betw_electrodes*(sim_paramsP->n_electrode_columns-1)/2.0 + i_col*sim_paramsP->dist_betw_electrodes;
      double loc_y = -sim_paramsP->dist_betw_electrodes*(sim_paramsP->n_electrode_rows-1)/2.0 + i_row*sim_paramsP->dist_betw_electrodes;

      //only add if that location does not exists in the locations loaded from file
      stringstream loc_sstr;
      loc_sstr << loc_x << "," << loc_y << ",0";
      string loc_string = loc_sstr.str();    
      unordered_map<string,int>::const_iterator it_loc_string = detector_locP->find(loc_string);
      if (it_loc_string == detector_locP->end()) {
        location *lP = new location(loc_x,loc_y,0);
        locationV.push_back(lP);
      }

      // go to next row when necessary
      if(i_col == sim_paramsP->n_electrode_columns-1) {
        i_row++;
        i_col=0;
      } else i_col++;
    }
  }
  // locations random case
  else {
    for (int i_detector=0; i_detector<sim_paramsP->n_detectors; i_detector++) {
      // pick rand coord for next location
      double loc_x = rand_uniform(-sim_paramsP->x_edge,sim_paramsP->x_edge,generatorP);
      double loc_y = rand_uniform(-sim_paramsP->y_edge,sim_paramsP->y_edge,generatorP);

      // add to location vector
      location *lP = new location(loc_x,loc_y,0);
      locationV.push_back(lP);
    }
  }

  // reset number of detectors
  sim_paramsP->n_detectors = locationV.size();
  
  // output info
  cout << sim_paramsP->n_detectors << " detectors to calculate" << endl;
  
  // save detector locations to file
  save_detector_locations_from_locV(out_file_prefix,&locationV);

  // create detectors
  for (int i_detector=0; i_detector<sim_paramsP->n_detectors; i_detector++) {
    cout << i_detector << " ";
    cout.flush();
    detector *dP = new detector(locationV[i_detector],neuronVP,noisenodeVP,generatorP,sim_paramsP);
    detectorVP->push_back(dP);
  }
  cout << endl;
}

void prep_noisenodes(vector<noisenode*> *noisenodeVP,std::default_random_engine* generatorP,params* sim_paramsP)
{
  int n_gauss_noise_nodes_x = floor((2*sim_paramsP->volume_radius)/sim_paramsP->dist_betw_gauss_noise_nodes)+1;
  int n_gauss_noise_nodes_y = floor((2*sim_paramsP->volume_radius+sim_paramsP->volume_cylinder_height)/sim_paramsP->dist_betw_gauss_noise_nodes)+1;
  
  for (int i_nn_x=0;i_nn_x<n_gauss_noise_nodes_x;i_nn_x++)
  {
    for (int i_nn_y=0;i_nn_y<n_gauss_noise_nodes_y;i_nn_y++)
    {
      double loc_x = -sim_paramsP->volume_radius/2 + i_nn_x*sim_paramsP->dist_betw_gauss_noise_nodes;
      double loc_y = -(sim_paramsP->volume_radius/2+sim_paramsP->volume_cylinder_height) + i_nn_y*sim_paramsP->dist_betw_gauss_noise_nodes;
      double loc_z = 0; 
      location *lP = new location(loc_x,loc_y,loc_z);
      noisenode *nP = new noisenode(lP,generatorP,sim_paramsP);
      noisenodeVP->push_back(nP);
    }
  }
  
}


void prep_neurons(vector<neuron*> *neuronVP,waveform template_w,std::default_random_engine* generatorP,params* sim_paramsP)
{
  // distribute neurons
  
  // first distribute in a cube
  double volume_box = pow(2 * sim_paramsP->volume_radius + sim_paramsP->volume_cylinder_height,3);
  int n_neurons_box = ceil(volume_box*sim_paramsP->density_neurons);
  
  MatrixXd poss_locs = rand_matrix_uniform(3,n_neurons_box,-(sim_paramsP->volume_radius + sim_paramsP->volume_cylinder_height/2),sim_paramsP->volume_radius + sim_paramsP->volume_cylinder_height/2,generatorP);

  int n_inactive = 0;
  //cout << poss_locs.size() << endl;
  
  if (sim_paramsP->log_normal_firing_distrib) {
    cout << "log normal distribution of firing rates" << endl;
    if (sim_paramsP->poisson_rate_lognorm_mean != 0)
    {
      cout << "mean " << sim_paramsP->poisson_rate_lognorm_mean << " var " << sim_paramsP->poisson_rate_lognorm_var << endl;
      sim_paramsP->poisson_rate_lognorm_mu = -(0.5)*log((sim_paramsP->poisson_rate_lognorm_var / pow(sim_paramsP->poisson_rate_lognorm_mean,2) + 1)/pow(sim_paramsP->poisson_rate_lognorm_mean,2));
      sim_paramsP->poisson_rate_lognorm_sigma = pow(2*log(sim_paramsP->poisson_rate_lognorm_mean)-2*sim_paramsP->poisson_rate_lognorm_mu,0.5);
    }
    cout << "mu " << sim_paramsP->poisson_rate_lognorm_mu << " sigma " << sim_paramsP->poisson_rate_lognorm_sigma << endl;
  } else {
    cout << "uniform distribution of firing rates" << endl;
  }

  // create all once for scoping reasons
  std::uniform_real_distribution<double> d_uniform_real(sim_paramsP->poisson_rate_min,sim_paramsP->poisson_rate_max);
  std::lognormal_distribution<double> d_lognormal(sim_paramsP->poisson_rate_lognorm_mu,sim_paramsP->poisson_rate_lognorm_sigma);
  std::bernoulli_distribution d_bernoulli(sim_paramsP->active_neuron_rate);

  for (int i_poss_neuron=0; i_poss_neuron<n_neurons_box; i_poss_neuron++) {
    if (pow(pow(poss_locs(0,i_poss_neuron),2) + pow(poss_locs(2,i_poss_neuron),2),0.5) <= sim_paramsP->volume_radius && (abs(poss_locs(1,i_poss_neuron)) < 0.5*sim_paramsP->volume_cylinder_height || pow(pow(poss_locs(0,i_poss_neuron),2) + pow(abs(poss_locs(1,i_poss_neuron))-0.5*sim_paramsP->volume_cylinder_height,2) + pow(poss_locs(2,i_poss_neuron),2),0.5) < sim_paramsP->volume_radius))
    {
      // get random orientation
      direction *dP = new direction(rand_matrix_uniform(3,1,-1,1,generatorP));
      location *lP = new location(poss_locs.col(i_poss_neuron));
      // get rate
      // init as 0
      double rate = 0;
      // determine if active neuron
      bool active_neuron = d_bernoulli(*generatorP);
      // if active, generate rate
      if (active_neuron) {
        if (sim_paramsP->log_normal_firing_distrib) {
          rate = d_lognormal(*generatorP);
        } else {
          rate = d_uniform_real(*generatorP);
        }
      } else {
        n_inactive++;
      }
      // create new wave
      waveform *wP = new waveform(template_w.jiggle_turns(),sim_paramsP->ratio_stdev_jiggle_x_draw,sim_paramsP->ratio_stdev_jiggle_y_draw,sim_paramsP->n_draws_per_neuron,generatorP,sim_paramsP);
      // create a new neuron
      neuron *nP = new neuron(rate,lP,dP,wP,generatorP,sim_paramsP);
      // bucket of neurons
      neuronVP->push_back(nP);
    }
  }
  cout << "n inactive: " << n_inactive << endl;
}

/**

void writeToCSVfile(string name, MatrixXd matrix)
{
  ofstream file(name.c_str());
  
  for(int  i = 0; i < matrix.rows(); i++){
    for(int j = 0; j < matrix.cols(); j++){
      string str = lexical_cast<std::string>(matrix(i,j));
      if(j+1 == matrix.cols()){
        file<<str;
      }else{
        file<<str<<',';
      }
    }
    file<<'\n';
  }
**/
/**

todo:
- move generator under params
- move all defaults into params files
 - do a good indices start at zero check / etc / search for 1

**/
