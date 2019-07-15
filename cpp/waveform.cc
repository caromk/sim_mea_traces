#include "waveform.h"
using namespace std;

waveform::waveform(MatrixXd template_turn_coord,double ratio_stdev_jiggle_x,double ratio_stdev_jiggle_y,int n_draw,std::default_random_engine* generatorP,params* sim_paramsP) {
  
  // initialize variables
  this->template_turn_coord = template_turn_coord;
  this->generatorP = generatorP;
  this->sim_paramsP = sim_paramsP;
  this->ratio_stdev_jiggle_x = ratio_stdev_jiggle_x;
  this->ratio_stdev_jiggle_y = ratio_stdev_jiggle_y;
  this->n_draw = n_draw;
  
  this->n_turns = template_turn_coord.row(0).size();
  
  // turn ratio_stdev_jiggle into stdev_jiggle by looking at the waves
  
  // for x, use minimum x step across beginning and end waveform
  double min_x_step = min((this->template_turn_coord.row(0).segment(4,n_turns-5) - this->template_turn_coord.row(0).segment(3,n_turns-5)).minCoeff(),(this->template_turn_coord.row(2).segment(4,n_turns-5) - this->template_turn_coord.row(2).segment(3,n_turns-5)).minCoeff());
  // for y, use max peak across beginning and end waveform
  // todo/note: this assumes positive side is larger
  double max_y_peak = max(this->template_turn_coord.row(1).maxCoeff(),this->template_turn_coord.row(3).maxCoeff());
  
  stdev_jiggle_x = this->ratio_stdev_jiggle_x*min_x_step;
  stdev_jiggle_y = this->ratio_stdev_jiggle_y*max_y_peak;
  
  // jiggle
  
  create_n_draws();
  
}

waveform::waveform(vector<MatrixXd> jiggled_turn_coord,params* sim_paramsP)
{
  this->jiggled_turn_coord = jiggled_turn_coord;
  this->sim_paramsP = sim_paramsP;
  this->n_draw = jiggled_turn_coord.size();
  this->n_turns = jiggled_turn_coord[0].cols();
  
  set_x_values();
}

void waveform::set_time_values() {
  set_x_values();

  //todo: this is totally hacky to put here, going out of order in the constructor could screw everything up, but...
  // readjust time in turn coords to be the correct length
  double max_time_spike_coord = max(template_turn_coord.row(0).maxCoeff(),template_turn_coord.row(2).maxCoeff());
  double spk_time_conversion = sim_paramsP->spk_length/max_time_spike_coord;
  template_turn_coord.row(0) *= spk_time_conversion;
  template_turn_coord.row(2) *= spk_time_conversion;
  stdev_jiggle_x *= spk_time_conversion;
}

void waveform::set_x_values() {
  x_values = VectorXd::Zero(sim_paramsP->n_pts_per_spike);
  for (int i_pt = 1; i_pt < sim_paramsP->n_pts_per_spike; i_pt++) {
    x_values(i_pt) = x_values(i_pt-1) + 1.0/sim_paramsP->sampling_rate;
  }
}

void waveform::create_n_draws() {
  // set up/fix time
  if(n_draw > 0) set_time_values();
  
  // jiggle
  for (int i=0; i<n_draw; i++) {
    jiggled_turn_coord.push_back(jiggle_turns());
  }
}

MatrixXd waveform::jiggle_turns() {
  
  MatrixXd new_template_turn_coord(template_turn_coord.rows(),template_turn_coord.cols());
  new_template_turn_coord = template_turn_coord;
  
  double stdev_cutoff = 2;
  
  // jiggle x coords
  new_template_turn_coord.row(0).segment(3,n_turns-5) += rand_matrix_norm_amp_limit(1,n_turns-5,0,stdev_jiggle_x,stdev_cutoff,generatorP);
  new_template_turn_coord.row(2).segment(3,n_turns-5) += rand_matrix_norm_amp_limit(1,n_turns-5,0,stdev_jiggle_x,stdev_cutoff,generatorP);
  // jiggle y coords
  new_template_turn_coord.row(1).segment(3,n_turns-5) += rand_matrix_norm_amp_limit(1,n_turns-5,0,stdev_jiggle_y,stdev_cutoff,generatorP);
  new_template_turn_coord.row(3).segment(3,n_turns-5) += rand_matrix_norm_amp_limit(1,n_turns-5,0,stdev_jiggle_y,stdev_cutoff,generatorP);
  
  return new_template_turn_coord;
}

// todo: pass a pointer #memory
MatrixXd waveform::waves_from_distance(MatrixXi *index_spikes_drawnP,double dist_from_normal) {
  
  int n_spks_selected = index_spikes_drawnP->size();
  
  MatrixXd waves_drawn(n_draw,sim_paramsP->n_pts_per_spike);
  MatrixXd waves_selected(n_spks_selected,sim_paramsP->n_pts_per_spike);
  
  for(int i_draw=0;i_draw<n_draw;i_draw++) {
    VectorXd curr_turns_x(n_turns);
    VectorXd curr_turns_y(n_turns);
    if (dist_from_normal==0) {
      curr_turns_x = jiggled_turn_coord[i_draw].row(0);
      curr_turns_y = jiggled_turn_coord[i_draw].row(1);
    }
    else if (dist_from_normal <= sim_paramsP->dist_to_fin) {
      curr_turns_x = jiggled_turn_coord[i_draw].row(0) + ((jiggled_turn_coord[i_draw].row(2)-jiggled_turn_coord[i_draw].row(0)) / sim_paramsP->dist_to_fin) * dist_from_normal;
      curr_turns_y = jiggled_turn_coord[i_draw].row(1) +  ((jiggled_turn_coord[i_draw].row(3)-jiggled_turn_coord[i_draw].row(1)) / sim_paramsP->dist_to_fin) * dist_from_normal;
    }
    else if (dist_from_normal <= sim_paramsP->dist_to_fin + sim_paramsP->dist_to_zero)
    {
      curr_turns_x = jiggled_turn_coord[i_draw].row(0);
      curr_turns_y = jiggled_turn_coord[i_draw].row(3) - (jiggled_turn_coord[i_draw].row(3) / (sim_paramsP->dist_to_zero - sim_paramsP->dist_to_fin)) * (dist_from_normal - sim_paramsP->dist_to_fin);
    }
    else {
      curr_turns_x = VectorXd::Zero(n_turns);
      curr_turns_y = VectorXd::Zero(n_turns);
    }

    //cout << endl << "curr_turns_x = " << curr_turns_x << endl;
    //cout << endl << "curr_turns_y = " << curr_turns_y << endl;

    waves_drawn.row(i_draw) = calc_wave(curr_turns_x,curr_turns_y);

    //cout << endl << "waves_draw current row = " << waves_drawn.row(i_draw) << endl;
  }
  
  for(int i_selected=0;i_selected < n_spks_selected;i_selected++)
  {
    waves_selected.row(i_selected) = waves_drawn.row((*index_spikes_drawnP)(i_selected));
  }
  
  return waves_selected;
}

VectorXd waveform::calc_wave(VectorXd curr_turns_x,VectorXd curr_turns_y) {
  VectorXd y_values(sim_paramsP->n_pts_per_spike);
  y_values=VectorXd::Zero(sim_paramsP->n_pts_per_spike);

  //cout << "n_pts_per_spike = " << sim_paramsP->n_pts_per_spike << endl;
  //cout << "n_turns = " << n_turns << endl;

  VectorXd x_values_to3 = x_values.array().pow(3);
  VectorXd x_values_to2 = x_values.array().pow(2);
  VectorXd x_values_to0(sim_paramsP->n_pts_per_spike); x_values_to0.fill(1);
  
  MatrixXd cubic_coeff = fit_cubics_n_turns(curr_turns_x,curr_turns_y);

  //cout << "cubic_coeff " << cubic_coeff << endl;

  int i_curve_start = 0;
  for(int i_turn=0;i_turn<n_turns-1;i_turn++) {
    int i_curve_end = ceil(sim_paramsP->sampling_rate*curr_turns_x(i_turn+1))-1;
    y_values.segment(i_curve_start,i_curve_end-i_curve_start+1) = cubic_coeff(i_turn,0)*x_values_to3.segment(i_curve_start,i_curve_end-i_curve_start+1) + cubic_coeff(i_turn,1)*x_values_to2.segment(i_curve_start,i_curve_end-i_curve_start+1) + cubic_coeff(i_turn,2) * x_values.segment(i_curve_start,i_curve_end-i_curve_start+1) + cubic_coeff(i_turn,3)*x_values_to0.segment(i_curve_start,i_curve_end-i_curve_start+1);
    i_curve_start = i_curve_end + 1;
  } 
  
  return y_values;

}

VectorXd waveform::fit_cubic_two_turns(double x1,double y1,double x2,double y2) {
  
  VectorXd coeff(4); coeff << 0.0,0.0,0.0,0.0;
  double k, m;
  
  if(x1 != 0 && x2 != 0) {
    k = (y1 - y2) / ((1.0/3.0)*(pow(x1,3)-pow(x2,3)) - (1.0/2.0)*(x1+x2)*(pow(x1,2)-pow(x2,2)) + (x1*x2)*(x1-x2));
    m = y1 - k*((1.0/3.0)*pow(x1,3) - (1.0/2.0)*(x1+x2)*pow(x1,2) + (x1*x2)*x1);
    
    coeff[0] = (1.0/3.0)*k;
    coeff[1] = -(1.0/2.0)*k*(x1+x2);
    coeff[2] = k*(x1*x2);
    coeff[3] = m;

    // check - matches matlab results
    /**
    std::cout << coeff[0] << std::endl;
    std::cout << coeff[1] << std::endl;
    std::cout << coeff[2] << std::endl;
    std::cout << coeff[3] << std::endl;
     **/
  }
  
  return coeff;
  
}


// MatrixXd cubic_coeff = fit_cubics_n_turns(this->template_turn_coord.row(0),this->template_turn_coord.row(1));

MatrixXd waveform::fit_cubics_n_turns(VectorXd turn_coord_x,VectorXd turn_coord_y)
{
  MatrixXd cubic_coeff(n_turns-1,4);
  
  for (int i_turn=0; i_turn < n_turns-1; i_turn++) {
    cubic_coeff.row(i_turn) = fit_cubic_two_turns(turn_coord_x(i_turn),turn_coord_y(i_turn),turn_coord_x(i_turn+1),turn_coord_y(i_turn+1));
  }

  /*
  // check - matches Matlab result
  
   for (int i = 0; i < n_turns-2; i++) {
   for (int j = 0; j < 4; j++) {
   cout << cubic_coeff(i,j) << " ";
   }
   cout << endl;
   }
  */
  
  return cubic_coeff;
    
}

vector<MatrixXd>* waveform::return_jiggled_turn_coord()
{
  return &jiggled_turn_coord;
}