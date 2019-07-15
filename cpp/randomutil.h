#ifndef _randomutil_
#define _randomutil_

#include <random>
#include <iostream>
#include "Eigen/Dense"
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;

inline MatrixXd rand_matrix_norm_amp_limit(int n_row,int n_col,double mean,double stdev,double stdev_cutoff,std::default_random_engine* generator)
{
  
  std::normal_distribution<double> distribution(mean,stdev);
  
  int n_num = n_row*n_col;
  MatrixXd num(n_row,n_col);
  
  double pos_cutoff = stdev*stdev_cutoff+mean;
  double neg_cutoff = -stdev*stdev_cutoff+mean;
//  std::cout << pos_cutoff << endl;
 // std::cout << neg_cutoff << endl;
  
  int i_num = 0;
  int i_row = 0;
  int i_col = 0;
  while (i_num < n_num) {
    double number = distribution(*generator);
//    std::cout << number << endl;
    
    if ((number>neg_cutoff)&&(number<pos_cutoff))
    {
//      std::cout << "row " << i_row << " col " << i_col << endl;
      num(i_row,i_col) = number;
      if (i_row < n_row-1) { i_row++; }
      else {
        i_col++;
        i_row = 0;
      }
      i_num++;
    }
  }
  
  return num;
}

inline MatrixXi rand_matrix_uniform_int(int n_row,int n_col,int min_val,int max_val,std::default_random_engine* generator) {
  
  std::uniform_real_distribution<double> distribution(min_val,max_val+1);
  
  MatrixXi num(n_row,n_col);
  
  for (int i_row=0; i_row<n_row; ++i_row) {
    for (int i_col=0; i_col<n_col; ++i_col) {
      double number = distribution(*generator);
      num(i_row,i_col) = floor(number);
    }
  }
  return num;
}

inline MatrixXd rand_matrix_uniform(int n_row,int n_col,double min_val,double max_val,std::default_random_engine* generator) {
  
  std::uniform_real_distribution<double> distribution(min_val,max_val);
  
  MatrixXd num(n_row,n_col);
  
  for (int i_row=0; i_row<n_row; ++i_row) {
    for (int i_col=0; i_col<n_col; ++i_col) {
      double number = distribution(*generator);
      num(i_row,i_col) = number;
    }
  }
  return num;
}

inline double rand_uniform(double min_val,double max_val,std::default_random_engine* generator) {
  
  std::uniform_real_distribution<double> distribution(min_val,max_val);
  double number = distribution(*generator);
  return number;
  
}

#endif
