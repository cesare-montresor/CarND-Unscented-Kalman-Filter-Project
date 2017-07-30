#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size() || estimations.size() == 0){
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }
  
  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  
  //calculate the mean
  rmse = rmse/estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}

VectorXd Tools::Polar2Cartesian(float distance,float angle,float velocity){
  VectorXd result=VectorXd(4);
  
  float px = distance * std::cos(angle);
  float py = distance * std::sin(angle);
  float vx = velocity * std::cos(angle);
  float vy = velocity * std::sin(angle);
  
  //cout<<"px: \t"<<(int)px<<"\tpy: \t"<<(int)py<<"distance: \t"<<distance<<"\tangle: \t"<<angle<<"\tvelocity: \t"<<velocity<<"\n";
  result << px, py, vx, vy;
  return result;
}

float Tools::NormalizeAngle(float angle){
  while ( angle > M_PI ){
    angle -= 2.*M_PI;
  }
  
  while ( angle < -M_PI ){
    angle += 2.*M_PI;
  }
  
  return angle;
}

