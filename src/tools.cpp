#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if (estimations.size() == 0) {
      return rmse;
  }
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here
  if (estimations.size() != ground_truth.size())  {
      return rmse;
    }
  // TODO: The creation of a temporal matrix doubles memory requirements
  // for estimates, bad longterm. They are already in estimations: rewrite to use it.
  MatrixXd residuals(4, estimations.size());
  for(int i=0; i < estimations.size(); ++i){
      VectorXd residual = estimations[i] - ground_truth[i];
      residuals.col(i) = residual.cwiseProduct(residual);
  }

  //calculate the mean
  rmse = residuals.rowwise().mean();
  //calculate the squared root
  rmse = rmse.cwiseSqrt();
  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float p_numerator = px * px + py * py;

  //check division by zero
  if (fabs(p_numerator) < 0.0000000001) {
      cout << "CalculateJacobian() - Error - Division by Zero" << endl;
      Hj.setZero();
      return Hj;
  }
  //compute the Jacobian matrix
  Hj(0,0) = Hj(2,2) = px / sqrt(p_numerator);
  Hj(0,1) = Hj(2,3) = py / sqrt(p_numerator);
  Hj(1,0) = -py / p_numerator;
  Hj(1,1) = px / p_numerator;
  Hj(2,0) = py * (vx * py - vy * px) / pow(p_numerator, 3/2);
  Hj(2,0) = px * (vy * px - vx * py) / pow(p_numerator, 3/2);
  return Hj;
}
