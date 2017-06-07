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
  VectorXd rmse = VectorXd::Zero(4);
  VectorXd residual = VectorXd::Zero(4);
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
  for(int i=0; i < estimations.size(); ++i){
      residual = estimations[i] - ground_truth[i];
      rmse += residual.cwiseProduct(residual);
  }

  //calculate the mean
  rmse *= 1.0 / estimations.size();
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
  MatrixXd Hj = MatrixXd::Zero(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //check division by zero
  if (fabs(c1) < 0.0000000001) {
      cout << "CalculateJacobian() - Error - Division by Zero" << endl;
      return Hj;
  }
  //compute the Jacobian matrix
  Hj(0,0) = Hj(2,2) = px / c2;
  Hj(0,1) = Hj(2,3) = py / c2;
  Hj(1,0) = - py / c1;
  Hj(1,1) = px / c1;
  Hj(2,0) = py * (vx * py - vy * px) / c3;
  Hj(2,1) = px * (vy * px - vx * py) / c3;
  return Hj;
}
