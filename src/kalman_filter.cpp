#include "kalman_filter.h"

#define EPS 0.0000001 // A small number

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = (H_ * P_ * Ht) + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred(z.rows());
  // temporal variables for readability
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  z_pred(0) = sqrt(px * px + py * py);
  // if rho = 0 (x=0, and y=0) we are at the origin. Use small rho instead
  z_pred(0) = std::max(z_pred(0), EPS);
  // same with px = 0
  if (fabs(px) < EPS) {
    px = EPS;
  }
  z_pred(1) = atan2(py, px);
  z_pred(2) = (px * vx + py * vy) / (z_pred(0));

  VectorXd y = z - z_pred;
  // correct y(1) to be between -pi and pi. We can get pi from 4*atan(1)
  y(1) = std::fmod(y(1), (4 * tan(1)));
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
