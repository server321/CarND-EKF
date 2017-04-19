#include "kalman_filter.h"
#define PI 3.141592653589793238463

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;  // object state
  P_ = P_in;  // object covariance matrix
  F_ = F_in;  // state transition matrix
  H_ = H_in;  // measurement matrix
  R_ = R_in;  // measurement covariance matrix
  Q_ = Q_in;  // process covariance matrix
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_ ;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_; // error calculation
  KF(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd z_pred = hx_;
  VectorXd y = z - z_pred;

  // if phi is not in the range (-pi, pi), put it in that range

  bool in_range = false;
  while (in_range == false) {
    if (y(1) > PI) {
        y(1) = y(1) - PI*2;
        }
    else if (y(1) < -PI) {
        y(1) = y(1) + PI*2;
        }
	else {
           in_range = true;
    }
  }

  // Calculations are essentially the same to the Update function
  KF(y);
}

// Universal update Kalman Filter step.
void KalmanFilter::KF(const VectorXd &y){
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  // New state
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
