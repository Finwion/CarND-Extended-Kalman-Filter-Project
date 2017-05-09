#include "kalman_filter.h"
#include <iostream>


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
  CompleteMeasurementUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  double px = x_(0); 
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
	VectorXd Hx(3);

  double rho = sqrt(pow(px, 2) + pow(py, 2));
  if (fabs(rho) < 0.0001)
  {
    std::cout << "UpdateEKF () - Error - Division by Zero" << std::endl;
    return;
  } 

  double phi = atan2(py, px);

  double rho_dot = (px*vx + py*vy)/(rho);
  Hx << rho, phi, rho_dot;
	VectorXd y = z - Hx;

  //Normalize y so that it is between -PI and PI
  if(y[1] < -M_PI){
    y[1] += 2*M_PI;
  }

  if(y[1] > M_PI) {
    y[1] -= 2*M_PI;
  } 
  CompleteMeasurementUpdate(y);
}


void KalmanFilter::CompleteMeasurementUpdate(const VectorXd &y) {
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}