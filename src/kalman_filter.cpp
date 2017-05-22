#include "kalman_filter.h"

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

void KalmanFilter::UpdateKF(const VectorXd &y) {

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

void KalmanFilter::Predict() {
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	UpdateKF(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  * TODO:
  * update the state by using Extended Kalman Filter equations
  **/
	// y = z - h(x`)
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	float rho = sqrt((px * px) + (py * py));
	float phi = 0.0;
	float rho_dot = 0.0;

	if (fabs(px) > 0.001)
	{
		phi = atan2(py, px);
	}

	if (fabs(rho) > 0.001)
	{
		rho_dot = (px * vx + py * vy) / rho;
	}
	
	VectorXd hX(3);
	hX << rho, phi , rho_dot;

	VectorXd y = z - hX;

	UpdateKF(y);
}


