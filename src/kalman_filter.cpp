#include "kalman_filter.h"
#include <math.h>
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

void KalmanFilter::UpdateKF(const VectorXd &y) {

	MatrixXd Ht = H_.transpose();
	MatrixXd P_Ht = P_ * Ht;

	MatrixXd S = H_ * P_Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_Ht * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::Predict() {
	MatrixXd Ft = F_.transpose();
	x_ = F_ * x_;
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	UpdateKF(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
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

	// Its ok if py is 0.0 woot woot
	// atan2 normalizes the angle from pi to -pi (http://www.cplusplus.com/reference/cmath/atan2/)
	phi = atan2(py, px);

	if (fabs(rho) > 0.0001)
	{
		rho_dot = (px * vx + py * vy) / rho;
	}
	else
	{
		// Saw a couple of students substituting with a very small value. Giving it a try
		rho_dot = (px * vx + py * vy) / 0.0001;
	}

	VectorXd hX(3);
	hX << rho, phi, rho_dot;
	VectorXd y = z - hX;

	// Additional normalization
	while (y(1)>M_PI)
	{
		y(1) -= 2 * M_PI;
	}
	while (y(1)<-M_PI)
	{
		y(1) += 2 * M_PI;
	}
	//std::cout << "EKF Update: " << y << std::endl;
	UpdateKF(y);
}


