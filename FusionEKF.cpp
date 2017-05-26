#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
	         0, 1, 0, 0;

  Hj_ << 0, 0, 0, 0,
	     0, 0, 0, 0,
	     0, 0, 0, 0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;

  //Tools to get Jacobian
  Tools tools;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    std::cout << "EKF: " << endl;

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 950, 0,
		0, 0, 0, 950;

	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
		0, 1, 0, 1,
		0, 0, 1, 0,
		0, 0, 0, 1;

	ekf_.Q_ = MatrixXd(4, 4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		float rho = measurement_pack.raw_measurements_(0);
		float phi = measurement_pack.raw_measurements_(1);
		float rho_dot = measurement_pack.raw_measurements_(2);

		float px = rho * cos(phi);
		float py = rho * sin(phi);
		float vx = rho_dot * cos(phi);
		float vy = rho_dot * sin(phi);

		ekf_.x_ << px, py, vx, vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = pow(dt, 2);
  float dt_3 = pow(dt, 3);
  float dt_4 = pow(dt, 4);

  // Update the F matrix:
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // set the covariance matrix
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
	  0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
	  dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
	  0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  float rho = measurement_pack.raw_measurements_[0];
	  float phi = measurement_pack.raw_measurements_[1];
	  float rho_dot = measurement_pack.raw_measurements_[2];
	  VectorXd z(3);
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.R_ = R_radar_;
	  z << rho, phi, rho_dot;
	  ekf_.UpdateEKF(z);
  } else {
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  VectorXd z(2);
	  z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
	  ekf_.Update(z);
  }

  // print the output
  //std::cout << "x_ = " << ekf_.x_ << endl;
  //std::cout << "P_ = " << ekf_.P_ << endl;
}
