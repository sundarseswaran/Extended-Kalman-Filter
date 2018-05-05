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

  // default values
  is_initialized_ = false;
  previous_timestamp_ = 0;
  ekf_.PI_ = atan(1)*4;

  // default noise component
  noise_ax = 9;
  noise_ay = 9;


  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  MatrixXd P_ = MatrixXd(4, 4);
  MatrixXd F_ = MatrixXd(4, 4);
  MatrixXd Q_ = MatrixXd(4, 4);

  VectorXd x_ = VectorXd(4);

  Hj_ = MatrixXd(3, 4);

  // covariance matrix for laser measurement
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  // covariance matrix for radar measurement
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //state covariance matrix P
  P_ <<  1, 0, 0,    0,
         0, 1, 0,    0,
         0, 0, 1000, 0,
         0, 0, 0,    1000;

  //the initial transition matrix F_
  F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

  // initialize EKF
  ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
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

    // step one.
    float px = 0;
    float py = 0;
    float vx = 0;
    float vy = 0;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      /**
        Convert radar from polar to cartesian coordinates and initialize state.
      */

      float rho = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float rhodot = measurement_pack.raw_measurements_[2];

      px = rho * cos(theta);
      py = rho * sin(theta);

      vx = rhodot * cos(theta);
      vy = rhodot * sin(theta);

//      if(fabs(px) < 0.0001) {
//        px = 1;
//        ekf_.P_(0,0) = 1000;
//      }
//
//      if(fabs(py) < 0.0001) {
//        py = 1;
//        ekf_.P_(1,1) = 1000;
//      }

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    ekf_.x_ << px, py, vx, vy;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // time delta - in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt2 = dt * dt;
  float dt3 = pow(dt,3);
  float dt4 = pow(dt,4);

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  ekf_.Q_ = Eigen::MatrixXd(4,4);

  ekf_.Q_ << dt4 * noise_ax / 4, 0, dt3 * noise_ax / 2, 0,
          0, dt4 * noise_ay / 4, 0, dt3 * noise_ay /2,
          dt3 * noise_ax / 2, 0, dt2 * noise_ax, 0,
          0, dt3 * noise_ay / 2, 0, dt2 * noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  Eigen::VectorXd rawMsmt = measurement_pack.raw_measurements_;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    if (rawMsmt(1) > ekf_.PI_) {
      rawMsmt(1) = rawMsmt(1)-2*ekf_.PI_;
    }

    if (rawMsmt(1) < -ekf_.PI_) {
      rawMsmt(1) = rawMsmt(1)+2*ekf_.PI_;
    }

    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;

//      cout << "RADAR" << endl;
//      cout << rawMsmt;

    ekf_.UpdateEKF(rawMsmt);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

//      cout << "LASER" << endl;
//      cout << rawMsmt;

    ekf_.Update(rawMsmt);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
