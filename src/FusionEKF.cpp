/* Copyright 2017 Antonia Reiter */
/* no obligations - feel free to copy/reuse/modify as you like*/

#include <iostream>
#include <string>
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"

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

  // measurement matrix for laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // measurement matrix for radar
  // --> has to be done dynamically as the calculation of the Jacobian
  //     depends on the state

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd::Zero(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

  // set the acceleration noise components
  //   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  noise_ax = 9.;
  noise_ay = 9.;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  // cout << "" << endl;
  // cout << "" << endl;
  // cout << "-------------------------------------------------------" << endl;
  // cout << " Processing new measurement pack for sensor:   " << measurement_pack.sensor_type_ << endl;  
  // cout << " Raw measurement: \n" << measurement_pack.raw_measurements_ << endl;

  string logstr = "";

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_ &&
      measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        ekf_.x_ = VectorXd(4);

    ekf_.x_ << 0, 0.5, 0, 0;
    return;
  }

  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
    */
    // first measurement
    ekf_.x_ = VectorXd(4);

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      // cout << "Polar coordinates: \n" <<
      // measurement_pack.raw_measurements_ << endl;
      logstr += tools.printVector("Init:Radar:",
                                  measurement_pack.raw_measurements_);
      ekf_.x_ = tools.Polar2Cartesian(measurement_pack.raw_measurements_);
      // cout << "\nis being converted to Cartesian coordinates: \n" <<
      // ekf_.x_ << endl;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      logstr += tools.printVector("Init:LASER:",
                                  measurement_pack.raw_measurements_);
      ekf_.x_ << measurement_pack.raw_measurements_[0],
      measurement_pack.raw_measurements_[1],
      0,
      0;
    }




    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //logstr += tools.printVector("Polar Radar:", measurement_pack.raw_measurements_);
    //logstr += tools.printVector("Cartesian Radar:", tools.Polar2Cartesian(measurement_pack.raw_measurements_));
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    //logstr += tools.printVector("Cartesian LASER:", measurement_pack.raw_measurements_);
  }

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;


  /*****************************************************************************
   *  1. Step: Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
     * Time is measured in seconds.
     * Update the process noise covariance matrix.
   */

  // 1. Modify the F matrix so that the time is integrated
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd::Zero(4, 4);  // don't forget to init the Matrix with 0 !
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd::Zero(4, 4);  // don't forget to init the Matrix with 0 !

  // calculate the Q-Matrix is described
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
            0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
            dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
            0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  // finally execute the predict step of the Kalman filter
  ekf_.Predict();

  // cout << "------ Predict finished ------ " << endl;

  /*****************************************************************************
   *  2. Step: Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = this->R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);  // for Radar use the Jacobian!!
    logstr += ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = this->R_laser_;
    ekf_.H_ = this->H_laser_;
    logstr += ekf_.Update(measurement_pack.raw_measurements_);
  }

  // cout << "------ Update finished ------ " << endl;

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
  cout << logstr << endl;
}
