/* Copyright 2017 Antonia Reiter */
/* no obligations - feel free to copy/reuse/modify as you like*/

#ifndef SRC_FUSIONEKF_H_
#define SRC_FUSIONEKF_H_

#include <vector>
#include <string>
#include <fstream>
#include "measurement_package.h"
#include "Eigen/Dense"
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
 public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;

 private:
  // check whether the tracking toolbox was initialized
  // or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  double previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::MatrixXd R_laser_;  // R Matrix for Laser
  Eigen::MatrixXd R_radar_;  // R Matrix for Radar
  Eigen::MatrixXd H_laser_;  // H Matrix for Laser
  Eigen::MatrixXd Hj_;       // H Matrix for Radar (--> will be Jacobian Matrix)

  double noise_ax;
  double noise_ay;
};

#endif   // SRC_FUSIONEKF_H_
