/* Copyright 2017 Antonia Reiter */
/* no obligations - feel free to copy/reuse/modify as you like*/
#ifndef SRC_UKF_H_
#define SRC_UKF_H_

#include <vector>
#include <string>
#include <fstream>
#include "measurement_package.h"
#include "Eigen/Dense"
#include "FusionEKF.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
 public:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  MatrixXd Xsig_aug_;

  ///* time when the state is true, in us
  double previous_timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  ///* Tools for helper functions
  Tools tools;

  ///* The "neutral" state vector for EKF and UKF: px, px, vx, vy
  VectorXd averaged_;

  ///* Time steps since the first measurement (actually an integer)
  double Time_Step_;

  ///* Threshold for the 95% mark of radar NIS
  double NIS_radar_threshold_;

  ///* Threshold for the 95% mark of laser NIS
  double NIS_laser_threshold_;

  ///* Number of occurences where the radar NIS is greater than 95% mark
  double NIS_radar_over_threshold_;

  ///* Number of occurences where the lidar NIS is greater than 95% mark
  double NIS_laser_over_threshold_;

  ///* Instance of the Extended Kalman Filter --> I want to have both in one SW
  FusionEKF ekf_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
   * Calculation of the agumented sigma points
   */
  void AugmentedSigmaPoints();

  /**
   * Prediction of the Sigma points
   * @param delta_t which is the time delta since the last prediction (in seconds)
   */  
  void SigmaPointPrediction(double delta_t);

  /**
   * Prediction of the Mean and Covariance Matrix
   */ 
  void PredictMeanAndCovariance();
};

#endif  // SRC_UKF_H_
