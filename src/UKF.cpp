#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;
  Time_Step_ = 0;
  NIS_radar_threshold_ = 7;
  NIS_laser_threshold_ = 7;
  NIS_radar_over_threshold_ = 0;
  NIS_laser_over_threshold_ = 0;


  ekf_ = FusionEKF();
  averaged_ = VectorXd(4);


  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  ///* --> State dimension is 5
  n_x_ = 5;
  x_ = VectorXd(n_x_);
  x_.fill(0.0);


  ///* Augmented state dimension because of 2 process noise parameters
  //   std_a_ and std_yawdd_
  n_aug_ = n_x_ + 2;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_;

    // set vector for weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i < 2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  ///* state covariance matrix
  P_ = MatrixXd(n_x_ , n_x_);
  P_.fill(0.0);

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  Xsig_aug_  = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.fill(0.0);

  NIS_radar_ = -1;
  NIS_laser_ = -1;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /***************************************************************************
   *  Initialization
   **************************************************************************/
  if (is_initialized_ == false) {
    // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    x_ << 0, 0.5, 0, 0, 0;
    averaged_ << 0, 0.5, 0, 0;

    // x_ << 1, 1, 1, 1, 0.1;
    // init covariance matrix
    P_ << 0.15,    0, 0, 0, 0,
             0, 0.15, 0, 0, 0,
             0,    0, 1, 0, 0,
             0,    0, 0, 1, 0,
             0,    0, 0, 0, 1;

     if (measurement_pack.sensor_type_ == MeasurementPackage::LASER
         && use_laser_) {
        x_(0) = measurement_pack.raw_measurements_(0);
        x_(1) = measurement_pack.raw_measurements_(1);
      } else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR
         && use_radar_) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        float ro = measurement_pack.raw_measurements_(0);
        float phi = measurement_pack.raw_measurements_(1);
        float ro_dot = measurement_pack.raw_measurements_(2);
        x_(0) = ro     * cos(phi);
        x_(1) = ro     * sin(phi);
      }

    previous_timestamp_ = measurement_pack.timestamp_;

    is_initialized_ = true;
    return;

  } else {
    /**
       Do nothing - and for heaven's sake don't set the measurement pack
                    data again for state vector x_ !!
                    otherwise you will get no Kalman Filter but a value jumper ;-) 
    **/
  }

  Time_Step_ += 1;
  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  // cout << "Delta_t = " << dt << endl;

  Prediction(dt);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(measurement_pack);
  } else {
    UpdateLidar(measurement_pack);
  }


  // Set useUKF = true in case you want to use the Unscented Kalman Filter
  // Set useUKF = false in case you want to use the Extenden Kalman Filter
  bool useUKF = true;

  if (useUKF == true) {
    // I've just moved this conversion for the main version to here
    // the reason is that I wanted to have a common SW for EKF and UKF
    double v   = x_(2);
    double yaw = x_(3);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    averaged_(0) = x_(0);
    averaged_(1) = x_(1);
    averaged_(2) = v1;
    averaged_(3) = v2;

  } else {
    /** Additionally run the the Extended Kalman Filter 
        in order to be able to switch from EKF to UKF
    **/
    ekf_.ProcessMeasurement(measurement_pack);

    averaged_(0) = ekf_.ekf_.x_(0);
    averaged_(1) = ekf_.ekf_.x_(1);
    averaged_(2) = ekf_.ekf_.x_(2);
    averaged_(3) = ekf_.ekf_.x_(3);
  }

  // cout << tools.printVector("Result status vector: " , averaged_) << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  AugmentedSigmaPoints();
  SigmaPointPrediction(delta_t);
  PredictMeanAndCovariance();
}




void UKF::AugmentedSigmaPoints() {
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create augmented mean state
  for (int i = 0; i < n_x_; i++) {
    x_aug(i) = x_(i);
  }

  x_aug(n_aug_-1) = 0;
  x_aug(n_aug_-2) = 0;


  // create augmented covariance matrix
  P_aug.setZero();

  for (int i = 0; i < n_x_; i++) {
    for (int j = 0; j < n_x_; j++) {
      P_aug(i, j) = P_(i, j);
    }
  }

  P_aug(n_aug_-1, n_aug_-1) = pow(std_yawdd_, 2);
  P_aug(n_aug_-2, n_aug_-2) = pow(std_a_, 2);

  // create square root matrix
  MatrixXd sigma_root = ((lambda_ + n_aug_) * P_aug).llt().matrixL();

  // create augmented sigma points
  Xsig_aug_.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i+1) = x_aug + sigma_root.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sigma_root.col(i);
  }

  // print result
  // std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
}

void UKF::SigmaPointPrediction(double delta_t) {
  // create matrix with predicted sigma points as columns

  // if derivative of psi == 0 --> we need different calculations!
  // avoid calculation by zero!

  for (int i = 0; i <  2 * n_aug_ + 1; i++) {
    double px = Xsig_aug_(0, i);
    double py = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double psi = Xsig_aug_(3, i);
    double psi_dot = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_psi = Xsig_aug_(6, i);

    if (fabs(psi_dot) < 0.001) {
      // px
      Xsig_pred_(0, i) = px;
      Xsig_pred_(0, i) += v*cos(psi)*delta_t;
      Xsig_pred_(0, i) += (1.0/2.0)*pow(delta_t, 2.0)*cos(psi)*nu_a;

      // py
      Xsig_pred_(1, i) = py;
      Xsig_pred_(1, i) += v*sin(psi)*delta_t;
      Xsig_pred_(1, i) += (1.0/2.0)*pow(delta_t, 2.0)*sin(psi)*nu_a;

      // v
      Xsig_pred_(2, i) = v;
      Xsig_pred_(2, i) +=  0.0;
      Xsig_pred_(2, i) +=  delta_t * nu_a;

      // psi
      Xsig_pred_(3, i) = psi;
      Xsig_pred_(3, i) += psi_dot*delta_t;
      Xsig_pred_(3, i) += (1.0/2.0)*pow(delta_t, 2.0)*nu_psi;

      // derivative psi
      Xsig_pred_(4, i) = psi_dot;
      Xsig_pred_(4, i) +=  0.0;
      Xsig_pred_(4, i) +=  delta_t * nu_psi;
    } else {
      // px
      Xsig_pred_(0, i) = px;
      Xsig_pred_(0, i) += (v/psi_dot)*(sin(psi + psi_dot*delta_t) - sin(psi));
      Xsig_pred_(0, i) += (1.0/2.0)*pow(delta_t, 2.0)*cos(psi)*nu_a;

      // py
      Xsig_pred_(1, i) = py;
      Xsig_pred_(1, i) += (v/psi_dot)*(-cos(psi + psi_dot*delta_t) + cos(psi));
      Xsig_pred_(1, i) += (1.0/2.0)*pow(delta_t, 2.0)*sin(psi)*nu_a;

      // v
      Xsig_pred_(2, i) = v;
      Xsig_pred_(2, i) +=  0.0;
      Xsig_pred_(2, i) +=  delta_t * nu_a;

      // psi
      Xsig_pred_(3, i) = psi;
      Xsig_pred_(3, i) += psi_dot*delta_t;
      Xsig_pred_(3, i) += (1.0/2.0)*pow(delta_t, 2.0)*nu_psi;

      // derivative psi
      Xsig_pred_(4, i) = psi_dot;
      Xsig_pred_(4, i) +=  0.0;
      Xsig_pred_(4, i) +=  delta_t * nu_psi;
    }
  }

  // print result
  // std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
}

void UKF::PredictMeanAndCovariance() {
  // initialize vector for predicted state
  x_.fill(0.0);

  // initialize covariance matrix for prediction
  P_.fill(0.0);


  // predicted state mean
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)< -M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  // print result
  // std::cout << "Predicted state" << std::endl;
  // std::cout << x << std::endl;
  // std::cout << "Predicted covariance matrix" << std::endl;
  // std::cout << P << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // transform sigma points into measurement space
  Zsig.fill(0.0);

  for (int i = 0; i <  2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double psi = Xsig_pred_(3, i);
    double psi_dot = Xsig_pred_(4, i);

    double rho = sqrt(pow(px, 2.0) + pow(py, 2.0));
    double phi = atan2(py, px);
    double rho_dot = (px*cos(psi)*v + py*sin(psi)*v)/rho;

    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    while (Zsig(1, i)> M_PI) Zsig(1, i)-=2.*M_PI;
    while (Zsig(1, i)< -M_PI) Zsig(1, i)+=2.*M_PI;

    Zsig(2, i) = rho_dot;
  }

  // calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // calculate measurement covariance matrix S
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (x_diff(1)> M_PI) x_diff(1)-=2.*M_PI;
    while (x_diff(1)< -M_PI) x_diff(1)+=2.*M_PI;

    S = S + weights_(i) * x_diff * x_diff.transpose();
  }

  S(0, 0) += pow(std_radr_, 2.0);
  S(1, 1) += pow(std_radphi_, 2.0);
  S(2, 2) += pow(std_radrd_, 2.0);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // z = rho, phi, rho_dot
    VectorXd z_diff = Zsig.col(i)-z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)< -M_PI) z_diff(1)+=2.*M_PI;

    // x = px,py,v,psi,psi_dot
    VectorXd x_diff = Xsig_pred_.col(i)-x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)< -M_PI) x_diff(3)+=2.*M_PI;
    Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
  }


  // calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  // update state mean and covariance matrix
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)< -M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  // calculate the radar NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  if (NIS_radar_ > NIS_radar_threshold_) {
    NIS_radar_over_threshold_ += 1;
    cout << "NIS_radar_over_threshold_ = " << NIS_radar_over_threshold_
         << " which is in percent: "
         << 100.0*(NIS_radar_over_threshold_/Time_Step_) << endl;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the lidar NIS.
  */

  // set measurement dimension, lidar can measure px and py
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // transform sigma points into measurement space
  Zsig.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 sigma points
    // measurement model
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  // calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  // calculate measurement covariance matrix S
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * x_diff * x_diff.transpose();
  }

  S(0, 0) += pow(std_laspx_, 2.0);
  S(1, 1) += pow(std_laspy_, 2.0);


  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i)-z_pred;

    VectorXd x_diff = Xsig_pred_.col(i)-x_;
    Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  // update state mean and covariance matrix
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

  // laser NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  if (NIS_laser_ > NIS_laser_threshold_) {
    NIS_laser_over_threshold_ += 1;

    cout << "NIS_laser_over_threshold_ = " << NIS_laser_over_threshold_
         << " which is in percent: "
         << 100.0*(NIS_laser_over_threshold_/Time_Step_) << endl;
  }

  // cout << "NIS_laser_ = " << NIS_laser_ << endl;
}
