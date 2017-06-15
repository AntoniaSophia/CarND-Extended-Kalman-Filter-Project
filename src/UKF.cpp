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

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  ##################Antonia: must be std_a_ or std_yawdd_ or both###################
  */

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
  VectorXd weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i < 2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  ///* state covariance matrix
  P_ = MatrixXd(n_x_ , n_x_);
  P_.fill(0.0);
  //##########################################################################
  //TODO: initialize P_ !!
  //##########################################################################

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);


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
    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    x_ << 0, 0.5, 0, 0, 0;

  //##########################################################################
  //TODO: initialize P_ !!
  //##########################################################################

    is_initialized_ = true;
  } else {
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);
      double phi_dot = measurement_pack.raw_measurements_(2);

      double px = rho*cos(phi);
      double py = rho*sin(phi);
      double v = phi_dot;

      x_ << px , py, v , 0, 0;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      double px = measurement_pack.raw_measurements_(0);
      double py = measurement_pack.raw_measurements_(1);

      x_ << px, py, 0, 0, 0;
    }
  }

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  Prediction(dt);

  cout << "Hello world!" << endl;
  exit(0);


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

  MatrixXd Xsig_aug  = MatrixXd(2*n_x_ +1, n_x_);
  MatrixXd Xsig_pred = MatrixXd(2*n_x_ +1, n_x_);
  // GenerateSigmaPoints(&Xsig);
  AugmentedSigmaPoints(&Xsig_aug);
  SigmaPointPrediction(delta_t, Xsig_aug, &Xsig_pred);
  PredictMeanAndCovariance(Xsig_pred, &x_, &P_);
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
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */
}



void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  // create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  MatrixXd sigma_root = ((lambda_ + n_x_) * P_).llt().matrixL();

  // calculate square root of P
  // MatrixXd A = P.llt().matrixL();
  Xsig.col(0) = x_;
  for (int i = 0; i < n_x_; i++) {
    Xsig.col(i+1) = x_ + sigma_root.col(i);
    Xsig.col(i+1+n_x_) = x_ - sigma_root.col(i);
  }


  cout << Xsig << endl;

  // set return value
  *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

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
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sigma_root.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sigma_root.col(i);
  }

  // print result
  // std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  // set return value
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(double delta_t, MatrixXd Xsig_aug, MatrixXd* Xsig_out) {
  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // if derivative of psi == 0 --> we need different calculations!
  // avoid calculation by zero!
  if (fabs(Xsig_aug(0, n_x_-1)) < 0.00001) {
    for (int i = 0; i <  2 * n_aug_ + 1; i++) {
      double px = Xsig_aug(0, i);
      double py = Xsig_aug(1, i);
      double v = Xsig_aug(2, i);
      double psi = Xsig_aug(3, i);
      double psi_dot = Xsig_aug(4, i);
      double nu_a = Xsig_aug(5, i);
      double nu_psi = Xsig_aug(6, i);

      // px
      Xsig_pred(0, i) = px;
      Xsig_pred(0, i) += v*cos(psi)*delta_t;
      Xsig_pred(0, i) += (1.0/2.0)*pow(delta_t, 2.0)*cos(psi)*nu_a;

      // py
      Xsig_pred(1, i) = py;
      Xsig_pred(1, i) += v*sin(psi)*delta_t;
      Xsig_pred(1, i) += (1.0/2.0)*pow(delta_t, 2.0)*sin(psi)*nu_a;

      // v
      Xsig_pred(2, i) = v;
      Xsig_pred(2, i) +=  0.0;
      Xsig_pred(2, i) +=  delta_t * nu_a;

      // psi
      Xsig_pred(3, i) = psi;
      Xsig_pred(3, i) += psi_dot*delta_t;
      Xsig_pred(3, i) += (1.0/2.0)*pow(delta_t, 2.0)*nu_psi;

      // derivative psi
      Xsig_pred(4, i) = psi_dot;
      Xsig_pred(4, i) +=  0.0;
      Xsig_pred(4, i) +=  delta_t * nu_psi;
    }

  } else {
    for (int i = 0; i <  2 * n_aug_ + 1; i++) {
      double px = Xsig_aug(0, i);
      double py = Xsig_aug(1, i);
      double v = Xsig_aug(2, i);
      double psi = Xsig_aug(3, i);
      double psi_dot = Xsig_aug(4, i);
      double nu_a = Xsig_aug(5, i);
      double nu_psi = Xsig_aug(6, i);

      // px
      Xsig_pred(0, i) = px;
      Xsig_pred(0, i) += (v/psi_dot) * (sin(psi + psi_dot*delta_t) - sin(psi));
      Xsig_pred(0, i) += (1.0/2.0)*pow(delta_t, 2.0)*cos(psi)*nu_a;

      // py
      Xsig_pred(1, i) = py;
      Xsig_pred(1, i) += (v/psi_dot) * (-cos(psi + psi_dot*delta_t) + cos(psi));
      Xsig_pred(1, i) += (1.0/2.0)*pow(delta_t, 2.0)*sin(psi)*nu_a;

      // v
      Xsig_pred(2, i) = v;
      Xsig_pred(2, i) +=  0.0;
      Xsig_pred(2, i) +=  delta_t * nu_a;

      // psi
      Xsig_pred(3, i) = psi;
      Xsig_pred(3, i) += psi_dot*delta_t;
      Xsig_pred(3, i) += (1.0/2.0)*pow(delta_t, 2.0)*nu_psi;

      // derivative psi
      Xsig_pred(4, i) = psi_dot;
      Xsig_pred(4, i) +=  0.0;
      Xsig_pred(4, i) +=  delta_t * nu_psi;
    }
  }

  // print result
  // std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  // set output result
  *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(MatrixXd Xsig_pred, VectorXd* x_out, MatrixXd* P_out) {
  
  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);

  // predicted state mean
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred.col(i);
  }

  // predicted state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)< -M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose();
  }

  // print result
  // std::cout << "Predicted state" << std::endl;
  // std::cout << x << std::endl;
  // std::cout << "Predicted covariance matrix" << std::endl;
  // std::cout << P << std::endl;

  // set output result
  *x_out = x;
  *P_out = P;
}