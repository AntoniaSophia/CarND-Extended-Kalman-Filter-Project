/* Copyright 2017 Udacity */
/* no obligations - feel free to copy/reuse/modify as you like*/
#ifndef SRC_MEASUREMENT_PACKAGE_H_
#define SRC_MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"

class MeasurementPackage {
 public:
  double timestamp_;

  enum SensorType {
    LASER,
    RADAR
  } sensor_type_;

  Eigen::VectorXd raw_measurements_;
};

#endif  // SRC_MEASUREMENT_PACKAGE_H_
