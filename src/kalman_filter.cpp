/* Copyright 2017 Antonia Reiter */

#include "kalman_filter.h"
#include <iostream>
#include <math.h>

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

void KalmanFilter::Predict() {
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;

}

std::string KalmanFilter::Update(const VectorXd &z) {
  std::string result="";
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  int16_t x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  return result;
}

std::string KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations

    One other important point when calculating y with radar sensor 
    data: the second value in the polar coordinate vector is the 
    angle ϕ. You'll need to make sure to normalize ϕ in the y vector 
    so that its angle is between −pi and pi; in other words, add or 
    subtract 2pi from ϕ until it is between −pi and pi.
  
  */
  std::string result="";

  // VectorXd z_pred = H_ * x_;

  // This is just a sanity check wether conversion from Cartesian to Polar is working as expected
  // result += tools.printVector("\n# x_ in Cartesian:",x_);
  // result += tools.printVector("\n# x_ in Polar:",tools.Cartesian2Polar(x_));
  // result += tools.printVector("\n# x_ again in Cartesian:",tools.Polar2Cartesian(tools.Cartesian2Polar(x_)));

  VectorXd ypolar = z - tools.Cartesian2Polar(x_);
  //result += tools.printVector("\n# Y-Vector in Polar:",ypolar);

  if (ypolar(1)>M_PI){
    while(ypolar(1)>M_PI){
      ypolar(1) -= 2*M_PI;

    }
  }else if(ypolar(1)<-M_PI){
    while(ypolar(1)<-M_PI){
      ypolar(1) += 2*M_PI;
    }
  } 

  //result += tools.printVector("\n# Y-Vector in Polar (corrected) :",y);

  //result += tools.printMatrix("\n# H_ :",H_)+"\n";
  MatrixXd Ht = H_.transpose();
  //result += tools.printMatrix("\n# Ht :",Ht)+"\n";
  MatrixXd S = H_ * P_ * Ht + R_;
  //result += tools.printMatrix("\n# S :",S)+"\n";
  MatrixXd Si = S.inverse();
 // result += tools.printMatrix("\n# Si :",Si)+"\n";
  MatrixXd PHt = P_ * Ht;
  //result += tools.printMatrix("\n# PHt :",PHt)+"\n";
  MatrixXd K = PHt * Si;
  //result += tools.printMatrix("\n# K :",K)+"\n";

  // new estimate
  //VectorXd y = tools.Polar2Cartesian(ypolar);

  x_ = x_ + (K * ypolar);
  //x_ = x_ + tools.Polar2Cartesian(K * ypolar);

  int16_t x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;  
  return result;
}