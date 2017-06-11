/* Copyright 2017 Antonia Reiter */
/* no obligations - feel free to copy/reuse/modify as you like*/

#ifndef SRC_TOOLS_H_
#define SRC_TOOLS_H_
#include <string>
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
 public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                         const vector<VectorXd> &ground_truth);

 /**
   * Compute the Jacobian Matrix
   * @param x_state: the Vector to calculate the Jacobian Matrix os
   */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

 /**
   * Convert Cartesian Coordinates into Polar
   * @param x_state: the Vector to convert to Polar coordinates
   */
  VectorXd Cartesian2Polar(const VectorXd& x_state);

 /**
   * Convert Polar Coordinates into Cartesian
   * @param x_state: the Vector to convert to Cartesian coordinates
   */
  VectorXd Polar2Cartesian(const VectorXd& x_state);

 /**
   * Pretty print Vector
   * @param x_state: the vector to be pretty printed
   */
  std::string printVector(std::string head, const VectorXd& x_state);

 /**
   * Pretty print Matrix object
   * @param x_state: the Matrix to be pretty printed
   */  
  std::string printMatrix(std::string head, const MatrixXd& x_state);
};

#endif /* SRC_TOOLS_H_ */
