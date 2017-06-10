#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include <string>

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
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

 /**
   * Compute the Jacobian Matrix
   * @param x_state The state
   */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

 /**
   * Convert Cartesian Coordinates into Polar
   * @param x_state The state
   */
  VectorXd Cartesian2Polar(const VectorXd& x_state);

 /**
   * Convert Polar Coordinates into Cartesian
   * @param x_state The state
   */
  VectorXd Polar2Cartesian(const VectorXd& x_state);

 /**
   * Pretty print Vector
   * @param x_state The vector
   */
  std::string printVector(std::string head, const VectorXd& x_state);
 
 /**
   * Pretty print Vector
   * @param x_state The vector
   */  
  std::string printMatrix(std::string head, const MatrixXd& x_state);
};

#endif /* TOOLS_H_ */
