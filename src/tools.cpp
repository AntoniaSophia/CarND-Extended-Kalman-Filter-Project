#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth){

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean

	rmse = rmse/estimations.size();

  cout << rmse << endl;

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3, 4);

	Hj << 0,0,0,0,0,0,0,0,0,0,0,0;

	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	//TODO: YOUR CODE HERE 
	double psum = pow(px, 2) + pow(py, 2);

	//check division by zero
	if (fabs(psum) < 0.0001) {
		std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
		return Hj;
	}

	double pmean = sqrt(psum);
	 
	Hj(0, 0) = px / pmean;
	Hj(0, 1) = py / pmean;

	Hj(1, 0) = -1 * (py / psum);
	Hj(1, 1) =  px / psum;

	Hj(2, 0) = py *(vx*py - vy*px) / (pow(psum, (3 / 2)));
	Hj(2, 1) = px *(vy*px - vx*py) / (pow(psum, (3 / 2)));

	Hj(2, 2) = Hj(0, 0);
	Hj(2, 3) = Hj(0, 1);

	return Hj;
}

VectorXd Tools::Cartesian2Polar(const VectorXd& x_state) {
  // result is map prediction in polar coordinates rho, phi, rho_dot
  VectorXd h_polar(3);
  h_polar << 0, 0, 0;

  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  double epsilon = 0.00001;

  if (fabs(px) < epsilon && fabs(py) < epsilon) {
  	cout << "EPSILON!!!!" << endl;
    return h_polar;
  }

  h_polar(0) = sqrt(pow(px, 2) + pow(py, 2));

  if(fabs(px) < 0.001) {
    px = 0.001;
  }

  h_polar(1) = atan2(py,px);
	if (fabs(h_polar(1)) < 0.001){
    h_polar(1) = 0.001;
  }
	//cout << "Tools::Cartesian2Polar --> Polar angle phi = " << h_polar(1) << endl;
  h_polar(2) = (px*vx + py*vy)/h_polar(0);
}

VectorXd Tools::Polar2Cartesian(const VectorXd& x_state) {
  VectorXd cartesian(4);
  cartesian << 0, 0, 0, 0;

  double rho = x_state(0);
  double phi = x_state(1);
	
	// cout << "Tools::Polar2Cartesian --> Polar angle phi = " << phi << endl;

  double rho_dot = x_state(2);

  double x_ein = cos(phi);
  double y_ein = sin(phi);
  cartesian(0) = rho*x_ein;
  cartesian(1) = rho*y_ein;
  cartesian(2) = rho_dot*x_ein;
  cartesian(3) = rho_dot*y_ein;

  return cartesian;
}

std::string Tools::printVector(std::string head, const VectorXd& x_state){
	stringstream  result;
	
	result << head;
	for (int i=0; i < x_state.size(); ++i){
		result << x_state(i) << " , ";
	}
	return result.str();
}

std::string Tools::printMatrix(std::string head, const MatrixXd& x_state){
	stringstream  result;
	
	result << head << ": ";

	result << x_state;

	return result.str();
}