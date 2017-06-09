#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	if (estimations.size() == 0) {
		return rmse;
	}

	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() != ground_truth.size()) {
		return rmse;
	}

	//accumulate squared residuals
	 

	//accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3, 4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE 
	float psum = pow(px, 2) + pow(py, 2);

	//check division by zero
	if (fabs(psum) < 0.0001) {
		std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
		return Hj;
	}

	float pmean = sqrt(psum);
	 
	Hj(0, 0) = px / pmean;
	Hj(0, 1) = py / pmean;

	Hj(1, 0) = -py / psum;
	Hj(1, 1) =  px / psum;

	Hj(2, 0) = py *(vx*py - vy*px) / (pow(psum, (3 / 2)));
	Hj(2, 1) = px *(vy*px - vx*py) / (pow(psum, (3 / 2)));

	Hj(2, 2) = Hj(0, 0);
	Hj(2, 3) = Hj(0, 1);

	return Hj;
}
