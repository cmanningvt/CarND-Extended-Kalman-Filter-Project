#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  // initialize variables
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	if ( estimations.size() == 0 ) {
        cout << "CalculateRMSE () - Error - Estimations vector size must be nonzero" << endl;
        return rmse;
    }
	//  * the estimation vector size should equal ground truth vector size
	if ( estimations.size() != ground_truth.size() ) {
        cout << "CalculateRMSE () - Error - Input vectors must have same size" << endl;
        return rmse;
    }

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        VectorXd residual = estimations[i] - ground_truth[i];
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
/*
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  // initialize variables
	MatrixXd Hj(3,4);
  
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	if ( px == 0 && py == 0 ) {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }
	
	//compute the Jacobian matrix
    Hj(0, 0) =  px                     / sqrt( px*px + py*py );
    Hj(0, 1) =  py                     / sqrt( px*px + py*py );
	Hj(0, 2) =  0;
	Hj(0, 3) =	0;
    Hj(1, 0) = -py                     /     ( px*px + py*py );
    Hj(1, 1) =  px                     /     ( px*px + py*py );
	Hj(1, 2) =  0;
	Hj(1, 3) =  0;
    Hj(2, 0) =  py * ( vx*py - vy*px ) /  pow( px*px + py*py , 3/2 );
    Hj(2, 1) =  px * ( vy*px - vx*py ) /  pow( px*px + py*py , 3/2 );
    Hj(2, 2) =  px                     / sqrt( px*px + py*py );
    Hj(2, 3) =  py                     / sqrt( px*px + py*py );

	return Hj;
}
*/

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	// Initiate Jacobian Matrix
	MatrixXd Hj(3, 4);

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px + py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy*px) / c3, px*(px*vy - py*vx) / c3, px / c2, py / c2;

	return Hj;
}