#pragma once
#ifndef __Points__
#define __Points__

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Points
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		Points() : num_points(0) {};
	virtual ~Points() {};

	int num_points;
	Eigen::MatrixXd pxyz;
	Eigen::MatrixXd norms;

private:

};

#endif