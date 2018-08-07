#pragma once
#ifndef __FixedList__
#define __FixedList__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class FixedList
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		FixedList() : when(0.0) {
		c1 = Eigen::VectorXd::Zero(6);
		c2 = Eigen::VectorXd::Zero(6);
		c3 = Eigen::VectorXd::Zero(6);
		c4 = Eigen::VectorXd::Zero(6);

		c1c2 = Eigen::VectorXd::Zero(6);
		c2c4 = Eigen::VectorXd::Zero(6);
		c4c3 = Eigen::VectorXd::Zero(6);
		c3c1 = Eigen::VectorXd::Zero(6);
	};
	virtual ~FixedList() {};

	double when;

	int c1i, c2i, c3i, c4i;

	Eigen::VectorXd c1;
	Eigen::VectorXd c2;
	Eigen::VectorXd c3;
	Eigen::VectorXd c4;

	// I didn't impliment this here yet
	Eigen::VectorXd c1c2;
	Eigen::VectorXd c2c4;
	Eigen::VectorXd c4c3;
	Eigen::VectorXd c3c1;

};

#endif