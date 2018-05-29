#pragma once
#ifndef __Constraints__
#define __Constraints__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Obstacles;

class Constraints
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		Constraints();
	virtual ~Constraints() {};

	std::vector<Eigen::VectorXd> constraintTable;

	void init(std::shared_ptr<Obstacles> obs);
	void update(std::shared_ptr<Obstacles> obs);
};

#endif