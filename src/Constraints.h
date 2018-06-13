#pragma once
#ifndef __Constraints__
#define __Constraints__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Mesh;
class Obstacles;

class Constraints
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		Constraints();
	virtual ~Constraints() {};

	Eigen::SparseMatrix<double> Aeq;
	Eigen::SparseMatrix<double> Aineq;

	Eigen::SparseVector<double> beq;
	Eigen::SparseVector<double> bineq;

	std::vector<Eigen::VectorXd> constraintTable;
	std::vector<Eigen::Vector2i> obsTable;

	void init(const std::shared_ptr<Obstacles> obs);
	void updateTable(const std::shared_ptr<Obstacles> obs);
	void fill(const Mesh& mesh, const std::shared_ptr<Obstacles> obs);
};

#endif