#pragma once
#ifndef __GeneralizedSolver__
#define __GeneralizedSolver__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class GeneralizedSolver
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		GeneralizedSolver();
	virtual ~GeneralizedSolver() {};

	enum Solver {
		NoSolver = 0, Mosek = 1, Gurobi = 2
	};

	int whichSolver;

	bool velocitySolve(const bool& fixedPoints, const bool& collisions,
		const Eigen::SparseMatrix<double>& MDK, const Eigen::VectorXd& b,
		const Eigen::SparseMatrix<double>& Aeq, const Eigen::VectorXd& beq,
		const Eigen::SparseMatrix<double>& Aineq, const Eigen::VectorXd& bineq,
		Eigen::VectorXd& v);

};

#endif