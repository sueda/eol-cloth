#pragma once
#ifndef __QuadProg_H__
#define __QuadProg_H__

#include <Eigen/Dense>
#include <Eigen/Sparse>

// Calling convention:
// https://www.mathworks.com/help/optim/ug/quadprog.html
class QuadProg
{
public:
	QuadProg() {}
	virtual ~QuadProg() {}

	virtual void setNumberOfVariables(int numVars) = 0;
	virtual void setNumberOfInequalities(int numIneq) = 0;
	virtual void setNumberOfEqualities(int numEq) = 0;

	virtual void setObjectiveMatrix(const Eigen::SparseMatrix<double> & mat) = 0;
	virtual void setObjectiveVector(const Eigen::VectorXd & vector) = 0;
	virtual void setObjectiveConstant(double constant) = 0;

	virtual void setLowerVariableBound(const Eigen::VectorXd & bounds) = 0;
	virtual void setUpperVariableBound(const Eigen::VectorXd & bounds) = 0;

	virtual void setInequalityMatrix(const Eigen::SparseMatrix<double> & mat) = 0;
	virtual void setInequalityVector(const Eigen::VectorXd & vector) = 0;

	virtual void setEqualityMatrix(const Eigen::SparseMatrix<double> & mat) = 0;
	virtual void setEqualityVector(const Eigen::VectorXd & vector) = 0;

	virtual bool solve() = 0;

	virtual Eigen::VectorXd getPrimalSolution() = 0;
	virtual Eigen::VectorXd getDualInequality() = 0;
	virtual Eigen::VectorXd getDualEquality() = 0;
	virtual Eigen::VectorXd getDualLower() = 0;
	virtual Eigen::VectorXd getDualUpper() = 0;
};

#endif
