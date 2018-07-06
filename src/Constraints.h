#pragma once
#ifndef __Constraints__
#define __Constraints__

#include <vector>
#include <memory>
#include <string>

#include "external\ArcSim\mesh.hpp"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#ifdef EOLC_ONLINE
class MatrixStack;
class Program;
#endif // EOLC_ONLINE

//class Mesh;
class Obstacles;

class Constraints
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		Constraints();
	virtual ~Constraints() {};

	Eigen::SparseMatrix<double> Aeq;
	Eigen::SparseMatrix<double> Aineq;

	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > drawAineq;

	Eigen::VectorXd beq;
	Eigen::VectorXd bineq;

	Eigen::MatrixXd constraintTable;

	void init(const std::shared_ptr<Obstacles> obs);
	void updateTable(const std::shared_ptr<Obstacles> obs);
	void fill(const Mesh& mesh, const std::shared_ptr<Obstacles> obs, double h);

#ifdef EOLC_ONLINE
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
#endif // EOLC_ONLINE
};

#endif