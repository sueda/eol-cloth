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
class FixedList;

class Constraints
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		Constraints();
	virtual ~Constraints() {};

	bool hasFixed;
	bool hasCollisions;

	Eigen::SparseMatrix<double> Aeq;
	Eigen::SparseMatrix<double> Aineq;

	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > drawAineq;
	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > drawAeq;

	Eigen::VectorXd beq;
	Eigen::VectorXd bineq;

	Eigen::MatrixXd constraintTable;

	void init(const std::shared_ptr<Obstacles> obs);
	void updateTable(const std::shared_ptr<Obstacles> obs);
	void fill(const Mesh& mesh, const std::shared_ptr<Obstacles> obs, const std::shared_ptr<FixedList> fs, double h, const bool& online);

#ifdef EOLC_ONLINE
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
#endif // EOLC_ONLINE
};

#endif