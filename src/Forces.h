#pragma once
#ifndef __Forces__
#define __Forces__

#include <vector>
#include <memory>
#include <string>

#include "Cloth.h"

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

class Forces
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		Forces() {};
	virtual ~Forces() {};

	Eigen::VectorXd f;
	Eigen::SparseMatrix<double> M;
	Eigen::SparseMatrix<double> MDK;

	void fill(const Mesh& mesh, const Material& mat, const Eigen::Vector3d& grav, double h);

#ifdef EOLC_ONLINE
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
#endif // EOLC_ONLINE
};

#endif