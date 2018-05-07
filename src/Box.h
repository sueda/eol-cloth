#pragma once
#ifndef __Box__
#define __Box__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include "BrenderManager.h"
#include "Brenderable.h"

// ArcSim
#include "mesh.hpp"

class Shape;
class Program;
class MatrixStack;
class Rigid;

class Box : public Brenderable
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Mesh mesh;
	std::vector<Mesh*> mesh_vec;
	std::shared_ptr<Rigid> rigid;

	Box();
	Box(const std::shared_ptr<Shape> shape);
	virtual ~Box();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void init();
	void step(double h);

	double thresh;

	Eigen::Vector3d dim;
	Eigen::Vector3d rot;
	Eigen::Vector3d x;  // position
	Eigen::Matrix4d E1;
	Eigen::Matrix4d E1inv;
	Eigen::VectorXd v;
	Eigen::MatrixXd adjoint;

	// Export
	int getBrenderCount() const;
	std::vector<std::string> getBrenderNames() const;
	void exportBrender(std::vector< std::shared_ptr< std::ofstream > > outfiles) const;

private:
	const std::shared_ptr<Shape> box;
};

#endif
