#pragma once
#ifndef __Box__
#define __Box__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include "BrenderManager.h"
#include "Brenderable.h"

#ifdef EOLC_ONLINE
class MatrixStack;
class Program;
#endif // EOLC_ONLINE

class Shape;
class Rigid;

class Box : public Brenderable
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	std::shared_ptr<Rigid> rigid;

	Box(const std::shared_ptr<Shape> shape);
	virtual ~Box();
	void step(double h);

#ifdef EOLC_ONLINE
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void init();
#endif // EOLC_ONLINE

	int num_points;
	int num_edges;

	Eigen::Vector3d dim;
	Eigen::Vector3d rot;
	Eigen::Vector3d x;  // position
	Eigen::Matrix4d E1;
	Eigen::Matrix4d E1inv;
	Eigen::VectorXd v;
	Eigen::MatrixXd adjoint;

	// These are used for constraints
	Eigen::MatrixXd faceNorms;
	Eigen::MatrixXi edgeFaces;


	// Export
	int getBrenderCount() const;
	std::vector<std::string> getBrenderNames() const;
	void exportBrender(std::vector< std::shared_ptr< std::ofstream > > outfiles) const;

private:
	void generateConstraints();

	const std::shared_ptr<Shape> boxShape;
};

#endif
