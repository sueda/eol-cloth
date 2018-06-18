#pragma once
#ifndef __Points__
#define __Points__

#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#ifdef EOLC_ONLINE
class MatrixStack;
class Program;
#endif // EOLC_ONLINE

class Points
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		Points() : num_points(0) {};
	virtual ~Points() {};

	int num_points;
	Eigen::MatrixXd pxyz;
	Eigen::MatrixXd norms;

#ifdef EOLC_ONLINE
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
#endif // EOLC_ONLINE

private:

};

#endif