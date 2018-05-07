#pragma once
#ifndef __Environmentobj__
#define __Environmentobj__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;

class Env_obj
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Env_obj();
	Env_obj(const std::shared_ptr<Shape> shape);
	virtual ~Env_obj();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;

	Eigen::Vector3d dim;
	Eigen::Vector3d x;  // position

private:
	const std::shared_ptr<Shape> env_obj;
};

#endif
#pragma once
