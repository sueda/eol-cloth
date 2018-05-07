#pragma once
#ifndef _RIGID_H_
#define _RIGID_H_

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 3, 6> Matrix3x6d;
typedef Eigen::Matrix<double, 6, 3> Matrix6x3d;

class Rigid
{
private:
	Rigid();

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	static Eigen::Matrix4d inverse(const Eigen::Matrix4d &E);
	static Matrix3x6d gamma(const Eigen::Vector3d &r); 
	static Matrix6d adjoint(const Eigen::Matrix4d &E);
	static Eigen::Matrix3d bracket3(const Eigen::Vector3d &a);
	static Eigen::Matrix4d bracket6(const Vector6d &a);
	static Eigen::Vector3d unbracket3(const Eigen::Matrix3d &A);
	static Vector6d unbracket6(const Eigen::Matrix4d &A);
	static Eigen::Matrix4d integrate(const Eigen::Matrix4d &E0, const Eigen::VectorXd &phi, double h);
	static Vector6d log(const Eigen::Matrix4d &A);
};

#endif
