#include "Rigid.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

#define EIGEN_DONT_ALIGN_STATICALLY
#define EIGEN_NO_STATIC_ASSERT
using namespace Eigen;

Matrix4d Rigid::inverse(const Matrix4d &E)
{
	Matrix4d Einv = Matrix4d::Identity();
	Matrix3d R = E.block<3,3>(0,0);
	Vector3d p = E.block<3,1>(0,3);
	Matrix3d Rt = R.transpose();
	Einv.block<3,3>(0,0) = Rt;
	Einv.block<3,1>(0,3) = -Rt * p;
	return Einv;
}

Matrix3x6d Rigid::gamma(const Eigen::Vector3d &r)
{
	Matrix3x6d G = Matrix3x6d::Zero();
	G.block<3,3>(0,0) = bracket3(r).transpose();
	G.block<3,3>(0,3) = Matrix3d::Identity();
	return G;
}

Matrix6d Rigid::adjoint(const Matrix4d &E)
{
	Matrix6d Ad = Matrix6d::Zero();
	Matrix3d R = E.block<3,3>(0,0);
	Vector3d p = E.block<3,1>(0,3);
	Ad.block(0,0,3,3) = R;
	Ad.block(3,0,3,3) = bracket3(p) * R;
	Ad.block(3,3,3,3) = R;
	return Ad;
}

Matrix3d Rigid::bracket3(const Vector3d &a)
{
	Matrix3d A = Matrix3d::Zero();
	A(0,1) = -a(2);
	A(0,2) =  a(1);
	A(1,0) =  a(2);
	A(1,2) = -a(0);
	A(2,0) = -a(1);
	A(2,1) =  a(0);
	return A;
}

Matrix4d Rigid::bracket6(const Vector6d &a)
{
	Matrix4d A = Matrix4d::Zero();
	A.block<3,3>(0,0) = bracket3(a.segment<3>(0));
	A.block<3,1>(0,3) = a.segment<3>(3);
	return A;
}

Vector3d Rigid::unbracket3(const Matrix3d &A)
{
	Vector3d a;
	a(0) = A(2,1);
	a(1) = A(0,2);
	a(2) = A(1,0);
	return a;
}

Vector6d Rigid::unbracket6(const Matrix4d &A)
{
	Vector6d a;
	a.segment<3>(0) = unbracket3(A.block<3,3>(0,0));
	a(3) = A(0,3);
	a(4) = A(1,3);
	a(5) = A(2,3);
	return a;
}

Matrix4d Rigid::integrate(const Matrix4d &E0, const VectorXd &phi, double h)
{
	Matrix3d I = Matrix3d::Identity();
	Vector3d w = phi.segment<3>(0);
	Vector3d v = phi.segment<3>(3);
	Matrix4d phib = Matrix4d::Identity();
	phib.block<3,1>(0,3) = h*v;
	double wlen = w.norm();
	if(wlen > 1e-10) {
		w /= wlen;
		v /= wlen;
		// Rodrigues formula %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		double wX = w(0);
		double wY = w(1);
		double wZ = w(2);
		double c = cos(wlen * h);
		double s = sin(wlen * h);
		double c1 = 1.0 - c;
		Matrix3d R;
		R << c + wX * wX * c1, -wZ * s + wX * wY * c1, wY * s + wX * wZ * c1,
		wZ * s + wX * wY * c1, c + wY * wY * c1, -wX * s + wY * wZ * c1,
		-wY * s + wX * wZ * c1, wX * s + wY * wZ * c1, c + wZ * wZ * c1;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		Matrix3d A = I - R;
		Vector3d cc = w.cross(v);
		Vector3d d = A * cc;
		double wv = w.dot(v);
		Vector3d p = (wv * wlen * h) * w + d;
		phib.block<3,3>(0,0) = R;
		phib.block<3,1>(0,3) = p;
		//cout << phib << endl;
	}
	return E0 * phib;
}

Vector6d Rigid::log(const Eigen::Matrix4d & A)
{
	Matrix3d R = A.block<3,3>(0,0);
	Vector3d p = A.block<3,1>(0,3);
	Vector6d phi;

	double cosTheta = 0.5*(R.trace() - 1);
	double theta = acos(cosTheta);

	if (cosTheta > 0.999999999999)
		theta = 0;
	else if (cosTheta < -0.999999999999)
		theta = M_PI;

	if (abs(theta) < 1e-8)  
	{
		phi.segment<3>(0) << 0.0, 0.0, 0.0;
		phi.segment<3>(3) = p;
	}
	else
	{
		double sinTheta = sin(theta);
		Matrix3d wBracket = theta / (2 * sinTheta)*(R - R.transpose());
		phi.segment<3>(0) = unbracket3(wBracket);
		Matrix3d V = Matrix3d::Identity() + ((1 - cosTheta) / (theta * theta) * wBracket) + ((theta - sinTheta) / (theta * theta * theta) * wBracket*wBracket);
		phi.segment<3>(3) = V.colPivHouseholderQr().solve(p);
	}

	return phi;
}


// #include <iostream>
// using namespace std;
// int main(int argc, char **argv)
// {
// 	Matrix4d E = Matrix4d::Identity();
// 	E <<    -0.1817,    0.7236,    0.6658,    2.7694,
// 			-0.6198,   -0.6100,    0.4938,   -1.3499,
// 			 0.7634,   -0.3230,    0.5594,    3.0349,
// 			      0,         0,         0,    1.0000;
// 	cout << E << endl << endl;
// 	cout << Rigid::inverse(E) << endl << endl;
// 	Vector3d r;
// 	r <<     0.7254, -0.0631, 0.7147;
// 	cout << Rigid::gamma(r) << endl << endl;
// 	cout << Rigid::adjoint(E) << endl << endl;
// 	Vector6d p;
// 	p <<    -0.2050, -0.1241, 1.4897, 1.4090, 1.4172, 0.6715;
// 	cout << Rigid::bracket6(p) << endl << endl;
// 	Vector3d ax;
// 	ax <<    -1.2075, 0.7172, 1.6302;
// 	double a = 0.4889;
// 	cout << Rigid::aaToMat(ax, a) << endl << endl;
// 	cout << Rigid::aaToMat(Vector3d(1, 0, 0), a) << endl << endl;
// 	cout << Rigid::aaToMat(Vector3d(0, 1, 0), a) << endl << endl;
// 	cout << Rigid::aaToMat(Vector3d(0, 0, 1), a) << endl << endl;
// }
