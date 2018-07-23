#include "UtilEOL.h"

#include "external\ArcSim\geometry.hpp"

#include <iostream>

#include <Eigen\SVD>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

MatrixXd deform_grad(const Face *f)
{
	MatrixXd Dx(3, 2);
	Dx(0, 0) = f->v[1]->node->x[0] - f->v[0]->node->x[0];
	Dx(0, 1) = f->v[2]->node->x[0] - f->v[0]->node->x[0];
	Dx(1, 0) = f->v[1]->node->x[1] - f->v[0]->node->x[1];
	Dx(1, 1) = f->v[2]->node->x[1] - f->v[0]->node->x[1];
	Dx(2, 0) = f->v[1]->node->x[2] - f->v[0]->node->x[2];
	Dx(2, 1) = f->v[2]->node->x[2] - f->v[0]->node->x[2];
	Matrix2d DX;
	DX(0, 0) = f->v[1]->u[0] - f->v[0]->u[0];
	DX(0, 1) = f->v[2]->u[0] - f->v[0]->u[0];
	DX(1, 0) = f->v[1]->u[1] - f->v[0]->u[1];
	DX(1, 1) = f->v[2]->u[1] - f->v[0]->u[1];
	return Dx * DX.inverse();
}

MatrixXd deform_grad_v(const Vert* v)
{
	double tot_ang = 0.0;
	Matrix3d Q = Matrix3d::Zero();
	Matrix2d P = Matrix2d::Zero();

	for (int f = 0; f < v->adjf.size(); f++) {
		Face* face = v->adjf[f];

		MatrixXd F = deform_grad(face);
		JacobiSVD<MatrixXd> svd(F, ComputeFullU | ComputeFullV);

		MatrixXd V3(2, 3);
		V3 << svd.matrixV(), Vector2d::Zero();

		MatrixXd Qx = svd.matrixU() * V3.transpose();
		Matrix3d Qxrot;
		Qxrot << Qx.col(0), Qx.col(1), Qx.block<3,1>(0,0).cross(Qx.block<3,1>(0,1));

		Matrix2d S2;
		S2 << svd.singularValues(), svd.singularValues();

		Matrix2d Px = svd.matrixV() * S2 * svd.matrixV().transpose();

		Q += incedent_angle(v, face) * Qxrot.log();
		
		P += incedent_angle(v, face) * Px;
		
		tot_ang += incedent_angle(v, face);
	}
	
	Q /= tot_ang;
	Q = Q.exp();

	P /= tot_ang;

	return Q.block<3, 2>(0, 0) * P;
}