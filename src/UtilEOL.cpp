#include "UtilEOL.h"

using namespace std;
using namespace Eigen;

MatrixXd deform_grad(const Face *f) {
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