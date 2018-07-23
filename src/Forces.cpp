#include "Forces.h"

#ifdef EOLC_ONLINE
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "glm/ext.hpp"

#include "online/Program.h"
#include "online/MatrixStack.h"
#endif // EOLC_ONLINE

#include "conversions.h"
#include "UtilEOL.h"
#include "ComputeBending.h"
#include "ComputeInertial.h"
#include "ComputeMembrane.h"
#include "ComputeInertiaEOL.h"

#include "external\ArcSim\util.hpp"

#include <iostream>
#include <utility>

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;

Matrix2d poldec(const Matrix2d& M) {
	double m11 = M(0, 0);
	double m12 = M(0, 1);
	double m21 = M(1, 0);
	double m22 = M(1, 1);
	double detM = m11*m22 - m12*m21;
	int sign = 1;
	if (detM < 0) {
		sign = -1;
	}
	else if (detM == 0) {
		sign = 0;
	}
	Matrix2d MM;
	MM << m22, -m21, -m12, m11;
	Matrix2d Q = M + sign*MM;
	double clen = Q.col(0).norm();
	Q = Q / clen;
	return Q;
}

void fillxMI(vector<T>& M_, vector<T>& MDK_, const Matrix3d& Mxx, const Matrix3d& Kxx, int index, const Vector2d& damping, double h)
{
	Matrix3d MDKxx = Mxx + damping(0) * h * Mxx + damping(1) * h * h * Kxx;
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			M_.push_back(T(index + j, index + k, Mxx(j, k)));
			MDK_.push_back(T(index + j, index + k, MDKxx(j, k)));
		}
	}
}

void fillxxMI(vector<T>& M_, vector<T>& MDK_, const Matrix3d& Mxx, const Matrix3d& Kxx, int i0, int i1, const Vector2d& damping, double h)
{
	Matrix3d MDKxx = Mxx + damping(0) * h * Mxx + damping(1) * h * h * Kxx;
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			//M_.push_back(T(i0 + j, i1 + k, Mxx(j, k)));
			//M_.push_back(T(i1 + k, i0 + j, Mxx(j, k)));
			MDK_.push_back(T(i0 + j, i1 + k, MDKxx(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, MDKxx(j, k)));
		}
	}
}

void fillXMI(vector<T>& M_, vector<T>& MDK_, const Matrix2d& MXX, const Matrix2d& KXX, int index, const Vector2d& damping, double h)
{
	Matrix2d MDKXX = MXX + damping(0) * h * MXX + damping(1) * h * h * KXX;
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
			M_.push_back(T(index + j, index + k, MXX(j, k)));
			MDK_.push_back(T(index + j, index + k, MDKXX(j, k)));
		}
	}
}

void fillXXMI(vector<T>& M_, vector<T>& MDK_, const Matrix2d& MXX, const Matrix2d& KXX, int i0, int i1, const Vector2d& damping, double h)
{
	Matrix2d MDKXX = MXX + damping(0) * h * MXX + damping(1) * h * h * KXX;
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
			//M_.push_back(T(i0 + j, i1 + k, MXX(j, k)));
			//M_.push_back(T(i1 + k, i0 + j, MXX(j, k)));
			MDK_.push_back(T(i0 + j, i1 + k, MDKXX(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, MDKXX(j, k)));
		}
	}
}

void fillXxMI(vector<T>& M_, vector<T>& MDK_, MatrixXd& MXx, MatrixXd& KXx, int i0, int i1, Vector2d& damping, double h)
{
	MatrixXd MDKXx = MXx + damping(0) * h * MXx + damping(1) * h * h * KXx;
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 3; k++) {
			//M_.push_back(T(i0 + j, i1 + k, MXx(j, k)));
			//M_.push_back(T(i1 + k, i0 + j, MXx(j, k)));
			MDK_.push_back(T(i0 + j, i1 + k, MDKXx(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, MDKXx(j, k)));
		}
	}
}

void fillEOLMembrane(const Face* face, VectorXd& fm, MatrixXd& Km)
{
	MatrixXd Fa, Fb, Fc;
	bool EOLA = face->v[0]->node->EoL; bool EOLB = face->v[1]->node->EoL; bool EOLC = face->v[2]->node->EoL;

	if (EOLA) Fa = deform_grad_v(face->v[0]);
	if (EOLB) Fb = deform_grad_v(face->v[1]);
	if (EOLC) Fc = deform_grad_v(face->v[2]);

	Vector3d fma = fm.segment<3>(0); Vector3d fmb = fm.segment<3>(3); Vector3d fmc = fm.segment<3>(6);

	Matrix3d Kmaa = Km.block<3, 3>(0, 0); Matrix3d Kmab = Km.block<3, 3>(0, 3); Matrix3d Kmac = Km.block<3, 3>(0, 6);
	Matrix3d Kmba = Km.block<3, 3>(3, 0); Matrix3d Kmbb = Km.block<3, 3>(3, 3); Matrix3d Kmbc = Km.block<3, 3>(3, 6);
	Matrix3d Kmca = Km.block<3, 3>(6, 0); Matrix3d Kmcb = Km.block<3, 3>(6, 3); Matrix3d Kmcc = Km.block<3, 3>(6, 6);

	int ja = 0; int jA = 3; int jb = 5; int jB = 8; int jc = 10; int jC = 13;

	fm.segment<3>(ja) = fma; fm.segment<3>(jb) = fmb; fm.segment<3>(jc) = fmc;

	Km.block<3, 3>(ja, ja) = Kmaa; Km.block<3, 3>(ja, jb) = Kmab; Km.block<3, 3>(ja, jc) = Kmac;
	/*Km.block<3, 3>(jb, ja) = Kmba;*/ Km.block<3, 3>(jb, jb) = Kmbb; Km.block<3, 3>(jb, jc) = Kmbc;
	/*Km.block<3, 3>(jc, ja) = Kmca; Km.block<3, 3>(jc, jb) = Kmcb;*/ Km.block<3, 3>(jc, jc) = Kmcc;

	if (EOLA) {
		fm.segment<2>(jA) = -Fa.transpose() * fma;
		Km.block<2, 2>(jA, jA) = Fa.transpose() * Kmaa * Fa;
	}
	if (EOLB) {
		fm.segment<2>(jB) = -Fb.transpose() * fmb;
		Km.block<2, 2>(jB, jB) = Fb.transpose() * Kmbb * Fb;
	}
	if (EOLC) {
		fm.segment<2>(jC) = -Fc.transpose() * fmc;
		Km.block<2, 2>(jC, jC) = Fc.transpose() * Kmcc * Fc;
	}

	if (EOLA && EOLB) {
		Km.block<2, 2>(jA, jB) = Fa.transpose() * Kmab * Fb;
	}
	if (EOLA && EOLC) {
		Km.block<2, 2>(jA, jC) = Fa.transpose() * Kmac * Fc;
	}
	if (EOLB && EOLC) {
		Km.block<2, 2>(jB, jC) = Fb.transpose() * Kmbc * Fc;
	}

	if (EOLA) {
		Km.block<2, 3>(jA, ja) = -Fa.transpose() * Kmaa;
		Km.block<2, 3>(jA, jb) = -Fa.transpose() * Kmab;
		Km.block<2, 3>(jA, jc) = -Fa.transpose() * Kmac;
	}
	if (EOLB) {
		Km.block<2, 3>(jB, ja) = -Fb.transpose() * Kmba;
		Km.block<2, 3>(jB, jb) = -Fb.transpose() * Kmbb;
		Km.block<2, 3>(jB, jc) = -Fb.transpose() * Kmbc;
	}
	if (EOLC) {
		Km.block<2, 3>(jC, ja) = -Fc.transpose() * Kmca;
		Km.block<2, 3>(jC, jb) = -Fc.transpose() * Kmcb;
		Km.block<2, 3>(jC, jc) = -Fc.transpose() * Kmcc;
	}
}

void faceBasedF(const Mesh& mesh, VectorXd& f, vector<T>& M_, vector<T>& MDK_, const Vector3d& grav, double h)
{
	for (int i = 0; i < mesh.faces.size(); i++) {
		Face* face = mesh.faces[i];

		double xa[3], xb[3], xc[3];
		double Xa[2], Xb[2], Xc[2];
		double g[3], PP[6], QQ[4];
		double Wi[1], Wm[1];

		const Material* mat = face->material;

		Vector3d txa = v2e(face->v[0]->node->x),
			txb = v2e(face->v[1]->node->x),
			txc = v2e(face->v[2]->node->x);
		Vector2d tXa = v322e(face->v[0]->u),
			tXb = v322e(face->v[1]->u),
			tXc = v322e(face->v[2]->u);

		Map<Vector3d>(xa, 3) = txa;
		Map<Vector3d>(xb, 3) = txb;
		Map<Vector3d>(xc, 3) = txc;
		Map<Vector2d>(Xa, 2) = tXa;
		Map<Vector2d>(Xb, 2) = tXb;
		Map<Vector2d>(Xc, 2) = tXc;

		Map<Vector3d>(g, grav.rows(), grav.cols()) = grav;

		MatrixXd Dxt(3, 2);
		MatrixXd DX(2, 2);
		Dxt << (txb - txa), (txc - txa);
		DX << (tXb - tXa), (tXc - tXa);
		Vector3d normm = (txb - txa).cross(txc - txa);
		Vector3d Pxm = (txb - txa) / (txb - txa).norm();
		Vector3d Pym = normm.cross(Pxm);
		Pym = Pym / Pym.norm();
		MatrixXd Pm(2, 3);
		Pm << Pxm.transpose(), Pym.transpose();
		MatrixXd Fm = Dxt * DX.inverse();
		Matrix2d Fbarm = Pm*Fm;
		Matrix2d Qm = poldec(Fbarm);

		Map<MatrixXd>(PP, Pm.rows(), Pm.cols()) = Pm;
		Map<Matrix2d>(QQ, Qm.rows(), Qm.cols()) = Qm;

		int aindex = face->v[0]->node->index * 3;
		int bindex = face->v[1]->node->index * 3;
		int cindex = face->v[2]->node->index * 3;
		int aindexX = mesh.nodes.size() * 3 + face->v[0]->node->EoL_index * 2;
		int bindexX = mesh.nodes.size() * 3 + face->v[1]->node->EoL_index * 2;
		int cindexX = mesh.nodes.size() * 3 + face->v[2]->node->EoL_index * 2;

		double fm[9], Km[81];

		VectorXd fme(15);
		MatrixXd Kme(15, 15);

		ComputeMembrane(xa, xb, xc, Xa, Xb, Xc, mat->e, mat->nu, PP, QQ, Wm, fm, Km);

		fme.segment<9>(0) = Map<VectorXd>(fm, 9);
		Kme.block<9, 9>(0, 0) = Map<MatrixXd>(Km, 9, 9);

		Vector2d damping(mat->dampingA, mat->dampingB);

		if (face->v[0]->node->EoL || face->v[1]->node->EoL || face->v[2]->node->EoL) {

			double fi[15], Mi[225];
			ComputeInertiaEOL(xa, xb, xc, Xa, Xb, Xc, g, mat->density, Wi, fi, Mi);
			VectorXd fie = Map<VectorXd>(fi, 15);
			MatrixXd Mie = Map<MatrixXd>(Mi, 15, 15);

			//MatrixXd F = deform_grad(face);
			//int ja = 0; int jA = 3; int jb = 5; int jB = 8; int jc = 10; int jC = 13;
			//Matrix3d Kmaa = Kme.block(0, 0, 3, 3); Matrix3d Kmab = Kme.block(0, 3, 3, 3); Matrix3d Kmac = Kme.block(0, 6, 3, 3);
			//Matrix3d Kmba = Kme.block(3, 0, 3, 3); Matrix3d Kmbb = Kme.block(3, 3, 3, 3); Matrix3d Kmbc = Kme.block(3, 6, 3, 3);
			//Matrix3d Kmca = Kme.block(6, 0, 3, 3); Matrix3d Kmcb = Kme.block(6, 3, 3, 3); Matrix3d Kmcc = Kme.block(6, 6, 3, 3);
			//Kme.block(ja, ja, 3, 3) = Kmaa; Kme.block(ja, jb, 3, 3) = Kmab; Kme.block(ja, jc, 3, 3) = Kmac;
			//Kme.block(jb, ja, 3, 3) = Kmba; Kme.block(jb, jb, 3, 3) = Kmbb; Kme.block(jb, jc, 3, 3) = Kmbc;
			//Kme.block(jc, ja, 3, 3) = Kmca; Kme.block(jc, jb, 3, 3) = Kmcb; Kme.block(jc, jc, 3, 3) = Kmcc;
			//Kme.block(ja, jA, 3, 2) = -Kmaa * F; Kme.block(ja, jB, 3, 2) = -Kmab * F; Kme.block(ja, jC, 3, 2) = -Kmac * F;
			//Kme.block(jb, jA, 3, 2) = -Kmba * F; Kme.block(jb, jB, 3, 2) = -Kmbb * F; Kme.block(jb, jC, 3, 2) = -Kmbc * F;
			//Kme.block(jc, jA, 3, 2) = -Kmca * F; Kme.block(jc, jB, 3, 2) = -Kmcb * F; Kme.block(jc, jC, 3, 2) = -Kmcc * F;
			//Kme.block(jA, ja, 2, 3) = -F.transpose() * Kmaa; Kme.block(jA, jb, 2, 3) = -F.transpose() * Kmab; Kme.block(jA, jc, 2, 3) = -F.transpose() * Kmac;
			//Kme.block(jB, ja, 2, 3) = -F.transpose() * Kmba; Kme.block(jB, jb, 2, 3) = -F.transpose() * Kmbb; Kme.block(jB, jc, 2, 3) = -F.transpose() * Kmbc;
			//Kme.block(jC, ja, 2, 3) = -F.transpose() * Kmca; Kme.block(jC, jb, 2, 3) = -F.transpose() * Kmcb; Kme.block(jC, jc, 2, 3) = -F.transpose() * Kmcc;
			//Kme.block(jA, jA, 2, 2) = F.transpose() * Kmaa * F; Kme.block(jA, jB, 2, 2) = F.transpose() * Kmab * F; Kme.block(jA, jC, 2, 2) = F.transpose() * Kmac * F;
			//Kme.block(jB, jA, 2, 2) = F.transpose() * Kmba * F; Kme.block(jB, jB, 2, 2) = F.transpose() * Kmbb * F; Kme.block(jB, jC, 2, 2) = F.transpose() * Kmbc * F;
			//Kme.block(jC, jA, 2, 2) = F.transpose() * Kmca * F; Kme.block(jC, jB, 2, 2) = F.transpose() * Kmcb * F; Kme.block(jC, jC, 2, 2) = F.transpose() * Kmcc * F;

			fillEOLMembrane(face, fme, Kme);

			f.segment<3>(aindex) += (fie.segment<3>(0) + fme.segment<3>(0));
			if (face->v[0]->node->EoL) f.segment<2>(aindexX) += (fie.segment<2>(3) + fme.segment<2>(3));

			f.segment<3>(bindex) += (fie.segment<3>(5) + fme.segment<3>(5));
			if (face->v[1]->node->EoL)f.segment<2>(bindexX) += (fie.segment<2>(8) + fme.segment<2>(8));

			f.segment<3>(cindex) += (fie.segment<3>(10) + fme.segment<3>(10));
			if (face->v[2]->node->EoL) f.segment<2>(cindexX) += (fie.segment<2>(13) + fme.segment<2>(13));

			//f.segment<3>(face->v[0]->node->index * 3) += (fie.segment<3>(0) + fme.segment<3>(0));
			//if (face->v[0]->node->EoL) f.segment<2>(aindexX) += -F.transpose() * (fie.segment<3>(0) + fme.segment<3>(0));

			//f.segment<3>(face->v[1]->node->index * 3) += (fie.segment<3>(5) + fme.segment<3>(3));
			//if (face->v[1]->node->EoL)f.segment<2>(bindexX) += -F.transpose() * (fie.segment<3>(5) + fme.segment<3>(3));

			//f.segment<3>(face->v[2]->node->index * 3) += (fie.segment<3>(10) + fme.segment<3>(6));
			//if (face->v[2]->node->EoL) f.segment<2>(cindexX) += -F.transpose() * (fie.segment<3>(10) + fme.segment<3>(6));

			cout << Mie << endl;

			// Diagonal x
			Matrix3d Mxx, Kxx;
			Mxx = Mie.block<3, 3>(0, 0); Kxx = Kme.block<3, 3>(0, 0);
			fillxMI(M_, MDK_, Mxx, Kxx, aindex, damping, h);
			Mxx = Mie.block<3, 3>(5, 5); Kxx = Kme.block<3, 3>(5, 5);
			fillxMI(M_, MDK_, Mxx, Kxx, bindex, damping, h);
			Mxx = Mie.block<3, 3>(10, 10); Kxx = Kme.block<3, 3>(10, 10);
			fillxMI(M_, MDK_, Mxx, Kxx, cindex, damping, h);

			// Off-Diagonal x
			Mxx = Mie.block<3, 3>(0, 5); Kxx = Kme.block<3, 3>(0, 5);
			fillxxMI(M_, MDK_, Mxx, Kxx, aindex, bindex, damping, h);
			Mxx = Mie.block<3, 3>(0, 10); Kxx = Kme.block<3, 3>(0, 10);
			fillxxMI(M_, MDK_, Mxx, Kxx, aindex, cindex, damping, h);
			Mxx = Mie.block<3, 3>(5, 10); Kxx = Kme.block<3, 3>(5, 10);
			fillxxMI(M_, MDK_, Mxx, Kxx, bindex, cindex, damping, h);

			// If has EOL componenent, Diagonal X
			Matrix2d MXX, KXX;
			if (face->v[0]->node->EoL) {
				MXX = Mie.block<2, 2>(3, 3); KXX = Kme.block<2, 2>(3, 3);
				fillXMI(M_, MDK_, MXX, KXX, aindexX, damping, h);
			}
			if (face->v[1]->node->EoL) {
				MXX = Mie.block<2, 2>(8, 8); KXX = Kme.block<2, 2>(8, 8);
				fillXMI(M_, MDK_, MXX, KXX, bindexX, damping, h);
			}
			if (face->v[2]->node->EoL) {
				MXX = Mie.block<2, 2>(13, 13); KXX = Kme.block<2, 2>(13, 13);
				fillXMI(M_, MDK_, MXX, KXX, cindexX, damping, h);
			}

			// If has EOL componenent, Off-Diagonal X
			if (face->v[0]->node->EoL && face->v[1]->node->EoL) {
				MXX = Mie.block<2, 2>(3, 8); KXX = Kme.block<2, 2>(3, 8);
				fillXXMI(M_, MDK_, MXX, KXX, aindexX, bindexX, damping, h);
			}
			if (face->v[0]->node->EoL && face->v[2]->node->EoL) {
				MXX = Mie.block<2, 2>(3, 13); KXX = Kme.block<2, 2>(3, 13);
				fillXXMI(M_, MDK_, MXX, KXX, aindexX, cindexX, damping, h);
			}
			if (face->v[1]->node->EoL && face->v[2]->node->EoL) {
				MXX = Mie.block<2, 2>(8, 13); KXX = Kme.block<2, 2>(8, 13);
				fillXXMI(M_, MDK_, MXX, KXX, bindexX, cindexX, damping, h);
			}

			// X-x values
			MatrixXd MXx, KXx;
			if (face->v[0]->node->EoL) {
				MXx = Mie.block<2, 3>(3, 0); KXx = Kme.block<2, 3>(3, 0);
				fillXxMI(M_, MDK_, MXx, KXx, aindexX, aindex, damping, h);
				MXx = Mie.block<2, 3>(3, 5); KXx = Kme.block<2, 3>(3, 5);
				fillXxMI(M_, MDK_, MXx, KXx, aindexX, bindex, damping, h);
				MXx = Mie.block<2, 3>(3, 10); KXx = Kme.block<2, 3>(3, 10);
				fillXxMI(M_, MDK_, MXx, KXx, aindexX, cindex, damping, h);
			}
			if (face->v[1]->node->EoL) {
				MXx = Mie.block<2, 3>(8, 0); KXx = Kme.block<2, 3>(8, 0);
				fillXxMI(M_, MDK_, MXx, KXx, bindexX, aindex, damping, h);
				MXx = Mie.block<2, 3>(8, 5); KXx = Kme.block<2, 3>(8, 5);
				fillXxMI(M_, MDK_, MXx, KXx, bindexX, bindex, damping, h);
				MXx = Mie.block<2, 3>(8, 10); KXx = Kme.block<2, 3>(8, 10);
				fillXxMI(M_, MDK_, MXx, KXx, bindexX, cindex, damping, h);
			}
			if (face->v[2]->node->EoL) {
				MXx = Mie.block<2, 3>(13, 0); KXx = Kme.block<2, 3>(13, 0);
				fillXxMI(M_, MDK_, MXx, KXx, cindexX, aindex, damping, h);
				MXx = Mie.block<2, 3>(13, 5); KXx = Kme.block<2, 3>(13, 5);
				fillXxMI(M_, MDK_, MXx, KXx, cindexX, bindex, damping, h);
				MXx = Mie.block<2, 3>(13, 10); KXx = Kme.block<2, 3>(13, 10);
				fillXxMI(M_, MDK_, MXx, KXx, cindexX, cindex, damping, h);
			}
		}
		else {

			double fi[9], Mi[81];
			ComputeInertial(xa, xb, xc, Xa, Xb, Xc, g, mat->density, Wi, fi, Mi);
			VectorXd fie = Map<VectorXd>(fi, 9);
			MatrixXd Mie = Map<MatrixXd>(Mi, 9, 9);

			f.segment<3>(aindex) += (fie.segment<3>(0) + fme.segment<3>(0));
			f.segment<3>(bindex) += (fie.segment<3>(3) + fme.segment<3>(3));
			f.segment<3>(cindex) += (fie.segment<3>(6) + fme.segment<3>(6));

			Matrix3d Mxx, Kxx;
			Mxx = Mie.block<3, 3>(0, 0); Kxx = Kme.block<3, 3>(0, 0);
			fillxMI(M_, MDK_, Mxx, Kxx, aindex, damping, h);
			Mxx = Mie.block<3, 3>(3, 3); Kxx = Kme.block<3, 3>(3, 3);
			fillxMI(M_, MDK_, Mxx, Kxx, bindex, damping, h);
			Mxx = Mie.block<3, 3>(6, 6); Kxx = Kme.block<3, 3>(6, 6);
			fillxMI(M_, MDK_, Mxx, Kxx, cindex, damping, h);

			Mxx = Mie.block<3, 3>(0, 3); Kxx = Kme.block<3, 3>(0, 3);
			fillxxMI(M_, MDK_, Mxx, Kxx, aindex, bindex, damping, h);
			Mxx = Mie.block<3, 3>(0, 6); Kxx = Kme.block<3, 3>(0, 6);
			fillxxMI(M_, MDK_, Mxx, Kxx, aindex, cindex, damping, h);
			Mxx = Mie.block<3, 3>(3, 6); Kxx = Kme.block<3, 3>(3, 6);
			fillxxMI(M_, MDK_, Mxx, Kxx, bindex, cindex, damping, h);
		}
	}
}

void fillxB(vector<T>& MDK_, const Matrix3d& Kxx, int index)
{
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			MDK_.push_back(T(index + j, index + k, Kxx(j, k)));
		}
	}
}

void fillxxB(vector<T>& MDK_, const Matrix3d& Kxx, int i0, int i1)
{
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			MDK_.push_back(T(i0 + j, i1 + k, Kxx(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, Kxx(j, k)));
		}
	}
}

void fillXB(vector<T>& MDK_, const Matrix2d& KXX, int index)
{
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
			MDK_.push_back(T(index + j, index + k, KXX(j, k)));
		}
	}
}

void fillXXB(vector<T>& MDK_, const Matrix2d& KXX, int i0, int i1)
{
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
			MDK_.push_back(T(i0 + j, i1 + k, KXX(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, KXX(j, k)));
		}
	}
}

void fillXxB(vector<T>& MDK_, const MatrixXd& KXx, int i0, int i1)
{
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 3; k++) {
			MDK_.push_back(T(i0 + j, i1 + k, KXx(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, KXx(j, k)));
		}
	}
}

void fillEOLMembrane(const Vert* v0, const Vert* v1, const Vert* v2, const Vert* v3, VectorXd& fb, MatrixXd& Kb)
{
	MatrixXd Fa, Fb, Fc, Fd;
	bool EOLA = v0->node->EoL; bool EOLB = v1->node->EoL; bool EOLC = v2->node->EoL; bool EOLD = v3->node->EoL;

	if (EOLA) Fa = deform_grad_v(v0);
	if (EOLB) Fb = deform_grad_v(v1);
	if (EOLC) Fc = deform_grad_v(v2);
	if (EOLD) Fd = deform_grad_v(v3);

	Vector3d fba = fb.segment<3>(0); Vector3d fbb = fb.segment<3>(3); Vector3d fbc = fb.segment<3>(6); Vector3d fbd = fb.segment<3>(9);

	Matrix3d Kbaa = Kb.block<3, 3>(0, 0); Matrix3d Kbab = Kb.block<3, 3>(0, 3); Matrix3d Kbac = Kb.block<3, 3>(0, 6); MatrixXd Kbad = Kb.block<3, 3>(0, 9);
	Matrix3d Kbba = Kb.block<3, 3>(3, 0); Matrix3d Kbbb = Kb.block<3, 3>(3, 3); Matrix3d Kbbc = Kb.block<3, 3>(3, 6); MatrixXd Kbbd = Kb.block<3, 3>(3, 9);
	Matrix3d Kbca = Kb.block<3, 3>(6, 0); Matrix3d Kbcb = Kb.block<3, 3>(6, 3); Matrix3d Kbcc = Kb.block<3, 3>(6, 6); MatrixXd Kbcd = Kb.block<3, 3>(6, 9);
	Matrix3d Kbda = Kb.block<3, 3>(9, 0); Matrix3d Kbdb = Kb.block<3, 3>(9, 3); Matrix3d Kbdc = Kb.block<3, 3>(9, 6); MatrixXd Kbdd = Kb.block<3, 3>(9, 9);

	int ja = 0; int jA = 3; int jb = 5; int jB = 8; int jc = 10; int jC = 13; int jd = 15; int jD = 18;

	fb.segment<3>(ja) = fba; fb.segment<3>(jb) = fbb; fb.segment<3>(jc) = fbc; fb.segment<3>(jd) = fbd;

	Kb.block<3, 3>(ja, ja) = Kbaa; Kb.block<3, 3>(ja, jb) = Kbab; Kb.block<3, 3>(ja, jc) = Kbac; Kb.block<3, 3>(ja, jd) = Kbad;
	/*Kd.block<3, 3>(jb, ja) = Kbba;*/ Kb.block<3, 3>(jb, jb) = Kbbb; Kb.block<3, 3>(jb, jc) = Kbbc; Kb.block<3, 3>(jb, jd) = Kbbd;
	/*Kd.block<3, 3>(jc, ja) = Kbca; Kb.block<3, 3>(jc, jb) = Kbcb;*/ Kb.block<3, 3>(jc, jc) = Kbcc; Kb.block<3, 3>(jc, jd) = Kbcd;
	/*Kd.block<3, 3>(jd, ja) = Kbda; Kb.block<3, 3>(jd, jb) = Kbdb; Kb.block<3, 3>(jd, jc) = Kbdc;*/ Kb.block<3, 3>(jd, jd) = Kbdd;

	if (EOLA) {
		fb.segment<2>(jA) = -Fa.transpose() * fba;
		Kb.block<2, 2>(jA, jA) = Fa.transpose() * Kbaa * Fa;
	}
	if (EOLB) {
		fb.segment<2>(jB) = -Fb.transpose() * fbb;
		Kb.block<2, 2>(jB, jB) = Fb.transpose() * Kbbb * Fb;
	}
	if (EOLC) {
		fb.segment<2>(jC) = -Fc.transpose() * fbc;
		Kb.block<2, 2>(jC, jC) = Fc.transpose() * Kbcc * Fc;
	}
	if (EOLD) {
		fb.segment<2>(jD) = -Fd.transpose() * fbd;
		Kb.block<2, 2>(jD, jD) = Fd.transpose() * Kbdd * Fd;
	}

	if (EOLA && EOLB) {
		Kb.block<2, 2>(jA, jB) = Fa.transpose() * Kbab * Fb;
	}
	if (EOLA && EOLC) {
		Kb.block<2, 2>(jA, jC) = Fa.transpose() * Kbac * Fc;
	}
	if (EOLA && EOLD) {
		Kb.block<2, 2>(jA, jD) = Fa.transpose() * Kbad * Fd;
	}
	if (EOLB && EOLC) {
		Kb.block<2, 2>(jB, jC) = Fb.transpose() * Kbbc * Fc;
	}
	if (EOLB && EOLD) {
		Kb.block<2, 2>(jB, jD) = Fb.transpose() * Kbbd * Fd;
	}
	if (EOLC && EOLD) {
		Kb.block<2, 2>(jC, jD) = Fc.transpose() * Kbcd * Fd;
	}

	if (EOLA) {
		Kb.block<2, 3>(jA, ja) = -Fa.transpose() * Kbaa;
		Kb.block<2, 3>(jA, jb) = -Fa.transpose() * Kbab;
		Kb.block<2, 3>(jA, jc) = -Fa.transpose() * Kbac;
		Kb.block<2, 3>(jA, jd) = -Fa.transpose() * Kbad;
	}
	if (EOLB) {
		Kb.block<2, 3>(jB, ja) = -Fb.transpose() * Kbba;
		Kb.block<2, 3>(jB, jb) = -Fb.transpose() * Kbbb;
		Kb.block<2, 3>(jB, jc) = -Fb.transpose() * Kbbc;
		Kb.block<2, 3>(jB, jd) = -Fb.transpose() * Kbbd;
	}
	if (EOLC) {
		Kb.block<2, 3>(jC, ja) = -Fc.transpose() * Kbca;
		Kb.block<2, 3>(jC, jb) = -Fc.transpose() * Kbcb;
		Kb.block<2, 3>(jC, jc) = -Fc.transpose() * Kbcc;
		Kb.block<2, 3>(jC, jd) = -Fc.transpose() * Kbcd;
	}
	if (EOLD) {
		Kb.block<2, 3>(jD, ja) = -Fd.transpose() * Kbda;
		Kb.block<2, 3>(jD, jb) = -Fd.transpose() * Kbdb;
		Kb.block<2, 3>(jD, jc) = -Fd.transpose() * Kbdc;
		Kb.block<2, 3>(jD, jd) = -Fd.transpose() * Kbdd;
	}
}

void edgeBasedF(const Mesh& mesh, const Material& mat, VectorXd& f, vector<T>& M_, vector<T>& MDK_, double h)
{
	for (int e = 0; e < mesh.edges.size(); e++) {
		if (mesh.edges[e]->adjf[0] == NULL || mesh.edges[e]->adjf[1] == NULL) {
			continue;
		}
		Edge* edge = mesh.edges[e];
		Face* f0 = edge->adjf[0], *f1 = edge->adjf[1];
		Node* n0 = edge->n[0], *n1 = edge->n[1];
		Vert* v0 = n0->verts[0], *v1 = n1->verts[0];

		Vert* v2 = get_other_vert(f0, v0, v1), *v3 = get_other_vert(f1, v0, v1);
		Node* n2 = v2->node, *n3 = v3->node;

		double xa[3], xb[3], xc[3], xd[3];
		double Xa[2], Xb[2], Xc[2], Xd[2];

		Vector3d txa = v2e(n0->x),
			txb = v2e(n1->x),
			txc = v2e(n2->x),
			txd = v2e(n3->x);
		Vector2d tXa = v322e(v0->u),
			tXb = v322e(v1->u),
			tXc = v322e(v2->u),
			tXd = v322e(v3->u);

		Map<Vector3d>(xa, 3) = txa;
		Map<Vector3d>(xb, 3) = txb;
		Map<Vector3d>(xc, 3) = txc;
		Map<Vector3d>(xd, 3) = txd;
		Map<Vector2d>(Xa, 2) = tXa;
		Map<Vector2d>(Xb, 2) = tXb;
		Map<Vector2d>(Xc, 2) = tXc;
		Map<Vector2d>(Xd, 2) = tXd;

		bool to_eolA = n0->EoL;
		bool to_eolB = n1->EoL;
		bool to_eolC = n2->EoL;
		bool to_eolD = n3->EoL;

		int aindex = n0->index * 3;
		int bindex = n1->index * 3;
		int cindex = n2->index * 3;
		int dindex = n3->index * 3;
		int aindexX = mesh.nodes.size() * 3 + n0->EoL_index * 2;
		int bindexX = mesh.nodes.size() * 3 + n1->EoL_index * 2;
		int cindexX = mesh.nodes.size() * 3 + n2->EoL_index * 2;
		int dindexX = mesh.nodes.size() * 3 + n3->EoL_index * 2;

		double Wb[1], fb[12], Kb[144];

		Vector2d damping(mat.dampingA, mat.dampingB);

		VectorXd fbe(20);
		MatrixXd Kbe(20, 20);

		ComputeBending(xa, xb, xc, xd, Xa, Xb, Xc, Xd, mat.beta, Wb, fb, Kb);

		fbe.segment<12>(0) = Map<VectorXd>(fb, 12);
		Kbe.block<12, 12>(0, 0) = Map<MatrixXd>(Kb, 12, 12);

		if (to_eolA || to_eolB || to_eolC || to_eolD) {
			MatrixXd F1 = deform_grad(f0);
			MatrixXd F2 = deform_grad(f1);
			MatrixXd F = (F1 + F2) / 2;

			int ja = 0; int jA = 3; int jb = 5; int jB = 8; int jc = 10; int jC = 13; int jd = 15; int jD = 18;
			Matrix3d Kbaa = Kbe.block(0, 0, 3, 3); Matrix3d Kbab = Kbe.block(0, 3, 3, 3); Matrix3d Kbac = Kbe.block(0, 6, 3, 3); Matrix3d Kbad = Kbe.block(0, 9, 3, 3);
			Matrix3d Kbba = Kbe.block(3, 0, 3, 3); Matrix3d Kbbb = Kbe.block(3, 3, 3, 3); Matrix3d Kbbc = Kbe.block(3, 6, 3, 3); Matrix3d Kbbd = Kbe.block(3, 9, 3, 3);
			Matrix3d Kbca = Kbe.block(6, 0, 3, 3); Matrix3d Kbcb = Kbe.block(6, 3, 3, 3); Matrix3d Kbcc = Kbe.block(6, 6, 3, 3); Matrix3d Kbcd = Kbe.block(6, 9, 3, 3);
			Matrix3d Kbda = Kbe.block(9, 0, 3, 3); Matrix3d Kbdb = Kbe.block(9, 3, 3, 3); Matrix3d Kbdc = Kbe.block(9, 6, 3, 3); Matrix3d Kbdd = Kbe.block(9, 9, 3, 3);
			Kbe.block(ja, ja, 3, 3) = Kbaa; Kbe.block(ja, jb, 3, 3) = Kbab; Kbe.block(ja, jc, 3, 3) = Kbac; Kbe.block(ja, jd, 3, 3) = Kbad;
			Kbe.block(jb, ja, 3, 3) = Kbba; Kbe.block(jb, jb, 3, 3) = Kbbb; Kbe.block(jb, jc, 3, 3) = Kbbc; Kbe.block(jb, jd, 3, 3) = Kbbd;
			Kbe.block(jc, ja, 3, 3) = Kbca; Kbe.block(jc, jb, 3, 3) = Kbcb; Kbe.block(jc, jc, 3, 3) = Kbcc; Kbe.block(jc, jd, 3, 3) = Kbcd;
			Kbe.block(jd, ja, 3, 3) = Kbda; Kbe.block(jd, jb, 3, 3) = Kbdb; Kbe.block(jd, jc, 3, 3) = Kbdc; Kbe.block(jd, jd, 3, 3) = Kbdd;
			Kbe.block(ja, jA, 3, 2) = -Kbaa * F; Kbe.block(ja, jB, 3, 2) = -Kbab * F; Kbe.block(ja, jC, 3, 2) = -Kbac * F; Kbe.block(ja, jD, 3, 2) = -Kbad * F;
			Kbe.block(jb, jA, 3, 2) = -Kbba * F; Kbe.block(jb, jB, 3, 2) = -Kbbb * F; Kbe.block(jb, jC, 3, 2) = -Kbbc * F; Kbe.block(jb, jD, 3, 2) = -Kbbd * F;
			Kbe.block(jc, jA, 3, 2) = -Kbca * F; Kbe.block(jc, jB, 3, 2) = -Kbcb * F; Kbe.block(jc, jC, 3, 2) = -Kbcc * F; Kbe.block(jc, jD, 3, 2) = -Kbcd * F;
			Kbe.block(jd, jA, 3, 2) = -Kbda * F; Kbe.block(jd, jB, 3, 2) = -Kbdb * F; Kbe.block(jd, jC, 3, 2) = -Kbdc * F; Kbe.block(jd, jD, 3, 2) = -Kbdd * F;
			Kbe.block(jA, ja, 2, 3) = -F.transpose() * Kbaa; Kbe.block(jA, jb, 2, 3) = -F.transpose() * Kbab; Kbe.block(jA, jc, 2, 3) = -F.transpose() * Kbac; Kbe.block(jA, jd, 2, 3) = -F.transpose() * Kbad;
			Kbe.block(jB, ja, 2, 3) = -F.transpose() * Kbba; Kbe.block(jB, jb, 2, 3) = -F.transpose() * Kbbb; Kbe.block(jB, jc, 2, 3) = -F.transpose() * Kbbc; Kbe.block(jB, jd, 2, 3) = -F.transpose() * Kbbd;
			Kbe.block(jC, ja, 2, 3) = -F.transpose() * Kbca; Kbe.block(jC, jb, 2, 3) = -F.transpose() * Kbcb; Kbe.block(jC, jc, 2, 3) = -F.transpose() * Kbcc; Kbe.block(jC, jd, 2, 3) = -F.transpose() * Kbcd;
			Kbe.block(jD, ja, 2, 3) = -F.transpose() * Kbda; Kbe.block(jD, jb, 2, 3) = -F.transpose() * Kbdb; Kbe.block(jD, jc, 2, 3) = -F.transpose() * Kbdc; Kbe.block(jD, jd, 2, 3) = -F.transpose() * Kbdd;
			Kbe.block(jA, jA, 2, 2) = F.transpose() * Kbaa * F; Kbe.block(jA, jB, 2, 2) = F.transpose() * Kbab * F; Kbe.block(jA, jC, 2, 2) = F.transpose() * Kbac * F; Kbe.block(jA, jD, 2, 2) = F.transpose() * Kbad * F;
			Kbe.block(jB, jA, 2, 2) = F.transpose() * Kbba * F; Kbe.block(jB, jB, 2, 2) = F.transpose() * Kbbb * F; Kbe.block(jB, jC, 2, 2) = F.transpose() * Kbbc * F; Kbe.block(jB, jD, 2, 2) = F.transpose() * Kbbd * F;
			Kbe.block(jC, jA, 2, 2) = F.transpose() * Kbca * F; Kbe.block(jC, jB, 2, 2) = F.transpose() * Kbcb * F; Kbe.block(jC, jC, 2, 2) = F.transpose() * Kbcc * F; Kbe.block(jC, jD, 2, 2) = F.transpose() * Kbcd * F;
			Kbe.block(jD, jA, 2, 2) = F.transpose() * Kbda * F; Kbe.block(jD, jB, 2, 2) = F.transpose() * Kbdb * F; Kbe.block(jD, jC, 2, 2) = F.transpose() * Kbdc * F; Kbe.block(jD, jD, 2, 2) = F.transpose() * Kbdd * F;

			f.segment<3>(aindex) += fbe.segment<3>(0);
			if (to_eolA) f.segment<2>(aindexX) += -F.transpose() * fbe.segment<3>(0);

			f.segment<3>(bindex) += fbe.segment<3>(3);
			if (to_eolB) f.segment<2>(bindexX) += -F.transpose() * fbe.segment<3>(3);

			f.segment<3>(cindex) += fbe.segment<3>(6);
			if (to_eolC) f.segment<2>(cindexX) += -F.transpose() * fbe.segment<3>(6);

			f.segment<3>(dindex) += fbe.segment<3>(9);
			if (to_eolD) f.segment<2>(dindexX) += -F.transpose() * fbe.segment<3>(9);

			Matrix3d Kxx;
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(0, 0);
			fillxB(MDK_, Kxx, aindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(5, 5);
			fillxB(MDK_, Kxx, bindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(10, 10);
			fillxB(MDK_, Kxx, cindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(15, 15);
			fillxB(MDK_, Kxx, dindex);

			Kxx = damping(1) * h * h * Kbe.block<3, 3>(0, 5);
			fillxxB(MDK_, Kxx, aindex, bindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(0, 10);
			fillxxB(MDK_, Kxx, aindex, cindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(0, 15);
			fillxxB(MDK_, Kxx, aindex, dindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(5, 10);
			fillxxB(MDK_, Kxx, bindex, cindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(5, 15);
			fillxxB(MDK_, Kxx, bindex, dindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(10, 15);
			fillxxB(MDK_, Kxx, cindex, dindex);

			Matrix2d KXX;
			if (to_eolA) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(3, 3);
				fillXB(MDK_, KXX, aindexX);
			}
			if (to_eolB) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(8, 8);
				fillXB(MDK_, KXX, bindexX);
			}
			if (to_eolC) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(13, 13);
				fillXB(MDK_, KXX, cindexX);
			}
			if (to_eolD) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(18, 18);
				fillXB(MDK_, KXX, dindexX);
			}

			if (to_eolA && to_eolB) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(3, 8);
				fillXXB(MDK_, KXX, aindexX, bindexX);
			}
			if (to_eolA && to_eolC) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(3, 13);
				fillXXB(MDK_, KXX, aindexX, cindexX);
			}
			if (to_eolA && to_eolD) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(3, 18);
				fillXXB(MDK_, KXX, aindexX, dindexX);
			}
			if (to_eolB && to_eolC) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(8, 13);
				fillXXB(MDK_, KXX, bindexX, cindexX);
			}
			if (to_eolB && to_eolD) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(8, 18);
				fillXXB(MDK_, KXX, bindexX, dindexX);
			}
			if (to_eolC && to_eolD) {
				KXX = damping(1) * h * h * Kbe.block<2, 2>(13, 18);
				fillXXB(MDK_, KXX, cindexX, dindexX);
			}

			MatrixXd KXx;
			if (to_eolA) {
				KXx = damping(1) * h * h * Kbe.block<2, 3>(3, 0);
				fillXxB(MDK_, KXx, aindexX, aindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(3, 5);
				fillXxB(MDK_, KXx, aindexX, bindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(3, 10);
				fillXxB(MDK_, KXx, aindexX, cindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(3, 15);
				fillXxB(MDK_, KXx, aindexX, dindex);
			}
			if (to_eolB) {
				KXx = damping(1) * h * h * Kbe.block<2, 3>(8, 0);
				fillXxB(MDK_, KXx, bindexX, aindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(8, 5);
				fillXxB(MDK_, KXx, bindexX, bindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(8, 10);
				fillXxB(MDK_, KXx, bindexX, cindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(8, 15);
				fillXxB(MDK_, KXx, bindexX, dindex);
			}
			if (to_eolC) {
				KXx = damping(1) * h * h * Kbe.block<2, 3>(13, 0);
				fillXxB(MDK_, KXx, cindexX, aindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(13, 5);
				fillXxB(MDK_, KXx, cindexX, bindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(13, 10);
				fillXxB(MDK_, KXx, cindexX, cindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(13, 15);
				fillXxB(MDK_, KXx, cindexX, dindex);
			}
			if (to_eolD) {
				KXx = damping(1) * h * h * Kbe.block<2, 3>(18, 0);
				fillXxB(MDK_, KXx, dindexX, aindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(18, 5);
				fillXxB(MDK_, KXx, dindexX, bindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(18, 10);
				fillXxB(MDK_, KXx, dindexX, cindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(18, 15);
				fillXxB(MDK_, KXx, dindexX, dindex);
			}
		}
		else {
			Matrix3d Kxx;
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(0, 0);
			fillxB(MDK_, Kxx, aindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(3, 3);
			fillxB(MDK_, Kxx, bindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(6, 6);
			fillxB(MDK_, Kxx, cindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(9, 9);
			fillxB(MDK_, Kxx, dindex);

			Kxx = damping(1) * h * h * Kbe.block<3, 3>(0, 3);
			fillxxB(MDK_, Kxx, aindex, bindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(0, 6);
			fillxxB(MDK_, Kxx, aindex, cindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(0, 9);
			fillxxB(MDK_, Kxx, aindex, dindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(3, 6);
			fillxxB(MDK_, Kxx, bindex, cindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(3, 9);
			fillxxB(MDK_, Kxx, bindex, dindex);
			Kxx = damping(1) * h * h * Kbe.block<3, 3>(6, 9);
			fillxxB(MDK_, Kxx, cindex, dindex);
		}
	}
}

void Forces::fill(const Mesh& mesh, const Material& mat, const Vector3d& grav, double h)
{
	f.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2);
	f.setZero();
	vector<T> M_;
	vector<T> MDK_;

	faceBasedF(mesh, f, M_, MDK_, grav, h);
	edgeBasedF(mesh, mat, f, M_, MDK_, h);

	M.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2, mesh.nodes.size() * 3 + mesh.EoL_Count * 2);
	MDK.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2, mesh.nodes.size() * 3 + mesh.EoL_Count * 2);

	M.setFromTriplets(M_.begin(), M_.end());
	MDK.setFromTriplets(MDK_.begin(), MDK_.end());
}