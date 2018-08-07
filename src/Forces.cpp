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

// We are using mass lumping so the mass matrix is simplified

void nodeBasedF(const Mesh& mesh, VectorXd& f, vector<T>& M_, const Vector3d& grav)
{
	for (int n = 0; n < mesh.nodes.size(); n++) {
		Node* node = mesh.nodes[n];
		Vert* vert = node->verts[0];

		int nindex = node->index * 3;

		double mass = 0.0;
		for (int f = 0; f < vert->adjf.size(); f++) {
			mass += (vert->adjf[f]->m * (1.0 / 3.0));
		}

		Vector3d mg = mass * grav;
		f.segment<3>(nindex) += mg;

		Matrix3d Mxx = Matrix3d::Identity() * mass;
		for (int j = 0; j < 3; j++) {
			M_.push_back(T(nindex + j, nindex + j, Mxx(j, j))); // We know it's a diagonal matrix
		}

		if (node->EoL) {
			int nindexX = mesh.nodes.size() * 3 + node->EoL_index * 2;

			MatrixXd F = deform_grad_v(vert);

			Vector2d mFg = -mass * F.transpose() * grav;
			f.segment<2>(nindexX) += mFg;

			Matrix2d MXX = -mass * F.transpose() * F;
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++) {
					M_.push_back(T(nindexX + j, nindexX + k, MXX(j, k)));
				}
			}

			MatrixXd MXx = -mass * F.transpose();
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(nindexX + j, nindex + k, MXx(j, k)));
					M_.push_back(T(nindex + k, nindexX + j, MXx(j, k)));
				}
			}
		}
	}
}

void fillxMI(vector<T>& MDK_, vector<T>& M_, const Matrix3d& Kxx, const Matrix3d& Mxx, int index, const Vector2d& damping, double h)
{
	Matrix3d MDKxx = Mxx + damping(1) * h * h * Kxx;
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			M_.push_back(T(index + j, index + k, Mxx(j, k)));
			MDK_.push_back(T(index + j, index + k, MDKxx(j, k)));
		}
	}
}

void fillxxMI(vector<T>& MDK_, vector<T>& M_, const Matrix3d& Kxx, const Matrix3d& Mxx, int i0, int i1, const Vector2d& damping, double h)
{
	Matrix3d MDKxx = Mxx + damping(1) * h * h * Kxx;
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			M_.push_back(T(i0 + j, i1 + k, Mxx(j, k)));
			M_.push_back(T(i1 + k, i0 + j, Mxx(j, k)));
			MDK_.push_back(T(i0 + j, i1 + k, MDKxx(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, MDKxx(j, k)));
		}
	}
}

void fillXMI(vector<T>& MDK_, vector<T>& M_, const Matrix2d& KXX, const Matrix2d& Mxx, int index, const Vector2d& damping, double h)
{
	Matrix2d MDKXX = Mxx + damping(1) * h * h * KXX;
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
			M_.push_back(T(index + j, index + k, Mxx(j, k)));
			MDK_.push_back(T(index + j, index + k, MDKXX(j, k)));
		}
	}
}

void fillXXMI(vector<T>& MDK_, vector<T>& M_, const Matrix2d& KXX, const Matrix2d& Mxx, int i0, int i1, const Vector2d& damping, double h)
{
	Matrix2d MDKXX = Mxx + damping(1) * h * h * KXX;
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
			M_.push_back(T(i0 + j, i1 + k, Mxx(j, k)));
			M_.push_back(T(i1 + k, i0 + j, Mxx(j, k)));
			MDK_.push_back(T(i0 + j, i1 + k, MDKXX(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, MDKXX(j, k)));
		}
	}
}

void fillXxMI(vector<T>& MDK_, vector<T>& M_, MatrixXd& KXx, const MatrixXd& Mxx, int i0, int i1, Vector2d& damping, double h)
{
	MatrixXd MDKXx = Mxx + damping(1) * h * h * KXx;
	for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 3; k++) {
			M_.push_back(T(i0 + j, i1 + k, Mxx(j, k)));
			M_.push_back(T(i1 + k, i0 + j, Mxx(j, k)));
			MDK_.push_back(T(i0 + j, i1 + k, MDKXx(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, MDKXx(j, k)));
		}
	}
}

void fillxXMI(vector<T>& MDK_, vector<T>& M_, MatrixXd& KxX, const MatrixXd& Mxx, int i0, int i1, Vector2d& damping, double h)
{
	MatrixXd MDKxX = Mxx + damping(1) * h * h * KxX;
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 2; k++) {
			M_.push_back(T(i0 + j, i1 + k, Mxx(j, k)));
			M_.push_back(T(i1 + k, i0 + j, Mxx(j, k)));
			MDK_.push_back(T(i0 + j, i1 + k, MDKxX(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, MDKxX(j, k)));
		}
	}
}

void fillEOLInertia(const Face* face, VectorXd& fi, MatrixXd& Mi)
{
	MatrixXd Fa, Fb, Fc;
	bool EOLA = face->v[0]->node->EoL; bool EOLB = face->v[1]->node->EoL; bool EOLC = face->v[2]->node->EoL;

	if (EOLA) Fa = deform_grad_v(face->v[0]);
	if (EOLB) Fb = deform_grad_v(face->v[1]);
	if (EOLC) Fc = deform_grad_v(face->v[2]);

	Vector3d fia = fi.segment<3>(0); Vector3d fib = fi.segment<3>(3); Vector3d fic = fi.segment<3>(6);

	Matrix3d Kiaa = Mi.block<3, 3>(0, 0); Matrix3d Kiab = Mi.block<3, 3>(0, 3); Matrix3d Kiac = Mi.block<3, 3>(0, 6);
	Matrix3d Kiba = Mi.block<3, 3>(3, 0); Matrix3d Kibb = Mi.block<3, 3>(3, 3); Matrix3d Kibc = Mi.block<3, 3>(3, 6);
	Matrix3d Kica = Mi.block<3, 3>(6, 0); Matrix3d Kicb = Mi.block<3, 3>(6, 3); Matrix3d Kicc = Mi.block<3, 3>(6, 6);

	int ja = 0; int jA = 3; int jb = 5; int jB = 8; int jc = 10; int jC = 13;

	fi.segment<3>(ja) = fia; fi.segment<3>(jb) = fib; fi.segment<3>(jc) = fic;

	Mi.block<3, 3>(ja, ja) = Kiaa; Mi.block<3, 3>(ja, jb) = Kiab; Mi.block<3, 3>(ja, jc) = Kiac;
	/*Mi.block<3, 3>(jb, ja) = Kiba;*/ Mi.block<3, 3>(jb, jb) = Kibb; Mi.block<3, 3>(jb, jc) = Kibc;
	/*Mi.block<3, 3>(jc, ja) = Kica; Mi.block<3, 3>(jc, jb) = Kicb;*/ Mi.block<3, 3>(jc, jc) = Kicc;

	if (EOLA) {
		fi.segment<2>(jA) = -Fa.transpose() * fia;
		Mi.block<2, 2>(jA, jA) = Fa.transpose() * Kiaa * Fa;
	}
	if (EOLB) {
		fi.segment<2>(jB) = -Fb.transpose() * fib;
		Mi.block<2, 2>(jB, jB) = Fb.transpose() * Kibb * Fb;
	}
	if (EOLC) {
		fi.segment<2>(jC) = -Fc.transpose() * fic;
		Mi.block<2, 2>(jC, jC) = Fc.transpose() * Kicc * Fc;
	}

	if (EOLA && EOLB) {
		Mi.block<2, 2>(jA, jB) = Fa.transpose() * Kiab * Fb;
	}
	if (EOLA && EOLC) {
		Mi.block<2, 2>(jA, jC) = Fa.transpose() * Kiac * Fc;
	}
	if (EOLB && EOLC) {
		Mi.block<2, 2>(jB, jC) = Fb.transpose() * Kibc * Fc;
	}

	if (EOLA) {
		Mi.block<2, 3>(jA, ja) = -Fa.transpose() * Kiaa;
		Mi.block<2, 3>(jA, jb) = -Fa.transpose() * Kiab;
		Mi.block<2, 3>(jA, jc) = -Fa.transpose() * Kiac;

		Mi.block<3, 2>(ja, jA) = -Kiaa * Fa;
	}
	if (EOLB) {
		//Mi.block<2, 3>(jB, ja) = -Fb.transpose() * Kiba;
		Mi.block<2, 3>(jB, jb) = -Fb.transpose() * Kibb;
		Mi.block<2, 3>(jB, jc) = -Fb.transpose() * Kibc;

		Mi.block<3, 2>(ja, jB) = -Kiab * Fb;
		Mi.block<3, 2>(jb, jB) = -Kibb * Fb;
	}
	if (EOLC) {
		//Mi.block<2, 3>(jC, ja) = -Fc.transpose() * Kica;
		//Mi.block<2, 3>(jC, jb) = -Fc.transpose() * Kicb;
		Mi.block<2, 3>(jC, jc) = -Fc.transpose() * Kicc;

		Mi.block<3, 2>(ja, jC) = -Kiac * Fc;
		Mi.block<3, 2>(jb, jC) = -Kibc * Fc;
		Mi.block<3, 2>(jc, jC) = -Kicc * Fc;
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

		Km.block<3, 2>(ja, jA) = -Kmaa * Fa;
	}
	if (EOLB) {
		//Km.block<2, 3>(jB, ja) = -Fb.transpose() * Kmba;
		Km.block<2, 3>(jB, jb) = -Fb.transpose() * Kmbb;
		Km.block<2, 3>(jB, jc) = -Fb.transpose() * Kmbc;

		Km.block<3, 2>(ja, jB) = -Kmab * Fb;
		Km.block<3, 2>(jb, jB) = -Kmbb * Fb;
	}
	if (EOLC) {
		//Km.block<2, 3>(jC, ja) = -Fc.transpose() * Kmca;
		//Km.block<2, 3>(jC, jb) = -Fc.transpose() * Kmcb;
		Km.block<2, 3>(jC, jc) = -Fc.transpose() * Kmcc;

		Km.block<3, 2>(ja, jC) = -Kmac * Fc;
		Km.block<3, 2>(jb, jC) = -Kmbc * Fc;
		Km.block<3, 2>(jc, jC) = -Kmcc * Fc;
	}
}

void faceBasedF(const Mesh& mesh, VectorXd& f, vector<T>& MDK_, vector<T>& M_, double h)
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
		double fi[9], Mi[81];
		double gg[3] = {0.0, 0.0, -9.81};

		VectorXd fme(15), fie(15);
		MatrixXd Kme(15, 15), Mie(15, 15);

		ComputeMembrane(xa, xb, xc, Xa, Xb, Xc, mat->e, mat->nu, PP, QQ, Wm, fm, Km);
		//ComputeInertiaEOL(xa, xb, xc, Xa, Xb, Xc, g, material.density, Wi, fi, Mi);
		ComputeInertial(xa, xb, xc, Xa, Xb, Xc, gg, mat->density, Wi, fi, Mi);

		fme.segment<9>(0) = Map<VectorXd>(fm, 9);
		Kme.block<9, 9>(0, 0) = Map<MatrixXd>(Km, 9, 9);
		fie.segment<9>(0) = Map<VectorXd>(fi, 9);
		Mie.block<9, 9>(0, 0) = Map<MatrixXd>(Mi, 9, 9);

		Vector2d damping(mat->dampingA, mat->dampingB);

		if (face->v[0]->node->EoL || face->v[1]->node->EoL || face->v[2]->node->EoL) {

			fillEOLInertia(face, fie, Mie);

			fillEOLMembrane(face, fme, Kme);

			f.segment<3>(aindex) += (fme.segment<3>(0) + fie.segment<3>(0));
			if (face->v[0]->node->EoL) f.segment<2>(aindexX) += (fme.segment<2>(3) + fie.segment<2>(3));

			f.segment<3>(bindex) += (fme.segment<3>(5) + fie.segment<3>(5));
			if (face->v[1]->node->EoL)f.segment<2>(bindexX) += (fme.segment<2>(8) + fie.segment<2>(8));

			f.segment<3>(cindex) += (fme.segment<3>(10) + fie.segment<3>(10));
			if (face->v[2]->node->EoL) f.segment<2>(cindexX) += (fme.segment<2>(13) + fie.segment<2>(13));

			// Diagonal x
			Matrix3d Mxx, Kxx;
			Kxx = Kme.block<3, 3>(0, 0); Mxx = Mie.block<3, 3>(0, 0);
			fillxMI(MDK_, M_, Kxx, Mxx, aindex, damping, h);
			Kxx = Kme.block<3, 3>(5, 5); Mxx = Mie.block<3, 3>(5, 5);
			fillxMI(MDK_, M_, Kxx, Mxx, bindex, damping, h);
			Kxx = Kme.block<3, 3>(10, 10); Mxx = Mie.block<3, 3>(10, 10);
			fillxMI(MDK_, M_, Kxx, Mxx, cindex, damping, h);

			// Off-Diagonal x
			Kxx = Kme.block<3, 3>(0, 5); Mxx = Mie.block<3, 3>(0, 5);
			fillxxMI(MDK_, M_, Kxx, Mxx, aindex, bindex, damping, h);
			Kxx = Kme.block<3, 3>(0, 10); Mxx = Mie.block<3, 3>(0, 10);
			fillxxMI(MDK_, M_, Kxx, Mxx, aindex, cindex, damping, h);
			Kxx = Kme.block<3, 3>(5, 10); Mxx = Mie.block<3, 3>(5, 10);
			fillxxMI(MDK_, M_, Kxx, Mxx, bindex, cindex, damping, h);

			// If has EOL componenent, Diagonal X
			Matrix2d MXX, KXX;
			if (face->v[0]->node->EoL) {
				KXX = Kme.block<2, 2>(3, 3); MXX = Mie.block<2, 2>(3, 3);
				fillXMI(MDK_, M_, KXX, MXX, aindexX, damping, h);
			}
			if (face->v[1]->node->EoL) {
				KXX = Kme.block<2, 2>(8, 8); MXX = Mie.block<2, 2>(8, 8);
				fillXMI(MDK_, M_, KXX, MXX, bindexX, damping, h);
			}
			if (face->v[2]->node->EoL) {
				KXX = Kme.block<2, 2>(13, 13); MXX = Mie.block<2, 2>(13, 13);
				fillXMI(MDK_, M_, KXX, MXX, cindexX, damping, h);
			}

			// If has EOL componenent, Off-Diagonal X
			if (face->v[0]->node->EoL && face->v[1]->node->EoL) {
				KXX = Kme.block<2, 2>(3, 8); MXX = Mie.block<2, 2>(3, 8);
				fillXXMI(MDK_, M_, KXX, MXX, aindexX, bindexX, damping, h);
			}
			if (face->v[0]->node->EoL && face->v[2]->node->EoL) {
				KXX = Kme.block<2, 2>(3, 13); MXX = Mie.block<2, 2>(3, 13);
				fillXXMI(MDK_, M_, KXX, MXX, aindexX, cindexX, damping, h);
			}
			if (face->v[1]->node->EoL && face->v[2]->node->EoL) {
				KXX = Kme.block<2, 2>(8, 13); MXX = Mie.block<2, 2>(8, 13);
				fillXXMI(MDK_, M_, KXX, MXX, bindexX, cindexX, damping, h);
			}

			// X-x values
			MatrixXd KXx, KxX, MXx, MxX;
			if (face->v[0]->node->EoL) {
				KXx = Kme.block<2, 3>(3, 0); MXx = Mie.block<2, 3>(3, 2);
				fillXxMI(MDK_, M_, KXx, MXx, aindexX, aindex, damping, h);
				KXx = Kme.block<2, 3>(3, 5); MXx = Mie.block<2, 3>(3, 5);
				fillXxMI(MDK_, M_, KXx, MXx, aindexX, bindex, damping, h);
				KXx = Kme.block<2, 3>(3, 10); MXx = Mie.block<2, 3>(3, 10);
				fillXxMI(MDK_, M_, KXx, MXx, aindexX, cindex, damping, h);

				//KxX = Kme.block<3, 2>(0, 3);
				//fillxXMI(MDK_, KxX, aindex, aindexX, damping, h);
			}
			if (face->v[1]->node->EoL) {
				//KXx = Kme.block<2, 3>(8, 0);
				//fillXxMI(MDK_, KXx, bindexX, aindex, damping, h);
				KXx = Kme.block<2, 3>(8, 5); MXx = Mie.block<2, 3>(8, 5);
				fillXxMI(MDK_, M_, KXx, MXx, bindexX, bindex, damping, h);
				KXx = Kme.block<2, 3>(8, 10); MXx = Mie.block<2, 3>(8, 10);
				fillXxMI(MDK_, M_, KXx, MXx, bindexX, cindex, damping, h);

				KxX = Kme.block<3, 2>(0, 8); MxX = Mie.block<3, 2>(0, 8);
				fillxXMI(MDK_, M_, KxX, MxX, aindex, bindexX, damping, h);
			}
			if (face->v[2]->node->EoL) {
				//KXx = Kme.block<2, 3>(13, 0);
				//fillXxMI(MDK_, KXx, cindexX, aindex, damping, h);
				//KXx = Kme.block<2, 3>(13, 5);
				//fillXxMI(MDK_, KXx, cindexX, bindex, damping, h);
				KXx = Kme.block<2, 3>(13, 10); MXx = Mie.block<2, 3>(13, 10);
				fillXxMI(MDK_, M_, KXx, MXx, cindexX, cindex, damping, h);

				KxX = Kme.block<3, 2>(0, 13); MxX = Mie.block<3, 2>(0, 13);
				fillxXMI(MDK_, M_, KxX, MxX, aindex, cindexX, damping, h);
				KxX = Kme.block<3, 2>(5, 13); MxX = Mie.block<3, 2>(5, 13);
				fillxXMI(MDK_, M_, KxX, MxX, bindex, cindexX, damping, h);
			}
		}
		else {

			f.segment<3>(aindex) += fme.segment<3>(0) + fie.segment<3>(0);
			f.segment<3>(bindex) += fme.segment<3>(3) + fie.segment<3>(3);
			f.segment<3>(cindex) += fme.segment<3>(6) + fie.segment<3>(6);

			Matrix3d Mxx, Kxx;
			Kxx = Kme.block<3, 3>(0, 0); Mxx = Mie.block<3, 3>(0, 0);
			fillxMI(MDK_, M_, Kxx, Mxx, aindex, damping, h);
			Kxx = Kme.block<3, 3>(3, 3); Mxx = Mie.block<3, 3>(3, 3);
			fillxMI(MDK_, M_, Kxx, Mxx, bindex, damping, h);
			Kxx = Kme.block<3, 3>(6, 6); Mxx = Mie.block<3, 3>(6, 6);
			fillxMI(MDK_, M_, Kxx, Mxx, cindex, damping, h);

			Kxx = Kme.block<3, 3>(0, 3); Mxx = Mie.block<3, 3>(0, 3);
			fillxxMI(MDK_, M_, Kxx, Mxx, aindex, bindex, damping, h);
			Kxx = Kme.block<3, 3>(0, 6); Mxx = Mie.block<3, 3>(0, 6);
			fillxxMI(MDK_, M_, Kxx, Mxx, aindex, cindex, damping, h);
			Kxx = Kme.block<3, 3>(3, 6); Mxx = Mie.block<3, 3>(3, 6);
			fillxxMI(MDK_, M_, Kxx, Mxx, bindex, cindex, damping, h);
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

void fillxXB(vector<T>& MDK_, const MatrixXd& KxX, int i0, int i1)
{
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 2; k++) {
			MDK_.push_back(T(i0 + j, i1 + k, KxX(j, k)));
			MDK_.push_back(T(i1 + k, i0 + j, KxX(j, k)));
		}
	}
}

void fillEOLBending(const Vert* v0, const Vert* v1, const Vert* v2, const Vert* v3, VectorXd& fb, MatrixXd& Kb)
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
		//Kb.block<2, 3>(jB, ja) = -Fb.transpose() * Kbba;
		Kb.block<2, 3>(jB, jb) = -Fb.transpose() * Kbbb;
		Kb.block<2, 3>(jB, jc) = -Fb.transpose() * Kbbc;
		Kb.block<2, 3>(jB, jd) = -Fb.transpose() * Kbbd;

		Kb.block<3, 2>(ja, jB) = -Kbab * Fb;
	}
	if (EOLC) {
		//Kb.block<2, 3>(jC, ja) = -Fc.transpose() * Kbca;
		//Kb.block<2, 3>(jC, jb) = -Fc.transpose() * Kbcb;
		Kb.block<2, 3>(jC, jc) = -Fc.transpose() * Kbcc;
		Kb.block<2, 3>(jC, jd) = -Fc.transpose() * Kbcd;

		Kb.block<3, 2>(ja, jC) = -Kbac * Fc;
		Kb.block<3, 2>(jb, jC) = -Kbbc * Fc;
	}
	if (EOLD) {
		//Kb.block<2, 3>(jD, ja) = -Fd.transpose() * Kbda;
		//Kb.block<2, 3>(jD, jb) = -Fd.transpose() * Kbdb;
		//Kb.block<2, 3>(jD, jc) = -Fd.transpose() * Kbdc;
		Kb.block<2, 3>(jD, jd) = -Fd.transpose() * Kbdd;

		Kb.block<3, 2>(ja, jD) = -Kbad * Fd;
		Kb.block<3, 2>(jb, jD) = -Kbbd * Fd;
		Kb.block<3, 2>(jc, jD) = -Kbcd * Fd;
	}
}

void edgeBasedF(const Mesh& mesh, const Material& mat, VectorXd& f, vector<T>& MDK_, double h)
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

			fillEOLBending(v0, v1, v2, v3, fbe, Kbe);
			
			f.segment<3>(aindex) += fbe.segment<3>(0);
			if (to_eolA) f.segment<2>(aindexX) += fbe.segment<2>(3);

			f.segment<3>(bindex) += fbe.segment<3>(5);
			if (to_eolB) f.segment<2>(bindexX) += fbe.segment<2>(8);

			f.segment<3>(cindex) += fbe.segment<3>(10);
			if (to_eolC) f.segment<2>(cindexX) += fbe.segment<2>(13);

			f.segment<3>(dindex) += fbe.segment<3>(15);
			if (to_eolD) f.segment<2>(dindexX) += fbe.segment<2>(18);

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

			MatrixXd KXx, KxX;
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
				//KXx = damping(1) * h * h * Kbe.block<2, 3>(8, 0);
				//fillXxB(MDK_, KXx, bindexX, aindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(8, 5);
				fillXxB(MDK_, KXx, bindexX, bindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(8, 10);
				fillXxB(MDK_, KXx, bindexX, cindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(8, 15);
				fillXxB(MDK_, KXx, bindexX, dindex);

				KxX = damping(1) * h * h * Kbe.block<3, 2>(0, 8);
				fillxXB(MDK_, KxX, aindex, bindexX);
			}
			if (to_eolC) {
				//KXx = damping(1) * h * h * Kbe.block<2, 3>(13, 0);
				//fillXxB(MDK_, KXx, cindexX, aindex);
				//KXx = damping(1) * h * h * Kbe.block<2, 3>(13, 5);
				//fillXxB(MDK_, KXx, cindexX, bindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(13, 10);
				fillXxB(MDK_, KXx, cindexX, cindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(13, 15);
				fillXxB(MDK_, KXx, cindexX, dindex);

				KxX = damping(1) * h * h * Kbe.block<3, 2>(0, 13);
				fillxXB(MDK_, KxX, aindex, cindexX);
				KxX = damping(1) * h * h * Kbe.block<3, 2>(5, 13);
				fillxXB(MDK_, KxX, bindex, cindexX);
			}
			if (to_eolD) {
				//KXx = damping(1) * h * h * Kbe.block<2, 3>(18, 0);
				//fillXxB(MDK_, KXx, dindexX, aindex);
				//KXx = damping(1) * h * h * Kbe.block<2, 3>(18, 5);
				//fillXxB(MDK_, KXx, dindexX, bindex);
				//KXx = damping(1) * h * h * Kbe.block<2, 3>(18, 10);
				//fillXxB(MDK_, KXx, dindexX, cindex);
				KXx = damping(1) * h * h * Kbe.block<2, 3>(18, 15);
				fillXxB(MDK_, KXx, dindexX, dindex);

				KxX = damping(1) * h * h * Kbe.block<3, 2>(0, 18);
				fillxXB(MDK_, KxX, aindex, dindexX);
				KxX = damping(1) * h * h * Kbe.block<3, 2>(5, 18);
				fillxXB(MDK_, KxX, bindex, dindexX);
				KxX = damping(1) * h * h * Kbe.block<3, 2>(10, 18);
				fillxXB(MDK_, KxX, cindex, dindexX);
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

	//nodeBasedF(mesh, f, M_, grav);
	faceBasedF(mesh, f, MDK_, M_, h);
	edgeBasedF(mesh, mat, f, MDK_, h);

	M.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2, mesh.nodes.size() * 3 + mesh.EoL_Count * 2);
	MDK.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2, mesh.nodes.size() * 3 + mesh.EoL_Count * 2);

	M.setFromTriplets(M_.begin(), M_.end());
	MDK.setFromTriplets(MDK_.begin(), MDK_.end());
}