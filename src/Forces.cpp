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

#include "external\ArcSim\util.hpp"

#include <iostream>
#include <utility>

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;

Matrix2d poldec(Matrix2d M) {
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

void faceBasedF(const Mesh& mesh, VectorXd& f, vector<T>& M_, vector<T>& MDK_, double h)
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

		double fi[9], fm[9], Mi[81], Km[81];

		VectorXd fie(15), fme(15);
		MatrixXd Mie(15, 15), Kme(15, 15);

		ComputeInertial(xa, xb, xc, Xa, Xb, Xc, g, mat->density, Wi, fi, Mi);
		ComputeMembrane(xa, xb, xc, Xa, Xb, Xc, mat->e, mat->nu, PP, QQ, Wm, fm, Km);

		fie.segment<9>(0) = Map<VectorXd>(fi, 9);
		fme.segment<9>(0) = Map<VectorXd>(fm, 9);
		Mie.block<9, 9>(0, 0) = Map<MatrixXd>(Mi, 9, 9);
		Kme.block<9, 9>(0, 0) = Map<MatrixXd>(Km, 9, 9);

		Vector2d damping(mat->dampingA, mat->dampingB);

		if (face->v[0]->node->EoL || face->v[1]->node->EoL || face->v[2]->node->EoL) {
			MatrixXd F = deform_grad(face);
			int ja = 0; int jA = 3; int jb = 5; int jB = 8; int jc = 10; int jC = 13;
			Matrix3d Kmaa = Kme.block(0, 0, 3, 3); Matrix3d Kmab = Kme.block(0, jb, 3, 3); Matrix3d Kmac = Kme.block(0, jc, 3, 3);
			Matrix3d Kmba = Kme.block(3, 0, 3, 3); Matrix3d Kmbb = Kme.block(3, 3, 3, 3); Matrix3d Kmbc = Kme.block(3, 6, 3, 3);
			Matrix3d Kmca = Kme.block(6, 0, 3, 3); Matrix3d Kmcb = Kme.block(6, 3, 3, 3); Matrix3d Kmcc = Kme.block(6, 6, 3, 3);
			Kme.block(ja, jA, 3, 2) = -Kmaa * F; Kme.block(ja, jB, 3, 2) = -Kmab * F; Kme.block(ja, jC, 3, 2) = -Kmac * F;
			Kme.block(jb, jA, 3, 2) = -Kmba * F; Kme.block(jb, jB, 3, 2) = -Kmbb * F; Kme.block(jb, jC, 3, 2) = -Kmbc * F;
			Kme.block(jc, jA, 3, 2) = -Kmca * F; Kme.block(jc, jB, 3, 2) = -Kmcb * F; Kme.block(jc, jC, 3, 2) = -Kmcc * F;
			Kme.block(jA, ja, 2, 3) = -F.transpose() * Kmaa; Kme.block(jA, jb, 2, 3) = -F.transpose() * Kmab; Kme.block(jA, jc, 2, 3) = -F.transpose() * Kmac;
			Kme.block(jB, ja, 2, 3) = -F.transpose() * Kmba; Kme.block(jB, jb, 2, 3) = -F.transpose() * Kmbb; Kme.block(jB, jc, 2, 3) = -F.transpose() * Kmbc;
			Kme.block(jC, ja, 2, 3) = -F.transpose() * Kmca; Kme.block(jC, jb, 2, 3) = -F.transpose() * Kmcb; Kme.block(jC, jc, 2, 3) = -F.transpose() * Kmcc;
			Kme.block(jA, jA, 2, 2) = F.transpose() * Kmaa * F; Kme.block(jA, jB, 2, 2) = F.transpose() * Kmab * F; Kme.block(jA, jC, 2, 2) = F.transpose() * Kmac * F;
			Kme.block(jB, jA, 2, 2) = F.transpose() * Kmba * F; Kme.block(jB, jB, 2, 2) = F.transpose() * Kmbb * F; Kme.block(jB, jC, 2, 2) = F.transpose() * Kmbc * F;
			Kme.block(jC, jA, 2, 2) = F.transpose() * Kmca * F; Kme.block(jC, jB, 2, 2) = F.transpose() * Kmcb * F; Kme.block(jC, jC, 2, 2) = F.transpose() * Kmcc * F;

			// x values
			Matrix3d Maa = Mie.block(0, 0, 3, 3);
			Matrix3d Kaa = Kme.block(0, 0, 3, 3);
			Matrix3d MDKaa = Maa + damping(0) * h * Maa + damping(1) * h * h * Kaa;
			f.segment<3>(face->v[0]->node->index * 3) += (fie.segment<3>(0) + fme.segment<3>(0));
			if (face->v[0]->node->EoL_index != -1) f.segment<2>(aindexX) += -F.transpose() * (fie.segment<3>(0) + fme.segment<3>(0));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(aindex + j, aindex + k, Maa(j, k)));
					MDK_.push_back(T(aindex + j, aindex + k, MDKaa(j, k)));
				}
			}

			Matrix3d Mbb = Mie.block(5, 5, 3, 3);
			Matrix3d Kbb = Kme.block(5, 5, 3, 3);
			Matrix3d MDKbb = Mbb + damping(0) * h * Mbb + damping(1) * h * h * Kbb;
			f.segment<3>(face->v[1]->node->index * 3) += (fie.segment<3>(5) + fme.segment<3>(5));
			if (face->v[1]->node->EoL_index != -1)f.segment<2>(bindexX) += -F.transpose() * (fie.segment<3>(5) + fme.segment<3>(5));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(bindex + j, bindex + k, Mbb(j, k)));
					MDK_.push_back(T(bindex + j, bindex + k, MDKbb(j, k)));
				}
			}

			Matrix3d Mcc = Mie.block(10, 10, 3, 3);
			Matrix3d Kcc = Kme.block(10, 10, 3, 3);
			Matrix3d MDKcc = Mcc + damping(0) * h * Mcc + damping(1) * h * h * Kcc;
			f.segment<3>(face->v[2]->node->index * 3) += (fie.segment<3>(10) + fme.segment<3>(10));
			if (face->v[2]->node->EoL_index != -1) f.segment<2>(cindexX) += -F.transpose() * (fie.segment<3>(10) + fme.segment<3>(10));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(cindex + j, cindex + k, Mcc(j, k)));
					MDK_.push_back(T(cindex + j, cindex + k, MDKcc(j, k)));
				}
			}

			Matrix3d Mab = Mie.block(0, 5, 3, 3);
			Matrix3d Kab = Kme.block(0, 5, 3, 3);
			Matrix3d MDKab = Mab + damping(0) * h * Mab + damping(1) * h * h * Kab;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(aindex + j, bindex + k, Mab(j, k)));
					M_.push_back(T(bindex + k, aindex + j, Mab(j, k)));
					MDK_.push_back(T(aindex + j, bindex + k, MDKab(j, k)));
					MDK_.push_back(T(bindex + k, aindex + j, MDKab(j, k)));
				}
			}

			Matrix3d Mac = Mie.block(0, 10, 3, 3);
			Matrix3d Kac = Kme.block(0, 10, 3, 3);
			Matrix3d MDKac = Mac + damping(0) * h * Mac + damping(1) * h * h * Kac;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(aindex + j, cindex + k, Mac(j, k)));
					M_.push_back(T(cindex + k, aindex + j, Mac(j, k)));
					MDK_.push_back(T(aindex + j, cindex + k, MDKac(j, k)));
					MDK_.push_back(T(cindex + k, aindex + j, MDKac(j, k)));
				}
			}

			Matrix3d Mbc = Mie.block(5, 10, 3, 3);
			Matrix3d Kbc = Kme.block(5, 10, 3, 3);
			Matrix3d MDKbc = Mbc + damping(0) * h * Mbc + damping(1) * h * h * Kbc;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(bindex + j, cindex + k, Mbc(j, k)));
					M_.push_back(T(cindex + k, bindex + j, Mbc(j, k)));
					MDK_.push_back(T(bindex + j, cindex + k, MDKbc(j, k)));
					MDK_.push_back(T(cindex + k, bindex + j, MDKbc(j, k)));
				}
			}

			// X values
			if (face->v[0]->node->EoL_index != -1) {
				Matrix2d MAA = Mie.block(3, 3, 2, 2);
				Matrix2d KAA = Kme.block(3, 3, 2, 2);
				Matrix2d MDKAA = MAA + damping(0) * h * MAA + damping(1) * h * h * KAA;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(aindexX + j, aindexX + k, MAA(j, k)));
						MDK_.push_back(T(aindexX + j, aindexX + k, MDKAA(j, k)));
					}
				}
			}

			if (face->v[1]->node->EoL_index != -1) {
				Matrix2d MBB = Mie.block(8, 8, 2, 2);
				Matrix2d KBB = Kme.block(8, 8, 2, 2);
				Matrix2d MDKBB = MBB + damping(0) * h * MBB + damping(1) * h * h * KBB;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(bindexX + j, bindexX + k, MBB(j, k)));
						MDK_.push_back(T(bindexX + j, bindexX + k, MDKBB(j, k)));
					}
				}
			}

			if (face->v[2]->node->EoL_index != -1) {
				Matrix2d MCC = Mie.block(13, 13, 2, 2);
				Matrix2d KCC = Kme.block(13, 13, 2, 2);
				Matrix2d MDKCC = MCC + damping(0) * h * MCC + damping(1) * h * h * KCC;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(cindexX + j, cindexX + k, MCC(j, k)));
						MDK_.push_back(T(cindexX + j, cindexX + k, MDKCC(j, k)));
					}
				}
			}

			if (face->v[0]->node->EoL_index != -1 && face->v[1]->node->EoL_index != -1) {
				Matrix2d MAB = Mie.block(3, 8, 2, 2);
				Matrix2d KAB = Kme.block(3, 8, 2, 2);
				Matrix2d MDKAB = MAB + damping(0) * h * MAB + damping(1) * h * h * KAB;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(aindexX + j, bindexX + k, MAB(j, k)));
						M_.push_back(T(bindexX + k, aindexX + j, MAB(j, k)));
						MDK_.push_back(T(aindexX + j, bindexX + k, MDKAB(j, k)));
						MDK_.push_back(T(bindexX + k, aindexX + j, MDKAB(j, k)));
					}
				}
			}

			if (face->v[0]->node->EoL_index != -1 && face->v[2]->node->EoL_index != -1) {
				Matrix2d MAC = Mie.block(3, 13, 2, 2);
				Matrix2d KAC = Kme.block(3, 13, 2, 2);
				Matrix2d MDKAC = MAC + damping(0) * h * MAC + damping(1) * h * h * KAC;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(aindexX + j, cindexX + k, MAC(j, k)));
						M_.push_back(T(cindexX + k, aindexX + j, MAC(j, k)));
						MDK_.push_back(T(aindexX + j, cindexX + k, MDKAC(j, k)));
						MDK_.push_back(T(cindexX + k, aindexX + j, MDKAC(j, k)));
					}
				}
			}

			if (face->v[1]->node->EoL_index != -1 && face->v[2]->node->EoL_index != -1) {
				Matrix2d MBC = Mie.block(8, 13, 2, 2);
				Matrix2d KBC = Kme.block(8, 13, 2, 2);
				Matrix2d MDKBC = MBC + damping(0) * h * MBC + damping(1) * h * h * KBC;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(bindexX + j, cindexX + k, MBC(j, k)));
						M_.push_back(T(cindexX + k, bindexX + j, MBC(j, k)));
						MDK_.push_back(T(bindexX + j, cindexX + k, MDKBC(j, k)));
						MDK_.push_back(T(cindexX + k, bindexX + j, MDKBC(j, k)));
					}
				}
			}

			// x-X values
			if (face->v[0]->node->EoL_index != -1) {
				MatrixXd MAa = Mie.block(3, 0, 2, 3);
				MatrixXd KAa = Kme.block(3, 0, 2, 3);
				MatrixXd MDKAa = MAa + damping(0) * h * MAa + damping(1) * h * h * KAa;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(aindexX + j, aindex + k, MAa(j, k)));
						M_.push_back(T(aindex + k, aindexX + j, MAa(j, k)));
						MDK_.push_back(T(aindexX + j, aindex + k, MDKAa(j, k)));
						MDK_.push_back(T(aindex + k, aindexX + j, MDKAa(j, k)));
					}
				}

				MatrixXd MAb = Mie.block(3, 5, 2, 3);
				MatrixXd KAb = Kme.block(3, 5, 2, 3);
				MatrixXd MDKAb = MAb + damping(0) * h * MAb + damping(1) * h * h * KAb;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(aindexX + j, bindex + k, MAb(j, k)));
						M_.push_back(T(bindex + k, aindexX + j, MAb(j, k)));
						MDK_.push_back(T(aindexX + j, bindex + k, MDKAb(j, k)));
						MDK_.push_back(T(bindex + k, aindexX + j, MDKAb(j, k)));
					}
				}

				MatrixXd MAc = Mie.block(3, 10, 2, 3);
				MatrixXd KAc = Kme.block(3, 10, 2, 3);
				MatrixXd MDKAc = MAc + damping(0) * h * MAc + damping(1) * h * h * KAc;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(aindexX + j, cindex + k, MAc(j, k)));
						M_.push_back(T(cindex + k, aindexX + j, MAc(j, k)));
						MDK_.push_back(T(aindexX + j, cindex + k, MDKAc(j, k)));
						MDK_.push_back(T(cindex + k, aindexX + j, MDKAc(j, k)));
					}
				}
			}

			if (face->v[1]->node->EoL_index != -1) {
				MatrixXd MBa = Mie.block(8, 0, 2, 3);
				MatrixXd KBa = Kme.block(8, 0, 2, 3);
				MatrixXd MDKBa = MBa + damping(0) * h * MBa + damping(1) * h * h * KBa;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(bindexX + j, aindex + k, MBa(j, k)));
						M_.push_back(T(aindex + k, bindexX + j, MBa(j, k)));
						MDK_.push_back(T(bindexX + j, aindex + k, MDKBa(j, k)));
						MDK_.push_back(T(aindex + k, bindexX + j, MDKBa(j, k)));
					}
				}

				MatrixXd MBb = Mie.block(8, 5, 2, 3);
				MatrixXd KBb = Kme.block(8, 5, 2, 3);
				MatrixXd MDKBb = MBb + damping(0) * h * MBb + damping(1) * h * h * KBb;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(bindexX + j, bindex + k, MBb(j, k)));
						M_.push_back(T(bindex + k, bindexX + j, MBb(j, k)));
						MDK_.push_back(T(bindexX + j, bindex + k, MDKBb(j, k)));
						MDK_.push_back(T(bindex + k, bindexX + j, MDKBb(j, k)));
					}
				}

				MatrixXd MBc = Mie.block(8, 10, 2, 3);
				MatrixXd KBc = Kme.block(8, 10, 2, 3);
				MatrixXd MDKBc = MBc + damping(0) * h * MBc + damping(1) * h * h * KBc;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(bindexX + j, cindex + k, MBc(j, k)));
						M_.push_back(T(cindex + k, bindexX + j, MBc(j, k)));
						MDK_.push_back(T(bindexX + j, cindex + k, MDKBc(j, k)));
						MDK_.push_back(T(cindex + k, bindexX + j, MDKBc(j, k)));
					}
				}
			}

			if (face->v[2]->node->EoL_index != -1) {
				MatrixXd MCa = Mie.block(13, 0, 2, 3);
				MatrixXd KCa = Kme.block(13, 0, 2, 3);
				MatrixXd MDKCa = MCa + damping(0) * h * MCa + damping(1) * h * h * KCa;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(cindexX + j, aindex + k, MCa(j, k)));
						M_.push_back(T(aindex + k, cindexX + j, MCa(j, k)));
						MDK_.push_back(T(cindexX + j, aindex + k, MDKCa(j, k)));
						MDK_.push_back(T(aindex + k, cindexX + j, MDKCa(j, k)));
					}
				}

				MatrixXd MCb = Mie.block(13, 5, 2, 3);
				MatrixXd KCb = Kme.block(13, 5, 2, 3);
				MatrixXd MDKCb = MCb + damping(0) * h * MCb + damping(1) * h * h * KCb;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(cindexX + j, bindex + k, MCb(j, k)));
						M_.push_back(T(bindex + k, cindexX + j, MCb(j, k)));
						MDK_.push_back(T(cindexX + j, bindex + k, MDKCb(j, k)));
						MDK_.push_back(T(bindex + k, cindexX + j, MDKCb(j, k)));
					}
				}

				MatrixXd MCc = Mie.block(13, 10, 2, 3);
				MatrixXd KCc = Kme.block(13, 10, 2, 3);
				MatrixXd MDKCc = MCc + damping(0) * h * MCc + damping(1) * h * h * KCc;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						M_.push_back(T(cindexX + j, cindex + k, MCc(j, k)));
						M_.push_back(T(cindex + k, cindexX + j, MCc(j, k)));
						MDK_.push_back(T(cindexX + j, cindex + k, MDKCc(j, k)));
						MDK_.push_back(T(cindex + k, cindexX + j, MDKCc(j, k)));
					}
				}
			}
		}
		else {
			Matrix3d Maa = Mie.block(0, 0, 3, 3);
			Matrix3d Kaa = Kme.block(0, 0, 3, 3);
			Matrix3d MDKaa = Maa + damping(0) * h * Maa + damping(1) * h * h * Kaa;
			f.segment<3>(aindex) += (fie.segment<3>(0) + fme.segment<3>(0));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(aindex + j, aindex + k, Maa(j, k)));
					MDK_.push_back(T(aindex + j, aindex + k, MDKaa(j, k)));
				}
			}

			Matrix3d Mbb = Mie.block(3, 3, 3, 3);
			Matrix3d Kbb = Kme.block(3, 3, 3, 3);
			Matrix3d MDKbb = Mbb + damping(0) * h * Mbb + damping(1) * h * h * Kbb;
			f.segment<3>(bindex) += (fie.segment<3>(3) + fme.segment<3>(3));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(bindex + j, bindex + k, Mbb(j, k)));
					MDK_.push_back(T(bindex + j, bindex + k, MDKbb(j, k)));
				}
			}

			Matrix3d Mcc = Mie.block(6, 6, 3, 3);
			Matrix3d Kcc = Kme.block(6, 6, 3, 3);
			Matrix3d MDKcc = Mcc + damping(0) * h * Mcc + damping(1) * h * h * Kcc;
			f.segment<3>(cindex) += (fie.segment<3>(6) + fme.segment<3>(6));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(cindex + j, cindex + k, Mcc(j, k)));
					MDK_.push_back(T(cindex + j, cindex + k, MDKcc(j, k)));
				}
			}

			Matrix3d Mab = Mie.block(0, 3, 3, 3);
			Matrix3d Kab = Kme.block(0, 3, 3, 3);
			Matrix3d MDKab = Mab + damping(0) * h * Mab + damping(1) * h * h * Kab;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(aindex + j, bindex + k, Mab(j, k)));
					M_.push_back(T(bindex + k, aindex + j, Mab(j, k)));
					MDK_.push_back(T(aindex + j, bindex + k, MDKab(j, k)));
					MDK_.push_back(T(bindex + k, aindex + j, MDKab(j, k)));
				}
			}

			Matrix3d Mac = Mie.block(0, 6, 3, 3);
			Matrix3d Kac = Kme.block(0, 6, 3, 3);
			Matrix3d MDKac = Mac + damping(0) * h * Mac + damping(1) * h * h * Kac;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(aindex + j, cindex + k, Mac(j, k)));
					M_.push_back(T(cindex + k, aindex + j, Mac(j, k)));
					MDK_.push_back(T(aindex + j, cindex + k, MDKac(j, k)));
					MDK_.push_back(T(cindex + k, aindex + j, MDKac(j, k)));
				}
			}

			Matrix3d Mbc = Mie.block(3, 6, 3, 3);
			Matrix3d Kbc = Kme.block(3, 6, 3, 3);
			Matrix3d MDKbc = Mbc + damping(0) * h * Mbc + damping(1) * h * h * Kbc;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(bindex + j, cindex + k, Mbc(j, k)));
					M_.push_back(T(cindex + k, bindex + j, Mbc(j, k)));
					MDK_.push_back(T(bindex + j, cindex + k, MDKbc(j, k)));
					MDK_.push_back(T(cindex + k, bindex + j, MDKbc(j, k)));
				}
			}
		}
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

			// x values
			Matrix3d K00 = damping(1) * h * h * Kbe.block(0, 0, 3, 3);
			f.segment<3>(aindex) += fbe.segment<3>(0);
			if (to_eolA) f.segment<2>(aindexX) += -F.transpose() * fbe.segment<3>(0);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(aindex + j, aindex + k, K00(j, k)));
				}
			}

			Matrix3d K11 = damping(1) * h * h * Kbe.block(5, 5, 3, 3);
			f.segment<3>(bindex) += fbe.segment<3>(5);
			if (to_eolB) f.segment<2>(bindexX) += -F.transpose() * fbe.segment<3>(5);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(bindex + j, bindex + k, K11(j, k)));
				}
			}

			Matrix3d K22 = damping(1) * h * h * Kbe.block(10, 10, 3, 3);
			f.segment<3>(cindex) += fbe.segment<3>(10);
			if (to_eolC) f.segment<2>(cindexX) += -F.transpose() * fbe.segment<3>(10);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(cindex + j, cindex + k, K22(j, k)));
				}
			}

			Matrix3d K33 = damping(1) * h * h * Kbe.block(15, 15, 3, 3);
			f.segment<3>(dindex) += fbe.segment<3>(15);
			if (to_eolD) f.segment<2>(dindexX) += -F.transpose() * fbe.segment<3>(15);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(dindex + j, dindex + k, K33(j, k)));
				}
			}

			Matrix3d K01 = damping(1) * h * h * Kbe.block(0, 5, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(aindex + j, bindex + k, K01(j, k)));
					MDK_.push_back(T(bindex + k, aindex + j, K01(j, k)));
				}
			}

			Matrix3d K02 = damping(1) * h * h * Kbe.block(0, 10, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(aindex + j, cindex + k, K02(j, k)));
					MDK_.push_back(T(cindex + k, aindex + j, K02(j, k)));
				}
			}

			Matrix3d K03 = damping(1) * h * h * Kbe.block(0, 15, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(aindex + j, dindex + k, K03(j, k)));
					MDK_.push_back(T(dindex + k, aindex + j, K03(j, k)));
				}
			}

			Matrix3d K12 = damping(1) * h * h * Kbe.block(5, 10, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(bindex + j, cindex + k, K12(j, k)));
					MDK_.push_back(T(cindex + k, bindex + j, K12(j, k)));
				}
			}

			Matrix3d K13 = damping(1) * h * h * Kbe.block(5, 15, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(bindex + j, dindex + k, K13(j, k)));
					MDK_.push_back(T(dindex + k, bindex + j, K13(j, k)));
				}
			}

			Matrix3d K23 = damping(1) * h * h * Kbe.block(10, 15, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(cindex + j, dindex + k, K23(j, k)));
					MDK_.push_back(T(dindex + k, cindex + j, K23(j, k)));
				}
			}

			// X values
			if (to_eolA) {
				Matrix2d K00X = damping(1) * h * h * Kbe.block(3, 3, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(aindexX + j, aindexX + k, K00X(j, k)));
					}
				}
			}

			if (to_eolB) {
				Matrix2d K11X = damping(1) * h * h * Kbe.block(8, 8, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(bindexX + j, bindexX + k, K11X(j, k)));
					}
				}
			}

			if (to_eolC) {
				Matrix2d K22X = damping(1) * h * h * Kbe.block(13, 13, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(cindexX + j, cindexX + k, K22X(j, k)));
					}
				}
			}

			if (to_eolD) {
				Matrix2d K33X = damping(1) * h * h * Kbe.block(18, 18, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(dindexX + j, dindexX + k, K33X(j, k)));
					}
				}
			}

			if (to_eolA && to_eolB) {
				Matrix2d K01X = damping(1) * h * h * Kbe.block(3, 8, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(aindexX + j, bindexX + k, K01X(j, k)));
						MDK_.push_back(T(bindexX + k, aindexX + j, K01X(j, k)));
					}
				}
			}

			if (to_eolA && to_eolC) {
				Matrix2d K02X = damping(1) * h * h * Kbe.block(3, 13, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(aindexX + j, cindexX + k, K02X(j, k)));
						MDK_.push_back(T(cindexX + k, aindexX + j, K02X(j, k)));
					}
				}
			}

			if (to_eolA && to_eolD) {
				Matrix2d K03X = damping(1) * h * h * Kbe.block(3, 18, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(aindexX + j, dindexX + k, K03X(j, k)));
						MDK_.push_back(T(dindexX + k, aindexX + j, K03X(j, k)));
					}
				}
			}

			if (to_eolB && to_eolC) {
				Matrix2d K12X = damping(1) * h * h * Kbe.block(8, 13, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(bindexX + j, cindexX + k, K12X(j, k)));
						MDK_.push_back(T(cindexX + k, bindexX + j, K12X(j, k)));
					}
				}
			}

			if (to_eolB && to_eolD) {
				Matrix2d K13X = damping(1) * h * h * Kbe.block(8, 18, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(bindexX + j, dindexX + k, K13X(j, k)));
						MDK_.push_back(T(dindexX + k, bindexX + j, K13X(j, k)));
					}
				}
			}

			if (to_eolC && to_eolD) {
				Matrix2d K23X = damping(1) * h * h * Kbe.block(13, 18, 2, 2);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(cindexX + j, dindexX + k, K23X(j, k)));
						MDK_.push_back(T(dindexX + k, cindexX + j, K23X(j, k)));
					}
				}
			}

			// x-X values
			if (to_eolA) {
				MatrixXd KAa = damping(1) * h * h * Kbe.block(3, 0, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(aindexX + j, aindex + k, KAa(j, k)));
						MDK_.push_back(T(aindex + k, aindexX + j, KAa(j, k)));
					}
				}

				MatrixXd KAb = damping(1) * h * h * Kbe.block(3, 5, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(aindexX + j, bindex + k, KAb(j, k)));
						MDK_.push_back(T(bindex + k, aindexX + j, KAb(j, k)));
					}
				}

				MatrixXd KAc = damping(1) * h * h * Kbe.block(3, 10, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(aindexX + j, cindex + k, KAc(j, k)));
						MDK_.push_back(T(cindex + k, aindexX + j, KAc(j, k)));
					}
				}

				MatrixXd KAd = damping(1) * h * h * Kbe.block(3, 15, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(aindexX + j, dindex + k, KAd(j, k)));
						MDK_.push_back(T(dindex + k, aindexX + j, KAd(j, k)));
					}
				}
			}

			if (to_eolB) {
				MatrixXd KBa = damping(1) * h * h * Kbe.block(8, 0, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(bindexX + j, aindex + k, KBa(j, k)));
						MDK_.push_back(T(aindex + k, bindexX + j, KBa(j, k)));
					}
				}

				MatrixXd KBb = damping(1) * h * h * Kbe.block(8, 5, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(bindexX + j, bindex + k, KBb(j, k)));
						MDK_.push_back(T(bindex + k, bindexX + j, KBb(j, k)));
					}
				}

				MatrixXd KBc = damping(1) * h * h * Kbe.block(8, 10, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(bindexX + j, cindex + k, KBc(j, k)));
						MDK_.push_back(T(cindex + k, bindexX + j, KBc(j, k)));
					}
				}

				MatrixXd KBd = damping(1) * h * h * Kbe.block(8, 15, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(bindexX + j, dindex + k, KBc(j, k)));
						MDK_.push_back(T(dindex + k, bindexX + j, KBc(j, k)));
					}
				}
			}

			if (to_eolC) {
				MatrixXd KCa = damping(1) * h * h * Kbe.block(13, 0, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(cindexX + j, aindex + k, KCa(j, k)));
						MDK_.push_back(T(aindex + k, cindexX + j, KCa(j, k)));
					}
				}

				MatrixXd KCb = damping(1) * h * h * Kbe.block(13, 5, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(cindexX + j, bindex + k, KCb(j, k)));
						MDK_.push_back(T(bindex + k, cindexX + j, KCb(j, k)));
					}
				}

				MatrixXd KCc = damping(1) * h * h * Kbe.block(13, 10, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(cindexX + j, cindex + k, KCc(j, k)));
						MDK_.push_back(T(cindex + k, cindexX + j, KCc(j, k)));
					}
				}

				MatrixXd KCd = damping(1) * h * h * Kbe.block(13, 15, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(cindexX + j, dindex + k, KCd(j, k)));
						MDK_.push_back(T(dindex + k, cindexX + j, KCd(j, k)));
					}
				}
			}

			if (to_eolD) {
				MatrixXd KDa = damping(1) * h * h * Kbe.block(18, 0, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(dindexX + j, aindex + k, KDa(j, k)));
						MDK_.push_back(T(aindex + k, dindexX + j, KDa(j, k)));
					}
				}

				MatrixXd KDb = damping(1) * h * h * Kbe.block(18, 5, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(dindexX + j, bindex + k, KDb(j, k)));
						MDK_.push_back(T(bindex + k, dindexX + j, KDb(j, k)));
					}
				}

				MatrixXd KDc = damping(1) * h * h * Kbe.block(18, 10, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(dindexX + j, cindex + k, KDc(j, k)));
						MDK_.push_back(T(cindex + k, dindexX + j, KDc(j, k)));
					}
				}

				MatrixXd KDd = damping(1) * h * h * Kbe.block(18, 15, 2, 3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 3; k++) {
						MDK_.push_back(T(dindexX + j, dindex + k, KDd(j, k)));
						MDK_.push_back(T(dindex + k, dindexX + j, KDd(j, k)));
					}
				}
			}
		}
		else {
			Matrix3d K00 = damping(1) * h * h * Kbe.block(0, 0, 3, 3);
			f.segment<3>(aindex) += fbe.segment<3>(0);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(aindex + j, aindex + k, K00(j, k)));
				}
			}

			Matrix3d K11 = damping(1) * h * h * Kbe.block(3, 3, 3, 3);
			f.segment<3>(bindex) += fbe.segment<3>(3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(bindex + j, bindex + k, K11(j, k)));
				}
			}

			Matrix3d K22 = damping(1) * h * h * Kbe.block(6, 6, 3, 3);
			f.segment<3>(cindex) += fbe.segment<3>(6);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(cindex + j, cindex + k, K22(j, k)));
				}
			}

			Matrix3d K33 = damping(1) * h * h * Kbe.block(9, 9, 3, 3);
			f.segment<3>(dindex) += fbe.segment<3>(9);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(dindex + j, dindex + k, K33(j, k)));
				}
			}

			Matrix3d K01 = damping(1) * h * h * Kbe.block(0, 3, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(aindex + j, bindex + k, K01(j, k)));
					MDK_.push_back(T(bindex + k, aindex + j, K01(j, k)));
				}
			}

			Matrix3d K02 = damping(1) * h * h * Kbe.block(0, 6, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(aindex + j, cindex + k, K02(j, k)));
					MDK_.push_back(T(cindex + k, aindex + j, K02(j, k)));
				}
			}

			Matrix3d K03 = damping(1) * h * h * Kbe.block(0, 9, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(aindex + j, dindex + k, K03(j, k)));
					MDK_.push_back(T(dindex + k, aindex + j, K03(j, k)));
				}
			}

			Matrix3d K12 = damping(1) * h * h * Kbe.block(3, 6, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(bindex + j, cindex + k, K12(j, k)));
					MDK_.push_back(T(cindex + k, bindex + j, K12(j, k)));
				}
			}

			Matrix3d K13 = damping(1) * h * h * Kbe.block(3, 9, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(bindex + j, dindex + k, K13(j, k)));
					MDK_.push_back(T(dindex + k, bindex + j, K13(j, k)));
				}
			}

			Matrix3d K23 = damping(1) * h * h * Kbe.block(6, 9, 3, 3);
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					MDK_.push_back(T(cindex + j, dindex + k, K23(j, k)));
					MDK_.push_back(T(dindex + k, cindex + j, K23(j, k)));
				}
			}
		}
	}
}

void Forces::fill(const Mesh& mesh, const Material& mat, double h)
{
	f.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2);
	f.setZero();
	vector<T> M_;
	vector<T> MDK_;

	faceBasedF(mesh, f, M_, MDK_, h);
	edgeBasedF(mesh, mat, f, M_, MDK_, h);

	M.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2, mesh.nodes.size() * 3 + mesh.EoL_Count * 2);
	MDK.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2, mesh.nodes.size() * 3 + mesh.EoL_Count * 2);

	M.setFromTriplets(M_.begin(), M_.end());
	MDK.setFromTriplets(MDK_.begin(), MDK_.end());
}