#include <iostream>

#include <fstream>

#include "Cloth.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include "Box.h"
#include "boxTriCollision.h"
#include "MatlabDebug.h"

// FEM calculations
#include "ComputeInertialLAG.h"
#include "ComputeMembraneLAG.h"
#include "ComputeBendingLAG.h"

// Mosek-libigl
#include "igl\mosek\mosek_quadprog.h"
#include "mosek.h"

// ArcSim
#include "mesh.hpp"
#include "vectors.hpp"
#include "util.hpp"
#include "io.hpp"
#include "geometry.hpp"
#include "remesh.hpp"

using namespace std;
using namespace Eigen;

void stepL(shared_ptr<Cloth> cloth, const Vector3d &grav, shared_ptr<Box> box, double t, double h)
{
	Mesh mesh = cloth->mesh;

	cloth->v.resize(mesh.nodes.size() * 3);
	f.resize(mesh.nodes.size() * 3);

	v.setZero();
	f.setZero();

	//double smallest_area = 1e8;
	//double smallest_aspect_ratio = 1e8;

	int collision = 0; // TODO:: Move

	bool fixed_points = false;  // TODO:: Restructure

	vector<T> M_;
	vector<T> MDK_;
	//vector<T> N_;
	N_.clear();
	int sparse_size = 0;

	// Box Setup
	collisions.clear();
	MatrixXd verts2(3, mesh.nodes.size());
	MatrixXi faces2(3, mesh.faces.size());

	for (int i = 0; i < mesh.nodes.size(); i++) {
		// For collision
		if (coll) verts2.col(i) = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);

		// Fill vertex vector
		v(3 * i) = mesh.nodes[i]->v[0];
		v(3 * i + 1) = mesh.nodes[i]->v[1];
		v(3 * i + 2) = mesh.nodes[i]->v[2];
		sparse_size += 3;

		int tmpe = 0;
		for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
			if (mesh.nodes[i]->adje[j]->preserve) tmpe++;
		}
		if (tmpe > 2) {
			cout << "to many" << endl;
		}
	}

	for (int i = 0; i < mesh.faces.size(); i++) {

		// For Collision
		if (coll) faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
		//if (area(mesh.faces[i]) < smallest_area) smallest_area = area(mesh.faces[i]);
		//if (aspect(mesh.faces[i]) < smallest_aspect_ratio) smallest_aspect_ratio = aspect(mesh.faces[i]);

		// Set up corrdaintes
		double xa[3];
		double xb[3];
		double xc[3];
		double Xa[2];
		double Xb[2];
		double Xc[2];

		// Fill temporary vectors 
		// TODO:: Convrsion functions
		// TODO:: Likely some negation problems with ordering of vertices
		Vector3d txa(3);
		Vector3d txb(3);
		Vector3d txc(3);
		Vector2d tXa(2);
		Vector2d tXb(2);
		Vector2d tXc(2);
		txa << mesh.faces[i]->v[0]->node->x[0], mesh.faces[i]->v[0]->node->x[1], mesh.faces[i]->v[0]->node->x[2];
		txb << mesh.faces[i]->v[1]->node->x[0], mesh.faces[i]->v[1]->node->x[1], mesh.faces[i]->v[1]->node->x[2];
		txc << mesh.faces[i]->v[2]->node->x[0], mesh.faces[i]->v[2]->node->x[1], mesh.faces[i]->v[2]->node->x[2];
		tXa << mesh.faces[i]->v[0]->u[0], mesh.faces[i]->v[0]->u[1];
		tXb << mesh.faces[i]->v[1]->u[0], mesh.faces[i]->v[1]->u[1];
		tXc << mesh.faces[i]->v[2]->u[0], mesh.faces[i]->v[2]->u[1];

		// Inertial
		double g[3];
		double Wi[1]; // Gravitational potential energy
		double fi[9]; // Gravity force vector
		double Mi[81]; //Inertial matrix
		Map<Vector3d>(xa, 3) = txa;
		Map<Vector3d>(xb, 3) = txb;
		Map<Vector3d>(xc, 3) = txc;
		Map<Vector2d>(Xa, 2) = tXa;
		Map<Vector2d>(Xb, 2) = tXb;
		Map<Vector2d>(Xc, 2) = tXc;
		Map<Vector3d>(g, grav.rows(), grav.cols()) = grav;
		ComputeInertiaLAG(xa, xb, xc, Xa, Xb, Xc, g, mesh.faces[i]->m, Wi, fi, Mi); // INVERSE
		VectorXd fie = Map<VectorXd>(fi, 9); // INVERSE
		MatrixXd Mie = Map<MatrixXd>(Mi, 9, 9);	// INVERSE



												// Membrane
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

		double PP[6];
		double QQ[4];
		Map<MatrixXd>(PP, Pm.rows(), Pm.cols()) = Pm;
		Map<Matrix2d>(QQ, Qm.rows(), Qm.cols()) = Qm;

		double Wm[1];
		double fm[9];
		double Km[81];
		ComputeMembraneLAG(xa, xb, xc, Xa, Xb, Xc, e, nu, PP, QQ, Wm, fm, Km);
		VectorXd fme = Map<VectorXd>(fm, 9); // INVERSE
		MatrixXd Kme = Map<MatrixXd>(Km, 9, 9); // INVERSE

		Matrix3d Maa = Mie.block(0, 0, 3, 3);
		Matrix3d Kaa = Kme.block(0, 0, 3, 3);
		Matrix3d MDKaa = Maa + damping(0) * h * Maa + damping(1) * h * h * Kaa;
		f.segment<3>(mesh.faces[i]->v[0]->node->index * 3) += (fie.segment<3>(0) + fme.segment<3>(0));
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				M_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[0]->node->index * 3 + k, Maa(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[0]->node->index * 3 + k, MDKaa(j, k)));
			}
		}

		Matrix3d Mbb = Mie.block(3, 3, 3, 3);
		Matrix3d Kbb = Kme.block(3, 3, 3, 3);
		Matrix3d MDKbb = Mbb + damping(0) * h * Mbb + damping(1) * h * h * Kbb;
		f.segment<3>(mesh.faces[i]->v[1]->node->index * 3) += (fie.segment<3>(3) + fme.segment<3>(3));
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				M_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + j, mesh.faces[i]->v[1]->node->index * 3 + k, Mbb(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + j, mesh.faces[i]->v[1]->node->index * 3 + k, MDKbb(j, k)));
			}
		}

		Matrix3d Mcc = Mie.block(6, 6, 3, 3);
		Matrix3d Kcc = Kme.block(6, 6, 3, 3);
		Matrix3d MDKcc = Mcc + damping(0) * h * Mcc + damping(1) * h * h * Kcc;
		f.segment<3>(mesh.faces[i]->v[2]->node->index * 3) += (fie.segment<3>(6) + fme.segment<3>(6));
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				M_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, Mcc(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, MDKcc(j, k)));
			}
		}

		Matrix3d Mab = Mie.block(0, 3, 3, 3);
		Matrix3d Kab = Kme.block(0, 3, 3, 3);
		Matrix3d MDKab = Mab + damping(0) * h * Mab + damping(1) * h * h * Kab;
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				M_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[1]->node->index * 3 + k, Mab(j, k)));
				M_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + k, mesh.faces[i]->v[0]->node->index * 3 + j, Mab(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[1]->node->index * 3 + k, MDKab(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + k, mesh.faces[i]->v[0]->node->index * 3 + j, MDKab(j, k)));
			}
		}

		Matrix3d Mac = Mie.block(0, 6, 3, 3);
		Matrix3d Kac = Kme.block(0, 6, 3, 3);
		Matrix3d MDKac = Mac + damping(0) * h * Mac + damping(1) * h * h * Kac;
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				M_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, Mac(j, k)));
				M_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + k, mesh.faces[i]->v[0]->node->index * 3 + j, Mac(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, MDKac(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + k, mesh.faces[i]->v[0]->node->index * 3 + j, MDKac(j, k)));
			}
		}

		Matrix3d Mbc = Mie.block(3, 6, 3, 3);
		Matrix3d Kbc = Kme.block(3, 6, 3, 3);
		Matrix3d MDKbc = Mbc + damping(0) * h * Mbc + damping(1) * h * h * Kbc;
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				M_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, Mbc(j, k)));
				M_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + k, mesh.faces[i]->v[1]->node->index * 3 + j, Mbc(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, MDKbc(j, k)));
				MDK_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + k, mesh.faces[i]->v[1]->node->index * 3 + j, MDKbc(j, k)));
			}
		}
	}

	//cout << "smallest area: " << smallest_area << endl;
	//cout << "smallest aspect ratio: " << smallest_aspect_ratio << endl;

	for (int i = 0; i < mesh.edges.size(); i++) {

		if (mesh.edges[i]->adjf[0] == NULL || mesh.edges[i]->adjf[1] == NULL) {
			continue;
		}

		// Set up corrdaintes
		double xa[3];
		double xb[3];
		double xc[3];
		double xd[3];
		double Xa[2];
		double Xb[2];
		double Xc[2];
		double Xd[2];

		// Keeps track of which face node so we don't have to traverse faces multiple times
		int aindex;
		int bindex;
		int cindex;
		int dindex;

		// Fill temporary vectors 
		// TODO:: This is a mess
		Vector3d txa(3);
		Vector3d txb(3);
		Vector3d txc(3);
		Vector3d txd(3);
		Vector2d tXa(2);
		Vector2d tXb(2);
		Vector2d tXc(2);
		Vector2d tXd(2);
		txa << mesh.edges[i]->n[0]->x[0], mesh.edges[i]->n[0]->x[1], mesh.edges[i]->n[0]->x[2];
		aindex = mesh.edges[i]->n[0]->index * 3;
		txb << mesh.edges[i]->n[1]->x[0], mesh.edges[i]->n[1]->x[1], mesh.edges[i]->n[1]->x[2];
		bindex = mesh.edges[i]->n[1]->index * 3;
		for (int j = 0; j < 3; j++) {
			if (mesh.edges[i]->adjf[0]->v[j]->node->x != mesh.edges[i]->n[0]->x && mesh.edges[i]->adjf[0]->v[j]->node->x != mesh.edges[i]->n[1]->x) {
				txc << mesh.edges[i]->adjf[0]->v[j]->node->x[0], mesh.edges[i]->adjf[0]->v[j]->node->x[1], mesh.edges[i]->adjf[0]->v[j]->node->x[2];
				cindex = mesh.edges[i]->adjf[0]->v[j]->node->index * 3;
			}
		}
		for (int j = 0; j < 3; j++) {
			if (mesh.edges[i]->adjf[1]->v[j]->node->x != mesh.edges[i]->n[0]->x && mesh.edges[i]->adjf[1]->v[j]->node->x != mesh.edges[i]->n[1]->x) {
				txd << mesh.edges[i]->adjf[1]->v[j]->node->x[0], mesh.edges[i]->adjf[1]->v[j]->node->x[1], mesh.edges[i]->adjf[1]->v[j]->node->x[2];
				dindex = mesh.edges[i]->adjf[1]->v[j]->node->index * 3;
			}
		}

		tXa << mesh.edges[i]->n[0]->verts[0]->u[0], mesh.edges[i]->n[0]->verts[0]->u[1];
		tXb << mesh.edges[i]->n[1]->verts[0]->u[0], mesh.edges[i]->n[1]->verts[0]->u[1];
		for (int j = 0; j < 3; j++) {
			if (mesh.edges[i]->adjf[0]->v[j]->u != mesh.edges[i]->n[0]->verts[0]->u && mesh.edges[i]->adjf[0]->v[j]->u != mesh.edges[i]->n[1]->verts[0]->u) {
				tXc << mesh.edges[i]->adjf[0]->v[j]->u[0], mesh.edges[i]->adjf[0]->v[j]->u[1];
			}
		}
		for (int j = 0; j < 3; j++) {
			if (mesh.edges[i]->adjf[1]->v[j]->u != mesh.edges[i]->n[0]->verts[0]->u && mesh.edges[i]->adjf[1]->v[j]->u != mesh.edges[i]->n[1]->verts[0]->u) {
				tXd << mesh.edges[i]->adjf[1]->v[j]->u[0], mesh.edges[i]->adjf[1]->v[j]->u[1];
			}
		}

		// Bending
		double Wb[1]; // Bending potential energy
		double fb[12]; // Bending force vector
		double Kb[144]; //Bending stiffness matrix
		Map<Vector3d>(xa, 3) = txa;
		Map<Vector3d>(xb, 3) = txb;
		Map<Vector3d>(xc, 3) = txc;
		Map<Vector3d>(xd, 3) = txd;
		Map<Vector2d>(Xa, 2) = tXa;
		Map<Vector2d>(Xb, 2) = tXb;
		Map<Vector2d>(Xc, 2) = tXc;
		Map<Vector2d>(Xd, 2) = tXd;
		ComputeBendingLAG(xa, xb, xc, xd, Xa, Xb, Xc, Xd, beta, Wb, fb, Kb);
		VectorXd fbe = Map<VectorXd>(fb, 12); // INVERSE
		MatrixXd Kbe = Map<MatrixXd>(Kb, 12, 12); // INVERSE

		Matrix3d K00 = h * h * Kbe.block(0, 0, 3, 3);
		f.segment<3>(aindex) += fbe.segment<3>(0);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(aindex + j, aindex + k, K00(j, k)));
			}
		}

		Matrix3d K11 = h * h * Kbe.block(3, 3, 3, 3);
		f.segment<3>(bindex) += fbe.segment<3>(3);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(bindex + j, bindex + k, K11(j, k)));
			}
		}

		Matrix3d K22 = h * h * Kbe.block(6, 6, 3, 3);
		f.segment<3>(cindex) += fbe.segment<3>(6);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(cindex + j, cindex + k, K22(j, k)));
			}
		}

		Matrix3d K33 = h * h * Kbe.block(9, 9, 3, 3);
		f.segment<3>(dindex) += fbe.segment<3>(9);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(dindex + j, dindex + k, K33(j, k)));
			}
		}

		Matrix3d K01 = h * h * Kbe.block(0, 3, 3, 3);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(aindex + j, bindex + k, K01(j, k)));
				MDK_.push_back(T(bindex + k, aindex + j, K01(j, k)));
			}
		}

		Matrix3d K02 = h * h * Kbe.block(0, 6, 3, 3);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(aindex + j, cindex + k, K02(j, k)));
				MDK_.push_back(T(cindex + k, aindex + j, K02(j, k)));
			}
		}

		Matrix3d K03 = h * h * Kbe.block(0, 9, 3, 3);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(aindex + j, dindex + k, K03(j, k)));
				MDK_.push_back(T(dindex + k, aindex + j, K03(j, k)));
			}
		}

		Matrix3d K12 = h * h * Kbe.block(3, 6, 3, 3);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(bindex + j, cindex + k, K12(j, k)));
				MDK_.push_back(T(cindex + k, bindex + j, K12(j, k)));
			}
		}

		Matrix3d K13 = h * h * Kbe.block(3, 9, 3, 3);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(bindex + j, dindex + k, K13(j, k)));
				MDK_.push_back(T(dindex + k, bindex + j, K13(j, k)));
			}
		}

		Matrix3d K23 = h * h * Kbe.block(6, 9, 3, 3);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				MDK_.push_back(T(cindex + j, dindex + k, K23(j, k)));
				MDK_.push_back(T(dindex + k, cindex + j, K23(j, k)));
			}
		}
	}

	//if (matlab_debug_collision) {
	//	mat_to_file(verts2, "verts2");
	//	VectorXi vvv(3);
	//	vvv << 1, 1, 1;
	//	mat_to_file(faces2.colwise() += vvv, "faces2");
	//}

	// Do Collision detection
	// INVERSE all below
	MatrixXi faces22;
	if (coll) {
		boxTriCollision(collisions, box->thresh, box->dim, box->E1, verts2, faces2);
		if (collisions.size() > 0) {
			for (int i = 0; i < collisions.size(); i++) {
				if (collisions[i]->count1 == 3 && collisions[i]->count2 == 1) {
					if (mesh.nodes[collisions[i]->verts2(0)]->has_coll_info) {
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
						collisions[i]->nor1(0) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0];
						collisions[i]->nor1(1) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1];
						collisions[i]->nor1(2) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2];
						collisions[i]->count1 = 2;
						collisions[i]->count2 = 2;
						//mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
						mesh.nodes[collisions[i]->verts2(0)]->was_on_pe = true;
						if (mesh.nodes[collisions[i]->verts2(0)]->preserve == true) mesh.nodes[collisions[i]->verts2(0)]->preserve = false;
						mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
						collision++;
						continue;
					}
					else {
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
						collision++;
					}
				}
				else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 1) {
					if (mesh.nodes[collisions[i]->verts2(0)]->has_coll_info) {
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
						collisions[i]->nor1(0) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0];
						collisions[i]->nor1(1) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1];
						collisions[i]->nor1(2) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2];
						if (mesh.nodes[collisions[i]->verts2(0)]->preserve == true) mesh.nodes[collisions[i]->verts2(0)]->preserve = false;
						mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
						collision++;
						continue;
					}
					else {
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
						collision++;
					}
				}
				else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 2) {
					if (mesh.nodes[collisions[i]->verts2(0)]->has_coll_info) {
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
						collisions[i]->nor1(0) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0];
						collisions[i]->nor1(1) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1];
						collisions[i]->nor1(2) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2];
						if (mesh.nodes[collisions[i]->verts2(0)]->preserve == true) mesh.nodes[collisions[i]->verts2(0)]->preserve = false;
						mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
						collision++;
						continue;
					}
					else {
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
						N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
						collision++;
					}
				}
				// These won't happen anymore on edge preserve
				else if (collisions[i]->count1 == 2 && collisions[i]->count2 == 2) {
					for (int j = 0; j < 2; j++) {
						N_.push_back(T(collision, collisions[i]->verts2(j) * 3, collisions[i]->nor1(0) * collisions[i]->weights2(j)));
						N_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 1, collisions[i]->nor1(1) * collisions[i]->weights2(j)));
						N_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 2, collisions[i]->nor1(2) * collisions[i]->weights2(j)));
					}
					//cout << "no 1" << endl;
					collision++;
				}
				else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 3) {
					for (int j = 0; j < 3; j++) {
						N_.push_back(T(collision, collisions[i]->verts2(j) * 3, collisions[i]->nor1(0) * collisions[i]->weights2(j)));
						N_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 1, collisions[i]->nor1(1) * collisions[i]->weights2(j)));
						N_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 2, collisions[i]->nor1(2) * collisions[i]->weights2(j)));
					}
					//cout << "no 2" << endl;
					collision++;
				}
			}
		}
	}

	if (collision == 0) {
		if (!fixon) {
			Eigen::SparseMatrix<double> M(sparse_size, sparse_size);
			Eigen::SparseMatrix<double> MDK(sparse_size, sparse_size);
			VectorXd b;
			VectorXd KKT_b;

			M.setFromTriplets(M_.begin(), M_.end());
			MDK.setFromTriplets(MDK_.begin(), MDK_.end());
			b = (M*v + h*f);

			ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
			cg.compute(MDK);
			v = cg.solve(b);
		}
		else {
			int holdcount = 0;
			for (int j = 0; j < fixed.rows(); j++) {
				if (fixed(j, 0) > 0) {
					if (j == 0) {
						//for (int k = 0; k < cols; k++) {
						//	MDK_.push_back(T(k * 3, sparse_size + (k * 3), fixed(j,0)));
						//	MDK_.push_back(T(k * 3 + 1, sparse_size + (k * 3) + 1, fixed(j, 1)));
						//	MDK_.push_back(T(k * 3 + 2, sparse_size + (k * 3) + 2, fixed(j, 2)));
						//	MDK_.push_back(T(sparse_size + (k * 3), k * 3, fixed(j, 0)));
						//	MDK_.push_back(T(sparse_size + (k * 3) + 1, k * 3 + 1, fixed(j, 1)));
						//	MDK_.push_back(T(sparse_size + (k * 3) + 2, k * 3 + 2, fixed(j, 2)));
						//	holdcount++;
						//}
					}
					else if (j == 4) {
						MDK_.push_back(T(0, sparse_size + (holdcount * 3), fixed(j, 0)));
						MDK_.push_back(T(1, sparse_size + (holdcount * 3) + 1, fixed(j, 1)));
						MDK_.push_back(T(2, sparse_size + (holdcount * 3) + 2, fixed(j, 2)));
						MDK_.push_back(T(sparse_size + (holdcount * 3), 0, fixed(j, 0)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 1, 1, fixed(j, 1)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 2, 2, fixed(j, 2)));
						holdcount++;
					}
					else if (j == 5) {
						MDK_.push_back(T((rows * (cols - 1)) * 3, sparse_size + (holdcount * 3), fixed(j, 0)));
						MDK_.push_back(T((rows * (cols - 1)) * 3 + 1, sparse_size + (holdcount * 3) + 1, fixed(j, 1)));
						MDK_.push_back(T((rows * (cols - 1)) * 3 + 2, sparse_size + (holdcount * 3) + 2, fixed(j, 2)));
						MDK_.push_back(T(sparse_size + (holdcount * 3), (rows * (cols - 1)) * 3, fixed(j, 0)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 1, (rows * (cols - 1)) * 3 + 1, fixed(j, 1)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 2, (rows * (cols - 1)) * 3 + 2, fixed(j, 2)));
						holdcount++;
					}
					else if (j == 6) {
						MDK_.push_back(T((rows * cols - 1) * 3, sparse_size + (holdcount * 3), fixed(j, 0)));
						MDK_.push_back(T((rows * cols - 1) * 3 + 1, sparse_size + (holdcount * 3) + 1, fixed(j, 1)));
						MDK_.push_back(T((rows * cols - 1) * 3 + 2, sparse_size + (holdcount * 3) + 2, fixed(j, 2)));
						MDK_.push_back(T(sparse_size + (holdcount * 3), (rows * cols - 1) * 3, fixed(j, 0)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 1, (rows * cols - 1) * 3 + 1, fixed(j, 1)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 2, (rows * cols - 1) * 3 + 2, fixed(j, 2)));
						holdcount++;
					}
					else if (j == 7) {
						MDK_.push_back(T((rows - 1) * 3, sparse_size + (holdcount * 3), fixed(j, 0)));
						MDK_.push_back(T((rows - 1) * 3 + 1, sparse_size + (holdcount * 3) + 1, fixed(j, 1)));
						MDK_.push_back(T((rows - 1) * 3 + 2, sparse_size + (holdcount * 3) + 2, fixed(j, 2)));
						MDK_.push_back(T(sparse_size + (holdcount * 3), (rows - 1) * 3, fixed(j, 0)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 1, (rows - 1) * 3 + 1, fixed(j, 1)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 2, (rows - 1) * 3 + 2, fixed(j, 2)));
						holdcount++;
					}
				}
			}

			Eigen::SparseMatrix<double> M(sparse_size, sparse_size);
			Eigen::SparseMatrix<double> MDK(sparse_size + (holdcount * 3), sparse_size + (holdcount * 3));
			VectorXd b;
			VectorXd KKT_b;

			M.setFromTriplets(M_.begin(), M_.end());
			MDK.setFromTriplets(MDK_.begin(), MDK_.end());
			b = (M*v + h*f);

			VectorXd append_zeros = VectorXd::Zero(holdcount * 3);
			int append_count = 0;
			for (int j = 0; j < fixed.rows(); j++) {
				if (fixed(j, 0) != -1) {
					if (j >= 4) {
						append_zeros(append_count * 3 + 0) = fixed(j, 3);
						append_zeros(append_count * 3 + 1) = fixed(j, 4);
						append_zeros(append_count * 3 + 2) = fixed(j, 5);
						append_count++;
					}
				}
			}

			KKT_b.resize(b.size() + (holdcount * 3));
			KKT_b << b, append_zeros;

			LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
			lscg.compute(MDK);
			v = lscg.solve(KKT_b);
		}
	}
	else {
		//cout << "Quadratic Solver" << endl;
		Eigen::SparseMatrix<double> M(sparse_size, sparse_size);
		Eigen::SparseMatrix<double> MDK(sparse_size, sparse_size);
		VectorXd b;

		M.setFromTriplets(M_.begin(), M_.end());
		MDK.setFromTriplets(MDK_.begin(), MDK_.end());
		b = -(M*v + h*f);

		Eigen::SparseMatrix<double> N(collision, sparse_size);
		N.setFromTriplets(N_.begin(), N_.end());
		igl::mosek::MosekData mosek_data;
		auto lc = Eigen::VectorXd(collision);
		auto uc = Eigen::VectorXd(collision);
		auto lx = Eigen::VectorXd(sparse_size);
		auto ux = Eigen::VectorXd(sparse_size);
		for (int i = 0; i < sparse_size; i++) {
			lx(i) = -MSK_INFINITY;
			ux(i) = +MSK_INFINITY;
		}
		// TODO:: Better way
		int append_count = 0;
		for (int j = 0; j < fixed.rows(); j++) {
			if (fixed(j, 0) != -1) {
				if (j == 4) {
					lx(0) = fixed(j, 3);
					ux(0) = fixed(j, 3);
					lx(1) = fixed(j, 4);
					ux(1) = fixed(j, 4);
					lx(2) = fixed(j, 5);
					ux(2) = fixed(j, 5);
				}
				else if (j == 5) {
					lx((rows * (cols - 1)) * 3) = fixed(j, 3);
					ux((rows * (cols - 1)) * 3) = fixed(j, 3);
					lx((rows * (cols - 1)) * 3 + 1) = fixed(j, 4);
					ux((rows * (cols - 1)) * 3 + 1) = fixed(j, 4);
					lx((rows * (cols - 1)) * 3 + 2) = fixed(j, 5);
					ux((rows * (cols - 1)) * 3 + 2) = fixed(j, 5);
				}
				else if (j == 6) {
					lx((rows * cols - 1) * 3) = fixed(j, 3);
					ux((rows * cols - 1) * 3) = fixed(j, 3);
					lx((rows * cols - 1) * 3 + 1) = fixed(j, 4);
					ux((rows * cols - 1) * 3 + 1) = fixed(j, 4);
					lx((rows * cols - 1) * 3 + 2) = fixed(j, 5);
					ux((rows * cols - 1) * 3 + 2) = fixed(j, 5);
				}
				else if (j == 7) {
					lx((rows - 1) * 3) = fixed(j, 3);
					ux((rows - 1) * 3) = fixed(j, 3);
					lx((rows - 1) * 3 + 1) = fixed(j, 4);
					ux((rows - 1) * 3 + 1) = fixed(j, 4);
					lx((rows - 1) * 3 + 2) = fixed(j, 5);
					ux((rows - 1) * 3 + 2) = fixed(j, 5);
				}
			}
		}
		for (int i = 0; i < collision; i++) {
			//lc(i) = -(0.1/h)* 0.001; // BAUMGARTE
			lc(i) = 0.0;
			uc(i) = +MSK_INFINITY;
		}

		if (matlab_debug_physics) {
			sparse_to_file_as_dense(MDK, "MDK");
			vec_to_file(b, "b");
			sparse_to_file_as_dense(N, "N");
			sparse_to_file_as_dense(M, "M");
			vec_to_file(lx, "lx");
			vec_to_file(ux, "ux");
		}

		//dfile.close();
		//SparseMatrix<double> spMat;
		//MatrixXd dMat;
		//dMat = MatrixXd(N);
		//cout << dMat << endl;
		//cout << "MDK: " << MDK.rows() << " x " << MDK.cols() << endl;
		//cout << "KKT_b: " << b.rows() << " x " << b.cols() << endl;
		//cout << "N: " << N << " x " << N.cols() << endl;
		//cout << "lc: " << lc.rows() << " x " << lc.cols() << endl;
		//cout << "lx: " << lx.rows() << " x " << lx.cols() << endl;
		//cout << "v: " << v.rows() << " x " << v.cols() << endl;
		////cout << "STEP" << endl;
		////cout << v << endl;
		//cout << N << endl;
		//cout << MDK.triangularView<Lower>() << endl;
		/*MatrixXd ttt;
		ttt = MatrixXd(MDK);
		cout << MDK << endl;
		Eigen::EigenSolver<MatrixXd> es(ttt);
		cout << "The eigenvalues of MDK are:" << endl << es.eigenvalues() << endl;*/
		igl::mosek::mosek_quadprog(MDK, b, 0, N, lc, uc, lx, ux, mosek_data, v);
		//v = -v; // INVERSE

		//if (matlab_debug_physics) {
		//	sparse_to_file_as_dense(MDK, "MDK");
		//	sparse_to_file_as_dense(b, "b");
		//	sparse_to_file_as_dense(N, "N");
		//	sparse_to_file_as_dense(lc, "lc");
		//	sparse_to_file_as_dense(uc, "uc");
		//	sparse_to_file_as_dense(lx, "lx");
		//	sparse_to_file_as_dense(ux, "ux");
		//	sparse_to_file_as_dense(v, "v");
		//}
	}

	for (int i = 0; i < mesh.nodes.size(); i++) {
		mesh.nodes[i]->v[0] = v(i * 3);
		mesh.nodes[i]->v[1] = v(i * 3 + 1);
		mesh.nodes[i]->v[2] = v(i * 3 + 2);
		mesh.nodes[i]->x = mesh.nodes[i]->x + h * mesh.nodes[i]->v;
	}
}