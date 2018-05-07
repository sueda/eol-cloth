#include <iostream>
#include <string>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <fstream>

#include "Cloth.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include "Box.h"
#include "boxTriCollision.h"
#include "MatlabDebug.h"
#include "ChronoTimer.h"
#include "Rigid.h"

// FEM calculations
#include "ComputeInertialLAG.h"
#include "ComputeMembraneLAG.h"
#include "ComputeBendingLAG.h"
#include "ComputeInertiaEOL.h"
#include "ComputeMembraneEOL.h"
#include "ComputeBendingEOL.h"

// Mosek-libigl
//#include "igl\mosek\mosek_quadprog.h"
#include "mosek.h"
#include "QuadProgMosek.h"
//#include "igl\active_set.h"

// Triangle-libigl
//#include "igl\triangle\cdt.h"
//#include "triangle_api.h"
//#include "triangle.h"
//#include "triangle_internal.h"
#include "Fade_2D.h"

// ArcSim
#include "mesh.hpp"
#include "vectors.hpp"
#include "util.hpp"
#include "io.hpp"
#include "geometry.hpp"
#include "remesh.hpp"

//extern "C" {
////#include "triangle.h"
//#include "triangle.c"
//}

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;
vector<T> N_;

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

Cloth::Cloth() :
	tmpbool(false),
	safety_margin(0.01),
	energyadd(1.0)
{

}

Cloth::~Cloth()
{
}

void Cloth::setup(const Vector3d &x00,
	const Vector3d &x01,
	const Vector3d &x10,
	const Vector3d &x11)
{
	Xmin = x00(0);
	Xmax = x11(0);
	Ymin = x00(1);
	Ymax = x11(1);

	// Create main grid particles
	int nVerts = rows*cols;
	for (int i = 0; i < rows; ++i) {
		double u = i / (rows - 1.0);
		Vector3d x0 = (1 - u)*x00 + u*x10;
		Vector3d x1 = (1 - u)*x01 + u*x11;
		for (int j = 0; j < cols; ++j) {
			double v = j / (cols - 1.0);
			Vector3d x = (1 - v)*x0 + v*x1;

			Vec3 ver(x[0], x[1], 0);
			Vec3 vel(0.0, 0.0, 0.0);
			//auto u = make_shared<Node>(n, Vec3(0));
			mesh.add(new Vert(ver, vel)); // TODO:: Shared pointer

		}
	}
	double dx = ((x01(0) - x00(0)) / (cols - 1)) / 2;
	double dy = ((x10(1) - x00(1)) / (rows - 1)) / 2;
	// Create grid center particles
	for (int i = 0; i < rows - 1; ++i) {
		for (int j = 0; j < cols - 1; ++j) {
			Vector3d x = x00;
			x(0) += (2 * j + 1) * dx;
			x(1) += (2 * i + 1) * dy;

			Vec3 ver(x[0], x[1], 0);
			Vec3 vel(0.0, 0.0, 0.0);
			//auto u = make_shared<Node>(n, Vec3(0));
			mesh.add(new Vert(ver, vel)); // TODO:: Shared pointer
		}
	}

	// Add nodes after texture verts
	for (int i = 0; i < rows; ++i) {
		double u = i / (rows - 1.0);
		Vector3d x0 = (1 - u)*x00 + u*x10;
		Vector3d x1 = (1 - u)*x01 + u*x11;
		for (int j = 0; j < cols; ++j) {
			double v = j / (cols - 1.0);
			Vector3d x = (1 - v)*x0 + v*x1;

			Vec3 n(x[0], x[1], x[2]);
			//auto u = make_shared<Node>(n, Vec3(0));
			mesh.add(new Node(n, n, Vec3(0), 0, 0, false)); // TODO:: Shared pointer

		}
	}
	for (int i = 0; i < rows - 1; ++i) {
		for (int j = 0; j < cols - 1; ++j) {
			Vector3d x = x00;
			x(0) += (2 * j + 1) * dx;
			x(1) += (2 * i + 1) * dy;

			Vec3 n(x[0], x[1], x[2]);
			//auto u = make_shared<Node>(n, Vec3(0));
			mesh.add(new Node(n, n, Vec3(0), 0, 0, false)); // TODO:: Shared pointer
		}
	}

	// Create FEs
	// k0 --  k0 + 1
	// |  \       /   |
	// |     kc0       |
	// |  /       \   |
	// k0 + cols   ----    k0 + cols + 1
	int center_index_cnt = 0;
	for (int i = 0; i < rows - 1; i++) {
		for (int j = 0; j < cols - 1; j++) {
			int k0 = (i * cols) + j; // upper right index
			int kc0 = ((rows * cols) + center_index_cnt); // center index
			center_index_cnt++;

			// Fill single mesh triangle for nodes and verts
			// TODO:: Enclose in functions
			// TODO:: Check fill order
			vector<Vert*> verts1;
			vector<Node*> nodes1;
			nodes1.push_back(mesh.nodes[k0]);
			verts1.push_back(mesh.verts[k0]);
			nodes1.push_back(mesh.nodes[k0 + 1]);
			verts1.push_back(mesh.verts[k0 + 1]);
			nodes1.push_back(mesh.nodes[kc0]);
			verts1.push_back(mesh.verts[kc0]);
			for (int v = 0; v < verts1.size(); v++)
				connect(verts1[v], nodes1[v]);
			vector<Face*> faces1 = triangulateARC(verts1);
			for (int f = 0; f < faces1.size(); f++)
				mesh.add(faces1[f]);

			vector<Vert*> verts2;
			vector<Node*> nodes2;
			nodes2.push_back(mesh.nodes[k0 + 1]);
			verts2.push_back(mesh.verts[k0 + 1]);
			nodes2.push_back(mesh.nodes[k0 + cols + 1]);
			verts2.push_back(mesh.verts[k0 + cols + 1]);
			nodes2.push_back(mesh.nodes[kc0]);
			verts2.push_back(mesh.verts[kc0]);
			for (int v = 0; v < verts2.size(); v++)
				connect(verts2[v], nodes2[v]);
			vector<Face*> faces2 = triangulateARC(verts2);
			for (int f = 0; f < faces2.size(); f++)
				mesh.add(faces2[f]);

			vector<Vert*> verts3;
			vector<Node*> nodes3;
			nodes3.push_back(mesh.nodes[k0 + cols + 1]);
			verts3.push_back(mesh.verts[k0 + cols + 1]);
			nodes3.push_back(mesh.nodes[k0 + cols]);
			verts3.push_back(mesh.verts[k0 + cols]);
			nodes3.push_back(mesh.nodes[kc0]);
			verts3.push_back(mesh.verts[kc0]);
			for (int v = 0; v < verts3.size(); v++)
				connect(verts3[v], nodes3[v]);
			vector<Face*> faces3 = triangulateARC(verts3);
			for (int f = 0; f < faces3.size(); f++)
				mesh.add(faces3[f]);

			vector<Vert*> verts4;
			vector<Node*> nodes4;
			nodes4.push_back(mesh.nodes[k0 + cols]);
			verts4.push_back(mesh.verts[k0 + cols]);
			nodes4.push_back(mesh.nodes[k0]);
			verts4.push_back(mesh.verts[k0]);
			nodes4.push_back(mesh.nodes[kc0]);
			verts4.push_back(mesh.verts[kc0]);
			for (int v = 0; v < verts4.size(); v++)
				connect(verts4[v], nodes4[v]);
			vector<Face*> faces4 = triangulateARC(verts4);
			for (int f = 0; f < faces4.size(); f++)
				mesh.add(faces4[f]);

			// Counter clockwise ordering but bugged
			/*vector<Vert*> verts1;
			vector<Node*> nodes1;
			nodes1.push_back(mesh.nodes[k0]);
			verts1.push_back(mesh.verts[k0]);
			nodes1.push_back(mesh.nodes[k0 + cols]);
			verts1.push_back(mesh.verts[k0 + cols]);
			nodes1.push_back(mesh.nodes[kc0]);
			verts1.push_back(mesh.verts[kc0]);
			for (int v = 0; v < verts1.size(); v++)
				connect(verts1[v], nodes1[v]);
			vector<Face*> faces1 = triangulate(verts1);
			for (int f = 0; f < faces1.size(); f++)
				mesh.add(faces1[f]);

			vector<Vert*> verts2;
			vector<Node*> nodes2;
			nodes2.push_back(mesh.nodes[k0 + cols]);
			verts2.push_back(mesh.verts[k0 + cols]);
			nodes2.push_back(mesh.nodes[k0 + cols + 1]);
			verts2.push_back(mesh.verts[k0 + cols + 1]);
			nodes2.push_back(mesh.nodes[kc0]);
			verts2.push_back(mesh.verts[kc0]);
			for (int v = 0; v < verts2.size(); v++)
				connect(verts2[v], nodes2[v]);
			vector<Face*> faces2 = triangulate(verts2);
			for (int f = 0; f < faces2.size(); f++)
				mesh.add(faces2[f]);

			vector<Vert*> verts3;
			vector<Node*> nodes3;
			nodes3.push_back(mesh.nodes[k0 + cols + 1]);
			verts3.push_back(mesh.verts[k0 + cols + 1]);
			nodes3.push_back(mesh.nodes[k0 + 1]);
			verts3.push_back(mesh.verts[k0 + 1]);
			nodes3.push_back(mesh.nodes[kc0]);
			verts3.push_back(mesh.verts[kc0]);
			for (int v = 0; v < verts3.size(); v++)
				connect(verts3[v], nodes3[v]);
			vector<Face*> faces3 = triangulate(verts3);
			for (int f = 0; f < faces3.size(); f++)
				mesh.add(faces3[f]);

			vector<Vert*> verts4;
			vector<Node*> nodes4;
			nodes4.push_back(mesh.nodes[k0 + 1]);
			verts4.push_back(mesh.verts[k0 + 1]);
			nodes4.push_back(mesh.nodes[k0]);
			verts4.push_back(mesh.verts[k0]);
			nodes4.push_back(mesh.nodes[kc0]);
			verts4.push_back(mesh.verts[kc0]);
			for (int v = 0; v < verts4.size(); v++)
				connect(verts4[v], nodes4[v]);
			vector<Face*> faces4 = triangulate(verts4);
			for (int f = 0; f < faces4.size(); f++)
				mesh.add(faces4[f]);*/
		}
	}

	for (int i = 0; i < mesh.faces.size(); i++) {
		mesh.faces[i]->material = &material;
	}

	// Further mesh setup
	mark_nodes_to_preserve(mesh);
	compute_ms_data(mesh);

	v.resize(mesh.nodes.size() * 3);
	f.resize(mesh.nodes.size() * 3);
	EOLverts.resize(mesh.nodes.size());
	EOLverts.setZero();

	// Build vertex buffers
	posBuf.clear();
	posBuf2D.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(mesh.nodes.size() * 3);
	posBuf2D.resize(mesh.nodes.size() * 3);
	norBuf.resize(mesh.nodes.size() * 3);
	texBuf.resize(mesh.nodes.size() * 2);
	eleBuf.resize(mesh.faces.size() * 3);

	updatePosNor();

	ClothTimer.push_back(make_shared<ChronoTimer>("CD2"));
	ClothTimer.push_back(make_shared<ChronoTimer>("VT"));
	ClothTimer.push_back(make_shared<ChronoTimer>("MatrixFill"));
	ClothTimer.push_back(make_shared<ChronoTimer>("VelocityIntegration"));
}

void Cloth::setup(const string &filename)
{
	load_cloth_obj(mesh, filename);

	for (int i = 0; i < mesh.faces.size(); i++) {
		mesh.faces[i]->material = &material;
	}

	// Further mesh setup
	mark_nodes_to_preserve(mesh);
	compute_ms_data(mesh);

	v_old.resize(mesh.nodes.size() * 3);
	v.resize(mesh.nodes.size() * 3);
	f.resize(mesh.nodes.size() * 3);

	// Build vertex buffers
	posBuf.clear();
	posBuf2D.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(mesh.nodes.size() * 3);
	posBuf2D.resize(mesh.nodes.size() * 3);
	norBuf.resize(mesh.nodes.size() * 3);
	texBuf.resize(mesh.nodes.size() * 2);
	eleBuf.resize(mesh.faces.size() * 3);

	updatePosNor();

	ClothTimer.push_back(make_shared<ChronoTimer>("CD2"));
	ClothTimer.push_back(make_shared<ChronoTimer>("VT"));
	ClothTimer.push_back(make_shared<ChronoTimer>("MatrixFill"));
	ClothTimer.push_back(make_shared<ChronoTimer>("VelocityIntegration"));
}

void Cloth::tare()
{

}

void Cloth::reset()
{
	//updatePosNor();
}

void Cloth::updatePosNor()
{	
	for (int i = 0; i < mesh.nodes.size(); i++) {
		Vec3 xm = mesh.nodes[i]->x;
		Vec3 Xm = mesh.nodes[i]->verts[0]->u;
		Vec3 nm = mesh.nodes[i]->n;
		Vec3 tm = mesh.nodes[i]->verts[0]->u;
		posBuf[3 * i + 0] = xm[0];
		posBuf[3 * i + 1] = xm[1];
		posBuf[3 * i + 2] = xm[2];
		posBuf2D[3 * i + 0] = Xm[0];
		posBuf2D[3 * i + 1] = Xm[1];
		posBuf2D[3 * i + 2] = Xm[2];
		norBuf[3 * i + 0] = nm[0];
		norBuf[3 * i + 1] = nm[1];
		norBuf[3 * i + 2] = nm[2];
		texBuf[2 * i + 0] = tm[0];
		texBuf[2 * i + 1] = tm[1];
	}
	for (int i = 0; i < mesh.faces.size(); i++) {
		eleBuf[3 * i + 0] = mesh.faces[i]->v[0]->node->index;
		eleBuf[3 * i + 1] = mesh.faces[i]->v[1]->node->index;
		eleBuf[3 * i + 2] = mesh.faces[i]->v[2]->node->index;
	}
	//cout << endl;
}

void Cloth::updateBuffers() 
{
	posBuf.clear();
	posBuf2D.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(mesh.nodes.size() * 3);
	posBuf2D.resize(mesh.nodes.size() * 3);
	norBuf.resize(mesh.nodes.size() * 3);
	texBuf.resize(mesh.nodes.size() * 2);
	eleBuf.resize(mesh.faces.size() * 3);

	// Update position and normal buffers
	updatePosNor();
}

// Lagrangian only physics step
void Cloth::stepL(double h, const Vector3d &grav, vector<shared_ptr<Box> > box, double t)
{
	v.resize(mesh.nodes.size() * 3);
	f.resize(mesh.nodes.size() * 3);

	v.setZero();
	f.setZero();

	int collision = 0; // TODO:: Move

	bool fixed_points = false;  // TODO:: Restructure

	vector<T> M_;
	vector<T> MDK_;
	vector<T> Aineq_;
	vector<T> Aeq_;
	//vector<T> N_;
	N_.clear();
	int sparse_size = 0;

	vector<int> ineq_move_index;
	vector<double> ineq_move;
	vector<int> eq_move_index;
	vector<double> eq_move;

	// Box Setup
	//collisions.clear();
	MatrixXd verts2(3,mesh.nodes.size());
	MatrixXi faces2(3,mesh.faces.size());

	MatrixXd x_X(mesh.nodes.size(), 5);
	VectorXi isEOL(mesh.nodes.size());

	ClothTimer[2]->tic();
	for (int i = 0; i < mesh.nodes.size(); i++) {

		isEOL(i) = 0;
		x_X(i, 0) = mesh.nodes[i]->x[0];
		x_X(i, 1) = mesh.nodes[i]->x[1];
		x_X(i, 2) = mesh.nodes[i]->x[2];
		x_X(i, 3) = mesh.nodes[i]->verts[0]->u[0];
		x_X(i, 4) = mesh.nodes[i]->verts[0]->u[1];

		// For collision
		if(coll) verts2.col(i) = Vector3d(mesh.nodes[i]->x[0],mesh.nodes[i]->x[1],mesh.nodes[i]->x[2]);

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
		if(coll) faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);

		// Set up corrdaintes
		double xa[3];
		double xb[3];
		double xc[3];
		double Xa[2];
		double Xb[2];
		double Xc[2];

		// Fill temporary vectors 
		// TODO:: Convrsion functions
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
		ComputeInertiaLAG(xa, xb, xc, Xa, Xb, Xc, g, material.density, Wi, fi, Mi); 
		VectorXd fie = Map<VectorXd>(fi,9);
		MatrixXd Mie = Map<MatrixXd>(Mi, 9,9);



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
		VectorXd fme = Map<VectorXd>(fm,9);
		MatrixXd Kme = Map<MatrixXd>(Km,9,9); 

		Matrix3d Maa = Mie.block(0, 0, 3, 3);
		Matrix3d Kaa = Kme.block(0, 0, 3, 3);
		Matrix3d MDKaa = Maa + damping(0) * h * Maa + damping(1) * h * h * Kaa;
		f.segment<3>(mesh.faces[i]->v[0]->node->index*3) += (fie.segment<3>(0) + fme.segment<3>(0));
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
		VectorXd fbe = Map<VectorXd>(fb, 12);
		MatrixXd Kbe = Map<MatrixXd>(Kb, 12, 12);

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
	ClothTimer[2]->toc();

	// Do Collision detection
	ClothTimer[0]->tic();
	VectorXi CisEOL;
	if (coll) {
		for (int b = 0; b < box.size(); b++) {
			collisions.clear();
			if(wire) boxTriCollisionHack2(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, CisEOL, false);
			else boxTriCollision(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, CisEOL, false);
			//if (wire && move_wire) {
			//	box[b]->E1(2, 3) += 1.0 * 2.5 * h;
			//	box[b]->x(2) += 1.0 * 2.5 * h;
			//}
			if (collisions.size() > 0) {
				for (int i = 0; i < collisions.size(); i++) {
					if (collisions[i]->count1 == 3 && collisions[i]->count2 == 1) {
						if (mesh.nodes[collisions[i]->verts2(0)]->has_coll_info) {
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
							//mesh.nodes[collisions[i]->verts2(0)]->coll_norm = Vec3(0.0, 0.0, 1.0);
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
							collision++;
							collisions[i]->nor1(0) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0];
							collisions[i]->nor1(1) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1];
							collisions[i]->nor1(2) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2];
							collisions[i]->count1 = 2;
							collisions[i]->count2 = 2;
							//mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
							if (mesh.nodes[collisions[i]->verts2(0)]->preserve == true) mesh.nodes[collisions[i]->verts2(0)]->preserve = false;
							if (mesh.nodes[collisions[i]->verts2(0)]->on_corner) mesh.nodes[collisions[i]->verts2(0)]->on_corner = false;
							mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
							collision++;
							continue;
						}
						else {
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3, -collisions[i]->nor1(0)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, -collisions[i]->nor1(1)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, -collisions[i]->nor1(2)));
							collision++;
							if (move_wire) {
								//collisions[i]->nor2.normalized() * box[b]->v.segment<3>(3).dot(collisions[i]->nor2.normalized());
								ineq_move_index.push_back(collision - 1);
								ineq_move.push_back(-box[b]->v.segment<3>(3).dot(collisions[i]->nor2.normalized()));
							}
						}
					}
					else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 1) {
						if (mesh.nodes[collisions[i]->verts2(0)]->has_coll_info) {
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
							collisions[i]->nor1(0) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0];
							collisions[i]->nor1(1) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1];
							collisions[i]->nor1(2) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2];
							if (mesh.nodes[collisions[i]->verts2(0)]->preserve == true) mesh.nodes[collisions[i]->verts2(0)]->preserve = false;
							if (mesh.nodes[collisions[i]->verts2(0)]->on_corner) mesh.nodes[collisions[i]->verts2(0)]->on_corner = false;
							mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
							collision++;
							continue;
						}
						else {
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3, -collisions[i]->nor1(0)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, -collisions[i]->nor1(1)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, -collisions[i]->nor1(2)));
							collision++;
						}
					}
					else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 2) {
						if (mesh.nodes[collisions[i]->verts2(0)]->has_coll_info) {
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, -mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
							collisions[i]->nor1(0) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0];
							collisions[i]->nor1(1) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1];
							collisions[i]->nor1(2) = mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2];
							if (mesh.nodes[collisions[i]->verts2(0)]->preserve == true) mesh.nodes[collisions[i]->verts2(0)]->preserve = false;
							if (mesh.nodes[collisions[i]->verts2(0)]->on_corner) mesh.nodes[collisions[i]->verts2(0)]->on_corner = false;
							mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
							collision++;
							continue;
						}
						else {
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
							//N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3, -collisions[i]->nor1(0)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, -collisions[i]->nor1(1)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, -collisions[i]->nor1(2)));
							collision++;
						}
					}
					// These won't happen anymore on edge preserve
					else if (collisions[i]->count1 == 2 && collisions[i]->count2 == 2) {
						if (wire) {
							Vec3 clothedge = normalize(mesh.nodes[collisions[i]->verts2(0)]->x - mesh.nodes[collisions[i]->verts2(1)]->x);
							Vec3 newconstraint = normalize(cross(clothedge,Vec3(1.0,0.0,0.0)));
							if (newconstraint[2] < 0.0) newconstraint *= -1.0;
							collisions[i]->nor2 = Vector3d(newconstraint[0], newconstraint[1], newconstraint[2]);
							collisions[i]->nor1 = Vector3d(newconstraint[0], newconstraint[1], newconstraint[2]);
						}
						for (int j = 0; j < 2; j++) {
							//N_.push_back(T(collision, collisions[i]->verts2(j) * 3, collisions[i]->nor1(0) * collisions[i]->weights2(j)));
							//N_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 1, collisions[i]->nor1(1) * collisions[i]->weights2(j)));
							//N_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 2, collisions[i]->nor1(2) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(j) * 3, -collisions[i]->nor2(0) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 1, -collisions[i]->nor2(1) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 2, -collisions[i]->nor2(2) * collisions[i]->weights2(j)));
						}
						//cout << "no 1" << endl;
						collision++;
						if (move_wire) {
							//collisions[i]->nor2.normalized() * box[b]->v.segment<3>(3).dot(collisions[i]->nor2.normalized());
							ineq_move_index.push_back(collision - 1);
							ineq_move.push_back(-box[b]->v.segment<3>(3).dot(collisions[i]->nor2.normalized()));
						}
					}
					else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 3) {
						for (int j = 0; j < 3; j++) {
							//N_.push_back(T(collision, collisions[i]->verts2(j) * 3, collisions[i]->nor1(0) * collisions[i]->weights2(j)));
							//N_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 1, collisions[i]->nor1(1) * collisions[i]->weights2(j)));
							//N_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 2, collisions[i]->nor1(2) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(j) * 3, -collisions[i]->nor2(0) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 1, -collisions[i]->nor2(1) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision, collisions[i]->verts2(j) * 3 + 2, -collisions[i]->nor2(2) * collisions[i]->weights2(j)));
						}
						//cout << "no 2" << endl;
						collision++;
					}
				}
			}
		}
	}

	ClothTimer[0]->toc();

	ClothTimer[3]->tic();
	if (collision == 0 && false) {
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
				if (fixed(j,0) > 0) {
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
						MDK_.push_back(T((cols * (rows - 1)) * 3, sparse_size + (holdcount * 3), fixed(j, 0)));
						MDK_.push_back(T((cols * (rows - 1)) * 3 + 1, sparse_size + (holdcount * 3) + 1, fixed(j, 1)));
						MDK_.push_back(T((cols * (rows - 1)) * 3 + 2, sparse_size + (holdcount * 3) + 2, fixed(j, 2)));
						MDK_.push_back(T(sparse_size + (holdcount * 3), (cols * (rows - 1)) * 3, fixed(j, 0)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 1, (cols * (rows - 1)) * 3 + 1, fixed(j, 1)));
						MDK_.push_back(T(sparse_size + (holdcount * 3) + 2, (cols * (rows - 1)) * 3 + 2, fixed(j, 2)));
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
			if (step > pull_step) {
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

		int aeqsize = 0;

		//Eigen::SparseMatrix<double> N(collision, sparse_size);
		//N.setFromTriplets(N_.begin(), N_.end());
		////igl::mosek::MosekData mosek_data;
		//auto lc = Eigen::VectorXd(collision);
		//auto uc = Eigen::VectorXd(collision);
		//auto lx = Eigen::VectorXd(sparse_size);
		//auto ux = Eigen::VectorXd(sparse_size);
		//for (int i = 0; i < sparse_size; i++) {
		//	lx(i) = -MSK_INFINITY;
		//	ux(i) = +MSK_INFINITY;
		//}
		// TODO:: Better way
		int append_count = 0;
		vector<double> beq_;
		double expofilA = 0.01;
		//for (int j = 0; j < fixed.rows(); j++) {
		//	if (fixed(j, 0) != -1) {
		//		if (j == 4) {
		//			//lx(0) = fixed(j, 3);
		//			//ux(0) = fixed(j, 3);
		//			//lx(1) = fixed(j, 4);
		//			//ux(1) = fixed(j, 4);
		//			//lx(2) = fixed(j, 5);
		//			//ux(2) = fixed(j, 5);
		//			if (step > pull_step) {
		//				beq_.push_back((1 - expofilA) * mesh.nodes[0]->v[0] + expofilA * fixed(j, 3));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[0]->v[1] + expofilA * fixed(j, 4));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[0]->v[2] + expofilA * fixed(j, 5));
		//			}
		//			Aeq_.push_back(T(aeqsize, 0, 1));
		//			aeqsize++;
		//			Aeq_.push_back(T(aeqsize, 1, 1));
		//			aeqsize++;
		//			Aeq_.push_back(T(aeqsize, 2, 1));
		//			aeqsize++;
		//		}
		//		else if (j == 5) {
		//			//lx((rows * (cols - 1)) * 3) = fixed(j, 3);
		//			//ux((rows * (cols - 1)) * 3) = fixed(j, 3);
		//			//lx((rows * (cols - 1)) * 3 + 1) = fixed(j, 4);
		//			//ux((rows * (cols - 1)) * 3 + 1) = fixed(j, 4);
		//			//lx((rows * (cols - 1)) * 3 + 2) = fixed(j, 5);
		//			//ux((rows * (cols - 1)) * 3 + 2) = fixed(j, 5);
		//			if (step > pull_step) {
		//				beq_.push_back((1 - expofilA) * mesh.nodes[cols * (rows - 1)]->v[0] + expofilA * fixed(j, 3));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[cols * (rows - 1)]->v[1] + expofilA * fixed(j, 4));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[cols * (rows - 1)]->v[2] + expofilA * fixed(j, 5));
		//			}
		//			Aeq_.push_back(T(aeqsize, (cols * (rows - 1)) * 3, 1));
		//			aeqsize++;
		//			Aeq_.push_back(T(aeqsize, (cols * (rows - 1)) * 3 + 1, 1));
		//			aeqsize++;
		//			Aeq_.push_back(T(aeqsize, (cols * (rows - 1)) * 3 + 2, 1));
		//			aeqsize++;
		//		}
		//		else if (j == 6) {
		//			//lx((rows * cols - 1) * 3) = fixed(j, 3);
		//			//ux((rows * cols - 1) * 3) = fixed(j, 3);
		//			//lx((rows * cols - 1) * 3 + 1) = fixed(j, 4);
		//			//ux((rows * cols - 1) * 3 + 1) = fixed(j, 4);
		//			//lx((rows * cols - 1) * 3 + 2) = fixed(j, 5);
		//			//ux((rows * cols - 1) * 3 + 2) = fixed(j, 5);
		//			if (step > pull_step) {
		//				beq_.push_back((1 - expofilA) * mesh.nodes[rows * cols - 1]->v[0] + expofilA * fixed(j, 3));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[rows * cols - 1]->v[1] + expofilA * fixed(j, 4));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[rows * cols - 1]->v[2] + expofilA * fixed(j, 5));
		//			}
		//			Aeq_.push_back(T(aeqsize, (rows * cols - 1) * 3, 1));
		//			aeqsize++;
		//			Aeq_.push_back(T(aeqsize, (rows * cols - 1) * 3 + 1, 1));
		//			aeqsize++;
		//			Aeq_.push_back(T(aeqsize, (rows * cols - 1) * 3 + 2, 1));
		//			aeqsize++;
		//		}
		//		else if (j == 7) {
		//			//lx((rows - 1) * 3) = fixed(j, 3);
		//			//ux((rows - 1) * 3) = fixed(j, 3);
		//			//lx((rows - 1) * 3 + 1) = fixed(j, 4);
		//			//ux((rows - 1) * 3 + 1) = fixed(j, 4);
		//			//lx((rows - 1) * 3 + 2) = fixed(j, 5);
		//			//ux((rows - 1) * 3 + 2) = fixed(j, 5);
		//			if (step > pull_step) {
		//				beq_.push_back((1 - expofilA) * mesh.nodes[rows - 1]->v[0] + expofilA * fixed(j, 3));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[rows - 1]->v[1] + expofilA * fixed(j, 4));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[rows - 1]->v[2] + expofilA * fixed(j, 5));
		//			}
		//			Aeq_.push_back(T(aeqsize, (rows - 1) * 3, 1));
		//			aeqsize++;
		//			Aeq_.push_back(T(aeqsize, (rows - 1) * 3 + 1, 1));
		//			aeqsize++;
		//			Aeq_.push_back(T(aeqsize, (rows - 1) * 3 + 2, 1));
		//			aeqsize++;
		//		}
		//	}
		//}
		//for (int i = 0; i < collision; i++) {
		//	//lc(i) = -(0.1/h)* 0.001; // BAUMGARTE
		//	lc(i) = 0.0;
		//	uc(i) = +MSK_INFINITY;
		//}

		// Pull sides
		for (int i = 0; i < mesh.nodes.size(); i++) {
			if (is_seam_or_boundary(mesh.nodes[i])) {
				//if (mesh.nodes[i]->preserve) continue;
				/*if (mesh.nodes[i]->verts[0]->u[0] == 0.0) {
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3, 1));
				aeqsize++;
				append_count++;
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 1, 1));
				aeqsize++;
				append_count++;
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 2, 1));
				aeqsize++;
				append_count++;
				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * fixed(4, 3));
				if (step < side_pull_change_step) {
				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * 0);
				}
				else {
				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * fixed(4,4));
				}
				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[2] + expofilA * fixed(4, 5));
				}
				if (mesh.nodes[i]->verts[0]->u[0] == 1.0) {
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3, 1));
				aeqsize++;
				append_count++;
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 1, 1));
				aeqsize++;
				append_count++;
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 2, 1));
				aeqsize++;
				append_count++;
				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * fixed(6, 3));
				if (step < side_pull_change_step) {
				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * 0);
				}
				else {
				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * fixed(4, 4));
				}
				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[2] + expofilA * fixed(4, 5));
				}*/
				if (mesh.nodes[i]->verts[0]->u[1] == 0.0) {
					Aeq_.push_back(T(aeqsize, i * 3, 1));
					aeqsize++;
					append_count++;
					Aeq_.push_back(T(aeqsize, i * 3 + 1, 1));
					aeqsize++;
					append_count++;
					Aeq_.push_back(T(aeqsize, i * 3 + 2, 1));
					aeqsize++;
					append_count++;
					if (step < side_pull_change_step) {
						beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * 0);
					}
					else {
						beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * fixed(4, 3));
					}
					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * fixed(4, 4));
					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[2] + expofilA * fixed(4, 5));
				}
				if (mesh.nodes[i]->verts[0]->u[1] == 1.0) {
					Aeq_.push_back(T(aeqsize, i * 3, 1));
					aeqsize++;
					append_count++;
					Aeq_.push_back(T(aeqsize, i * 3 + 1, 1));
					aeqsize++;
					append_count++;
					Aeq_.push_back(T(aeqsize, i * 3 + 2, 1));
					aeqsize++;
					append_count++;
					if (step < side_pull_change_step) {
						beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * 0);
					}
					else {
						beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * fixed(5, 3));
					}
					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * fixed(5, 4));
					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[2] + expofilA * fixed(5, 5));
				}
			}
		}

		Eigen::SparseMatrix<double> Aineq(collision, sparse_size);
		Aineq.setFromTriplets(Aineq_.begin(), Aineq_.end());
		Eigen::SparseMatrix<double> Aeq(aeqsize, sparse_size);
		Aeq.setFromTriplets(Aeq_.begin(), Aeq_.end());
		VectorXd bineq(collision);
		VectorXd beq(aeqsize);
		bineq.setZero();
		beq.setZero();

		for (int i = 0; i < ineq_move_index.size(); i++) {
			bineq(ineq_move_index[i]) = ineq_move[i];
		}

		for (int i = 0; i < eq_move_index.size(); i++) {
			beq(eq_move_index[i]) = eq_move[i];
		}

		if (step > pull_step) {
			for (int i = 0; i < beq_.size(); i++) {
				beq(i) = beq_[i];
			}
		}

		if (matlab_debug_physics) {
			ofstream ofs;
			ofs.open("test.m", ofstream::out | ofstream::trunc);
			ofs.close();
			double_to_file(h, "h");
			double_to_file(grav(2), "grav");
			double_to_file(material.density, "rho");
			double_to_file(e, "e");
			double_to_file(nu, "nu");
			mat_to_file(x_X, "x_X");
			vec_to_file(isEOL, "isEol");
			VectorXi vvv(3);
			vvv << 1, 1, 1;
			mat_to_file(faces2.colwise() += vvv, "faces");
			sparse_to_file_as_sparse_m(MDK, "MDK");
			vec_to_file(b, "b");
			//sparse_to_file_as_dense(N, "N");
			//sparse_to_file_as_dense(M, "M");
			//sparse_to_file_as_sparse_m(N, "N");
			sparse_to_file_as_sparse_m(M, "M");
			vec_to_file(f, "f");
			//vec_to_file(lc, "lc");
			//vec_to_file(uc, "uc");
			//vec_to_file(lx, "lx");
			//vec_to_file(ux, "ux");
			vec_to_file(v_old, "v_last_step");
			vec_to_file(v, "v_transfer");
			sparse_to_file_as_sparse_m(Aineq, "Aineq");
			sparse_to_file_as_sparse_m(Aeq, "Aeq");
			vec_to_file(bineq, "bineq");
			vec_to_file(beq, "beq");
		}

		//igl::mosek::mosek_quadprog(MDK, b, 0, N, lc, uc, lx, ux, mosek_data, v);
		QuadProgMosek *program = new QuadProgMosek();
		double inf = std::numeric_limits<double>::infinity();
		program->setNumberOfVariables(sparse_size);
		program->setNumberOfInequalities(collision);
		program->setNumberOfEqualities(aeqsize);
		program->setObjectiveMatrix(MDK);
		program->setObjectiveVector(b);
		program->setInequalityMatrix(Aineq);
		program->setInequalityVector(bineq);
		program->setEqualityMatrix(Aeq);
		program->setEqualityVector(beq);
		//program->setParamInt(MSK_IPAR_LOG, 0);

		bool success = program->solve();

		v = program->getPrimalSolution();

		if (matlab_debug_physics) {
		//	sparse_to_file_as_dense(MDK, "MDK");
		//	sparse_to_file_as_dense(b, "b");
		//	sparse_to_file_as_dense(N, "N");
		//	sparse_to_file_as_dense(lc, "lc");
		//	sparse_to_file_as_dense(uc, "uc");
		//	sparse_to_file_as_dense(lx, "lx");
		//	sparse_to_file_as_dense(ux, "ux");
		//	sparse_to_file_as_dense(v, "v");
			vec_to_file(v, "v_solved");
		}
	}
	ClothTimer[3]->toc();

	for (int i = 0; i < mesh.nodes.size(); i++) {
		mesh.nodes[i]->v[0] = v(i * 3);
		mesh.nodes[i]->v[1] = v(i * 3 + 1);
		mesh.nodes[i]->v[2] = v(i * 3 + 2);
		mesh.nodes[i]->x = mesh.nodes[i]->x + h * mesh.nodes[i]->v;
	}

	if (export_timings) {
		ClothTimer[0]->export_csv();
		ClothTimer[1]->export_csv();
		ClothTimer[2]->export_csv();
		ClothTimer[3]->export_csv();
		EOL_outputter(mesh.nodes.size(), mesh.faces.size(), mesh.nodes.size(), mesh.faces.size());
	}
}

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

//double calc_edge_weight(Node* n)
//{
//	Node* node0;
//	Node* node1;
//	bool first = true;
//	for (int i = 0; i < n->adje.size(); i++) {
//		if (n->adje[i]->preserve && first) {
//			if (n->adje[i]->n[0] != n) {
//				node0 = n->adje[i]->n[0];
//			}
//			else {
//				node0 = n->adje[i]->n[1];
//			}
//			first = false;
//		}
//		else if (n->adje[i]->preserve && !first) {
//			if (n->adje[i]->n[0] != n) {
//				node1 = n->adje[i]->n[0];
//			}
//			else {
//				node1 = n->adje[i]->n[1];
//			}
//			break;
//		}
//	}
//	if(node0->on_corner) n->which_edge = node1->which_edge;
//	else n->which_edge = node0->which_edge; // MOVE?
//	double d1 = sqrt(pow(node1->x[0] - node0->x[0], 2) + pow(node1->x[1] - node0->x[1], 2) + pow(node1->x[2] - node0->x[2], 2));
//	double d2 = sqrt(pow(node1->x[0] - n->x[0], 2) + pow(node1->x[1] - n->x[1], 2) + pow(node1->x[2] - n->x[2], 2));
//	double ratio = d2 / d1;
//	return abs((ratio * (node1->verts[0]->egde_weight - node0->verts[0]->egde_weight)) - node1->verts[0]->egde_weight);
//}

// EoL physics step
bool Cloth::stepEoL(double h, const Vector3d &grav, vector<shared_ptr<Box> > box, double t, bool equality)
{
	//igl::mosek::MosekData mosek_data;

	int collision = 0; // TODO:: Move

	int eol_edges = 0;

	int eols = 0;
	int LAG_constraints = 0;
	int EOL_constraints = 0;
	vector<double> fill_lc;
	vector<double> fill_uc;

	bool fixed_points = false;  // TODO:: Restructure

	vector<T> M_;
	vector<T> MDK_;
	vector<T> Aineq_;
	vector<T> Aeq_;
	int aineqsize = 0;
	int aeqsize = 0;

	//vector<T> N_;
	N_.clear();
	int sparse_size = 0;

	// Box Setup
	//collisions.clear();
	EOLverts.resize(mesh.nodes.size());
	EOLverts.setZero();
	MatrixXd verts2(3, mesh.nodes.size());
	MatrixXi faces2(3, mesh.faces.size());
	vector<int> ineq_move_index;
	vector<double> ineq_move;
	vector<int> eq_move_index;
	vector<double> eq_move;

	for (int i = 0; i < mesh.nodes.size(); i++) {
		//mesh.nodes[i]->eol_index = -1;
		//mesh.nodes[i]->EoL = false;
		//mesh.nodes[i]->EoL_secondary = false;
		mesh.nodes[i]->preserve_once = false;
		verts2.col(i) = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
		if (mesh.nodes[i]->EoL) {
			if (points) {
				// Remove boundary EOLs hack
				if (mesh.verts[i]->u[0] < 0.01) {
					mesh.nodes[i]->EoL = false;
					mesh.nodes[i]->preserve = false;
					mesh.nodes[i]->EoL_secondary = false;
					mesh.nodes[i]->on_preserved_edge = false;
					mesh.nodes[i]->on_corner = false;
					mesh.nodes[i]->eol_index = -1;
					mesh.nodes[i]->verts[0]->v = Vec3(0); // Just in case
					mesh.nodes[i]->on_preserved_edge = false;
					continue;
				}
			}
			EOLverts(i) = 1;
			mesh.nodes[i]->eol_index = eols;
			eols++;
			if (wire) {
				fill_lc.push_back(0.0);
				fill_uc.push_back(+MSK_INFINITY);
				LAG_constraints++;
				fill_lc.push_back(0.0);
				fill_uc.push_back(0.0);
				LAG_constraints++;
				Aineq_.push_back(T(collision + aineqsize, i * 3, 0.0)); // Inequality
				Aineq_.push_back(T(collision + aineqsize, i * 3 + 1, 0.0));
				Aineq_.push_back(T(collision + aineqsize, i * 3 + 2, -1.0));
				aineqsize++;
				// Edge movement hack
				if (move_wire) {
					ineq_move_index.push_back(collision + aineqsize - 1);
					ineq_move.push_back(-1.0);
				}
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3, 0.0));
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 1, -1.0));
				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 2, 0.0));
				aeqsize++;
				mesh.nodes[i]->n = Vec3(0.0, 0.0, 1.0);
			}
		}
		else {
			//if (mesh.verts[i]->u[0] > 0.05 && mesh.verts[i]->u[0] < 1-0.05 && mesh.verts[i]->u[1] > 0.05 && mesh.verts[i]->u[1] < 1-0.05) {
				int preservenum = 0;
				Node* node0;
				Node* node1;
				for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
					if (mesh.nodes[i]->adje[j]->preserve) preservenum++;
				}
				if (preservenum == 2) {
					
					//mesh.nodes[i]->verts[0]->egde_weight[0] = calc_edge_weight(mesh.nodes[i]);
					cout << "ArcSim point added. Index: " << i << "Edge: " << mesh.nodes[i]->which_edge << " Weight: " << mesh.nodes[i]->verts[0]->egde_weight << endl;
					mesh.nodes[i]->EoL = true;
					mesh.nodes[i]->EoL_secondary = true;
					mesh.nodes[i]->on_preserved_edge = true;
					EOLverts(i) = 1;
					mesh.nodes[i]->eol_index = eols;
					eols++;
					if (wire) {
						fill_lc.push_back(0.0);
						fill_uc.push_back(+MSK_INFINITY);
						LAG_constraints++;
						fill_lc.push_back(0.0);
						fill_uc.push_back(0.0);
						LAG_constraints++;
						Aineq_.push_back(T(collision + aineqsize, i * 3, 0.0)); // Inequality
						Aineq_.push_back(T(collision + aineqsize, i * 3 + 1, 0.0));
						Aineq_.push_back(T(collision + aineqsize, i * 3 + 2, -1.0));
						aineqsize++;
						// Edge movement hack
						if (move_wire) {
							ineq_move_index.push_back(collision + aineqsize - 1);
							ineq_move.push_back(-1.0);
						}
						Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3, 0.0));
						Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 1, -1.0));
						Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 2, 0.0));
						aeqsize++;
						mesh.nodes[i]->n = Vec3(0.0, 0.0, 1.0);
					}
				}
				else {
					mesh.nodes[i]->eol_index = -1;
					mesh.nodes[i]->verts[0]->v = Vec3(0); // Just in case
					mesh.nodes[i]->on_preserved_edge = false;
				}
				//if (wire) {
				//	if (get_angle(Vec3(0.0, 0.0, 1.0), mesh.nodes[i]->n) > 0.5) {
				//		if (mesh.nodes[i]->x[1] > 0.515 && dot(mesh.nodes[i]->n, Vec3(0.0, 1.0, 0.0)) < 0.0) {
				//			mesh.nodes[i]->x[1] = 0.505;
				//			Aineq_.push_back(T(collision + aineqsize, i * 3, -mesh.nodes[i]->n[0])); // Inequality
				//			Aineq_.push_back(T(collision + aineqsize, i * 3 + 1, -mesh.nodes[i]->n[1]));
				//			Aineq_.push_back(T(collision + aineqsize, i * 3 + 2, -mesh.nodes[i]->n[2]));
				//			aineqsize++;
				//		}
				//		else if (mesh.nodes[i]->x[1] < 0.505 && dot(mesh.nodes[i]->n, Vec3(0.0, -1.0, 0.0)) < 0.0) {
				//			mesh.nodes[i]->x[1] = 0.505;
				//			Aineq_.push_back(T(collision + aineqsize, i * 3, -mesh.nodes[i]->n[0])); // Inequality
				//			Aineq_.push_back(T(collision + aineqsize, i * 3 + 1, -mesh.nodes[i]->n[1]));
				//			Aineq_.push_back(T(collision + aineqsize, i * 3 + 2, -mesh.nodes[i]->n[2]));
				//			aineqsize++;
				//		}
				//	}
				//}
			//}
			//else {
			//	mesh.nodes[i]->eol_index = -1;
			//}
		}
	}
	for (int i = 0; i < mesh.faces.size(); i++) {
		faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}

	//MatrixXd verts1(3, 5);
	//verts1.col(0) = Vector3d(0.51, 0.51, -0.1);
	//verts1.col(1) = Vector3d(0.21, 0.21, -0.08);
	//verts1.col(2) = Vector3d(0.3, 0.51, -0.1);
	//verts1.col(3) = Vector3d(0.71, 0.21, -0.1);
	//verts1.col(4) = Vector3d(0.51, 0.3, -0.12);
	//MatrixXd norms1(3, 5);
	//norms1.col(0) = Vector3d(0.0, 0.0, 1.0);
	//norms1.col(1) = Vector3d(0.0, 0.0, 1.0);
	//norms1.col(2) = Vector3d(0.0, 0.0, 1.0);
	//norms1.col(3) = Vector3d(0.0, 0.0, 1.0);
	//norms1.col(4) = Vector3d(0.0, 0.0, 1.0);
	/*vector<int> ineq_move_index;
	vector<double> ineq_move;*/

	ClothTimer[0]->tic();
	for (int b = 0; b < box.size(); b++) {
		collisions.clear();
		if(!wire && !points) boxTriCollision(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, EOLverts, true);
		else if (points) {
			pointTriCollision(collisions, box[b]->thresh, verts1, norms1, verts2, faces2, true);
			//pointTriCollisionBoxWalls(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, EOLverts, true);
			boxTriCollision(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, EOLverts, false);
			if(move_wire) verts1 += norms1 * pmove1 * h;
		}
		if (wire && move_wire) {
			box[b]->E1(2, 3) += 1.0 * 2.5 * h;
			box[b]->x(2) += 1.0 * 2.5 * h;
		}
		if (collisions.size() > 0) {
			//if (matlab_debug_physics) {
			//	double_to_file(box->thresh, "threshold");
			//	vec_to_file(box->dim, "whd1");
			//	mat_to_file(box->E1, "E1");
			//	mat_to_file(verts2, "verts2");
			//	VectorXi vvv(3);
			//	vvv << 1, 1, 1;
			//	mat_to_file(faces2.colwise() += vvv, "faces2");
			//}
			shared_ptr<Rigid> rigid;
			for (int i = 0; i < collisions.size(); i++) {
				if (collisions[i]->count1 == 3 && collisions[i]->count2 == 1) {
					if (mesh.nodes[collisions[i]->verts2(0)]->EoL) {
						if (mesh.nodes[collisions[i]->verts2(0)]->on_corner) {
							Vec3 nor_ave = Vec3(0.0, 0.0, 0.0);
							for (int j = 0; j < mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size(); j++) {
								nor_ave += (mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf[j]->n * incedent_angle(mesh.nodes[collisions[i]->verts2(0)]->verts[0], mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf[j]));
							}
							nor_ave[0] = nor_ave[0] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
							nor_ave[1] = nor_ave[1] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
							nor_ave[2] = nor_ave[2] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
							//normalize(nor_ave);
							//collisions[i]->nor1 = Vector3d(nor_ave[0], nor_ave[1], nor_ave[2]).normalized();
							//collisions[i]->nor1 = Vector3d(0.5774, -0.5774, 0.5774);
							//Vector3d fnxt = collisions[i]->nor1.cross(Vector3d(1.0, 0.0, 0.0));
							//Vector3d fnxb = collisions[i]->nor1.cross(fnxt);


							if (!points) collisions[i]->nor1 = Vector3d(mesh.nodes[collisions[i]->verts2(0)]->nor_ave[0], mesh.nodes[collisions[i]->verts2(0)]->nor_ave[1], mesh.nodes[collisions[i]->verts2(0)]->nor_ave[2]);
							if (points) collisions[i]->nor1 = Vector3d(0.0, 0.0, 1.0);
							//if(!points) collisions[i]->nor1 = Vector3d(0.0, 0.0, 1.0); // Up ineq hack
							Vector3d fnxt;
							if (fabs(collisions[i]->nor1(0)) - 1.0 > 1e-3) {
								fnxt << 1.0, 0.0, 0.0;
							}
							else {
								fnxt << 0.0, 1.0, 0.0;
							}
							Vector3d fnxb = (fnxt.cross(collisions[i]->nor1)).normalized();
							fnxt = (collisions[i]->nor1.cross(fnxb)).normalized();
							collisions[i]->tan1 = fnxt;
							collisions[i]->bin1 = fnxb;

							if (equality) {
								Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0))); // Equality hack
								Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
								Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
								aeqsize++;
							}
							else {
								Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3, -collisions[i]->nor1(0))); // Inequality
								Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3 + 1, -collisions[i]->nor1(1)));
								Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3 + 2, -collisions[i]->nor1(2)));
								aineqsize++;
							}
							// Point movement hack
							if (points && move_wire) {
								ineq_move_index.push_back(collision + aineqsize - 1);
								ineq_move.push_back(-pmove1);
								//ineq_move.push_back(-1.0);
							}
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3, fnxt(0)));
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, fnxt(1)));
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, fnxt(2)));
							aeqsize++;
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3, fnxb(0)));
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, fnxb(1)));
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, fnxb(2)));
							aeqsize++;
							if(friction) mesh.nodes[collisions[i]->verts2(0)]->coll_case = 3;
							mesh.nodes[collisions[i]->verts2(0)]->lookout_for = true;

							if (move_wire && !points) {
								// Moving rigid body
								Vector4d xl = box[b]->E1inv*Vector4d(collisions[i]->pos1(0), collisions[i]->pos1(1), collisions[i]->pos1(2), 1.0);
								//MatrixXd gamma = rigid->gamma(xl.segment<3>(0));
								//Matrix3d theta = box[b]->E1.block<3, 3>(0, 0);
								//Vector3d xdot = theta * gamma * box[b]->v;
								//cout << xdot << endl;
								Matrix4d et = rigid->integrate(box[b]->E1, box[b]->v, h);
								Vector4d xdot = ((et*xl) - (box[b]->E1*xl)) / h;
								double temp = collisions[i]->nor1.transpose()*xdot.segment<3>(0);
								if (equality) {
									eq_move_index.push_back(aeqsize + EOL_constraints - 3);
									eq_move.push_back(temp);
								}
								else {
									ineq_move_index.push_back(collision + aineqsize - 1);
									ineq_move.push_back(-temp);
								}
								eq_move_index.push_back(aeqsize + EOL_constraints - 2);
								eq_move_index.push_back(aeqsize + EOL_constraints - 1);
								temp = fnxt.transpose()*xdot.segment<3>(0);
								eq_move.push_back(temp);
								temp = fnxb.transpose()*xdot.segment<3>(0);
								eq_move.push_back(temp);
							}

							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
							//fill_lc.push_back(0.0);
							//fill_uc.push_back(0.0);
							//LAG_constraints++;
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, collisions[i]->bin1(0)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, collisions[i]->bin1(1)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, collisions[i]->bin1(2)));
							//fill_lc.push_back(0.0);
							//fill_uc.push_back(0.0);
							//LAG_constraints++;
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, collisions[i]->tan1(0)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, collisions[i]->tan1(1)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, collisions[i]->tan1(2)));
							//fill_lc.push_back(0.0);
							//fill_uc.push_back(0.0);
							//LAG_constraints++;
						}
						else {
							Vec3 nor_ave = Vec3(0.0, 0.0, 0.0);
							for (int j = 0; j < mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size(); j++) {
								nor_ave += (mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf[j]->n * incedent_angle(mesh.nodes[collisions[i]->verts2(0)]->verts[0], mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf[j]));
							}
							nor_ave[0] = nor_ave[0] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
							nor_ave[1] = nor_ave[1] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
							nor_ave[2] = nor_ave[2] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
							//normalize(nor_ave);
							//collisions[i]->nor1 = Vector3d(nor_ave[0], nor_ave[1], nor_ave[2]).normalized();
							// Average Face normals	
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
							fill_lc.push_back(0.0);
							fill_uc.push_back(+MSK_INFINITY);
							LAG_constraints++;
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
							//fill_lc.push_back(0.0);
							//fill_uc.push_back(0.0);
							//LAG_constraints++;
							//collisions[i]->nor1 = Vector3d(0.0, 0.0, 1.0); // Up ineq hack
							Vector3d fnxt = (collisions[i]->nor1.cross(collisions[i]->tan1)).normalized();
							//fnxt = Vector3d(0.0, -1.0, 0.0); // Up ineq hack
							collisions[i]->bin1 = fnxt;
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, fnxt(0)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, fnxt(1)));
							//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, fnxt(2)));
							fill_lc.push_back(0.0);
							fill_uc.push_back(0.0);
							LAG_constraints++;
							
							if (equality) {
								Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0))); // Equality hack
								Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
								Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
								aeqsize++;
							}
							else {
								Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3, -collisions[i]->nor1(0))); // Inequality
								Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3 + 1, -collisions[i]->nor1(1)));
								Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3 + 2, -collisions[i]->nor1(2)));
								aineqsize++;
							}
							// Edge movement hack
							//if (wire) {
							//	ineq_move_index.push_back(collision + aineqsize - 1);
							//	ineq_move.push_back(-1.0);
							//}
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3, fnxt(0)));
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, fnxt(1)));
							Aeq_.push_back(T(aeqsize + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, fnxt(2)));
							aeqsize++;
							//if (mesh.nodes[collisions[i]->verts2(0)]->preserve == true) mesh.nodes[collisions[i]->verts2(0)]->preserve = false; // Was this needed?
							if (friction) {
								mesh.nodes[collisions[i]->verts2(0)]->coll_case = 2;
								mesh.nodes[collisions[i]->verts2(0)]->nor_ave = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(1));
							}
							mesh.nodes[collisions[i]->verts2(0)]->lookout_for = true; // For lifting off EOLS

							if (move_wire) {
								// Moving rigid body
								Vector4d xl = box[b]->E1inv*Vector4d(collisions[i]->pos1(0), collisions[i]->pos1(1), collisions[i]->pos1(2), 1.0);
								//MatrixXd gamma = rigid->gamma(xl.segment<3>(0));
								//Matrix3d theta = box[b]->E1.block<3, 3>(0, 0);
								//Vector3d xdot = theta * gamma * box[b]->v;
								//cout << xdot << endl;
								Matrix4d et = rigid->integrate(box[b]->E1, box[b]->v, h);
								Vector4d xdot = ((et*xl) - (box[b]->E1*xl)) / h;
								double temp = collisions[i]->nor1.transpose()*xdot.segment<3>(0);
								if (equality) {
									eq_move_index.push_back(aeqsize + EOL_constraints - 2);
									eq_move.push_back(temp);
								}
								else {
									ineq_move_index.push_back(collision + aineqsize - 1);
									ineq_move.push_back(-temp);
								}
								eq_move_index.push_back(aeqsize + EOL_constraints - 1);
								temp = fnxt.transpose()*xdot.segment<3>(0);
								eq_move.push_back(temp);
							}
						}

						//Vec3 nor_ave = Vec3(0.0, 0.0, 0.0);
						//for (int j = 0; j < mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size(); j++) {
						//	nor_ave += mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf[j]->n;
						//}
						//nor_ave[0] = nor_ave[0] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
						//nor_ave[1] = nor_ave[1] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
						//nor_ave[2] = nor_ave[2] * (1.0 / mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size());
						//normalize(nor_ave);
						////nor_ave = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(2));
						//nor_ave = mesh.nodes[collisions[i]->verts2(0)]->coll_norm;
						//mesh.nodes[collisions[i]->verts2(0)]->nor_ave = nor_ave;
						//// Average Face normals	
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, nor_ave[0]));
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, nor_ave[1]));
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, nor_ave[2]));
						//fill_lc.push_back(0.0);
						//fill_uc.push_back(0.0);
						//LAG_constraints++;
						//// Tangent
						////Vector3d ave = Vector3d(nor_ave[0], nor_ave[1], nor_ave[2]);
						////Vector3d b1 = collisions[i]->tan1.cross(ave);
						//Vec3 b1 = cross(mesh.nodes[collisions[i]->verts2(0)]->tan, nor_ave);
						//mesh.nodes[collisions[i]->verts2(0)]->b1 = Vec3(b1[0], b1[1], b1[2]);
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, b1[0]));
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, b1[1]));
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, b1[2]));
						//fill_lc.push_back(0.0);
						//fill_uc.push_back(0.0);
						//LAG_constraints++;
						////
						////N_.push_back(T(collision, collisions[i]->verts2(0) * 3, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[0]));
						////N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[1]));
						////N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, mesh.nodes[collisions[i]->verts2(0)]->coll_norm[2]));
						////
						//mesh.nodes[collisions[i]->verts2(0)]->EoL = true;
						//mesh.nodes[collisions[i]->verts2(0)]->eol_index = eols;
						//eols++;
						////for (int j = 0; j < mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf.size(); j++) {
						////	for (int k = 0; k < 3; k++) {
						////		if (!mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf[j]->v[k]->node->EoL_secondary) {
						////			mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf[j]->v[k]->node->EoL_secondary = true;
						////			mesh.nodes[collisions[i]->verts2(0)]->verts[0]->adjf[j]->v[k]->node->eol_index = eols;
						////			eols++;
						////		}
						////	}
						////}
						//if (mesh.nodes[collisions[i]->verts2(0)]->preserve == true) mesh.nodes[collisions[i]->verts2(0)]->preserve = false;
						//mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = false;
					}
					else {
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
						//N_.push_back(T(collision + LAG_constraints + EOL_constraints, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
						if (!wire) {
							fill_lc.push_back(0.0);
							fill_uc.push_back(+MSK_INFINITY);
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3, -collisions[i]->nor1(0)));
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3 + 1, -collisions[i]->nor1(1)));
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(0) * 3 + 2, -collisions[i]->nor1(2)));
							if (friction) {
								mesh.nodes[collisions[i]->verts2(0)]->coll_case = 1;
								mesh.nodes[collisions[i]->verts2(0)]->nor_ave = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(1));
							}
							collision++;
							if (move_wire) {
								// Moving rigid body
								Vector4d xl = box[b]->E1inv*Vector4d(collisions[i]->pos1(0), collisions[i]->pos1(1), collisions[i]->pos1(2), 1.0);
								//MatrixXd gamma = rigid->gamma(xl.segment<3>(0));
								//Matrix3d theta = box[b]->E1.block<3, 3>(0, 0);
								//Vector3d xdot = theta * gamma * box[b]->v;
								Matrix4d et = rigid->integrate(box[b]->E1, box[b]->v, h);
								Vector4d xdot = ((et*xl) - (box[b]->E1*xl)) / h;
								ineq_move_index.push_back(collision + aineqsize - 1);
								double temp = collisions[i]->nor1.transpose()*xdot.segment<3>(0);
								ineq_move.push_back(-temp);
							}
						}
						// Point movement hack
						//ineq_move_index.push_back(collision + aineqsize - 1);
						//ineq_move.push_back(-1.0);
					}
				}
				if (points) {
					if (collisions[i]->count1 == 2 && collisions[i]->count2 == 2) {
						for (int j = 0; j < 2; j++) {
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(j) * 3, -collisions[i]->nor2(0) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(j) * 3 + 1, -collisions[i]->nor2(1) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(j) * 3 + 2, -collisions[i]->nor2(2) * collisions[i]->weights2(j)));
						}
						collision++;
					}
					else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 3) {
						for (int j = 0; j < 3; j++) {
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(j) * 3, -collisions[i]->nor2(0) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(j) * 3 + 1, -collisions[i]->nor2(1) * collisions[i]->weights2(j)));
							Aineq_.push_back(T(collision + aineqsize, collisions[i]->verts2(j) * 3 + 2, -collisions[i]->nor2(2) * collisions[i]->weights2(j)));
						}
						collision++;
					}
				}
			}
		}
	}
	ClothTimer[0]->toc();

	// For margin collisions
	//for (int i = 0; i < collisions_passed.size(); i++) {
	//	N_.push_back(T(collision, collisions[i]->verts2(0) * 3, collisions[i]->nor1(0)));
	//	N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 1, collisions[i]->nor1(1)));
	//	N_.push_back(T(collision, collisions[i]->verts2(0) * 3 + 2, collisions[i]->nor1(2)));
	//	collision++;
	//}

	v_old.resize((mesh.nodes.size() * 3) + (eols * 2));
	v.resize((mesh.nodes.size() * 3) + (eols * 2));
	f.resize((mesh.nodes.size() * 3) + (eols * 2));

	v_old.setZero();
	v.setZero();
	f.setZero();

	MatrixXd x_X(mesh.nodes.size(), 5);
	VectorXi isEoL(mesh.nodes.size());
	isEoL.setZero();

	MatrixXd regulizeV = MatrixXd::Zero(mesh.nodes.size()*3 + eols * 2, mesh.nodes.size()*3 + eols * 2);

	if (step == 412) {
		cout << endl;
	}

	ClothTimer[1]->tic();
	for (int i = 0; i < mesh.nodes.size(); i++) {

		x_X(i, 0) = mesh.nodes[i]->x[0];
		x_X(i, 1) = mesh.nodes[i]->x[1];
		x_X(i, 2) = mesh.nodes[i]->x[2];
		x_X(i, 3) = mesh.nodes[i]->verts[0]->u[0];
		x_X(i, 4) = mesh.nodes[i]->verts[0]->u[1];

		// If EOL but wasn't detected by the CD we need to turn it off
		if (!wire) {
			if (mesh.nodes[i]->EoL && !mesh.nodes[i]->lookout_for) {
				mesh.nodes[i]->EoL = false;
				mesh.nodes[i]->EoL_secondary = false;
				mesh.nodes[i]->on_preserved_edge = false;
				mesh.nodes[i]->on_corner = false;
				mesh.nodes[i]->coll_case = -1;
				for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
					mesh.nodes[i]->adje[j]->preserve = false;
				}
				mesh.nodes[i]->preserve_once = true;
				mesh.nodes[i]->preserve = true;
				Aeq_.push_back(T(aeqsize + EOL_constraints, mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2, 1));
				aeqsize++;
				Aeq_.push_back(T(aeqsize + EOL_constraints, mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1, 1));
				aeqsize++;
				cout << "EOL: " << i << " removed" << endl;
			}
			mesh.nodes[i]->lookout_for = false;
		}	

		bool found = false;
		double how_close = 1e-6;
		int closest = -1;
		for (int j = 0; j < last_mesh.nodes.size(); j++) {
			if (unsigned_vv_distance(mesh.nodes[i]->verts[0]->u, last_mesh.nodes[j]->verts[0]->u) < how_close) {
				how_close = unsigned_vv_distance(mesh.nodes[i]->verts[0]->u, last_mesh.nodes[j]->verts[0]->u);
				closest = j;

				//v_old(3 * i) = last_mesh.nodes[j]->v[0];
				//v_old(3 * i + 1) = last_mesh.nodes[j]->v[1];
				//v_old(3 * i + 2) = last_mesh.nodes[j]->v[2];
				//if (last_mesh.nodes[j]->EoL) {
				//	v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) = last_mesh.nodes[j]->verts[0]->v[0];
				//	v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) = last_mesh.nodes[j]->verts[0]->v[1];
				//}
				found = true;
				break;
			}
		}
		if (found) {
			if (mesh.nodes[i]->preserve_once && !mesh.nodes[i]->EoL) {
				MatrixXd F = MatrixXd::Zero(3, 2);
				Vector3d v = Vector3d(mesh.nodes[i]->v[0], mesh.nodes[i]->v[1], mesh.nodes[i]->v[2]);;
				Vector2d V = Vector2d(mesh.nodes[i]->verts[0]->v[0], mesh.nodes[i]->verts[0]->v[1]);
				for (int j = 0; j < mesh.nodes[i]->verts[0]->adjf.size(); j++) {
					F += incedent_angle(mesh.nodes[i]->verts[0], mesh.nodes[i]->verts[0]->adjf[j]) * deform_grad(mesh.nodes[i]->verts[0]->adjf[j]);
				}
				is_seam_or_boundary(mesh.nodes[i]) ? F *= (1 / M_PI) : F *= (1 / (2 * M_PI));
				//F /= mesh.nodes[i]->verts[0]->adjf.size();
				Vector3d newV = v - F*V;
				v_old(3 * i) = newV(0);
				v_old(3 * i + 1) = newV(1);
				v_old(3 * i + 2) = newV(2);
				v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) = 0.0;
				v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) = 0.0;
			}
			else {
				v_old(3 * i) = last_mesh.nodes[closest]->v[0];
				v_old(3 * i + 1) = last_mesh.nodes[closest]->v[1];
				v_old(3 * i + 2) = last_mesh.nodes[closest]->v[2];
				if (last_mesh.nodes[closest]->EoL) {
					v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) = last_mesh.nodes[closest]->verts[0]->v[0];
					v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) = last_mesh.nodes[closest]->verts[0]->v[1];
				}
			}
		}

		if (!found && !mesh.nodes[i]->EoL) {
			// Get new LAG world velocities
			Face* old_face = get_enclosing_face(last_mesh, Vec2(mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->verts[0]->u[1]));
			Vec3 bary = get_barycentric_coords(Vec2(mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->verts[0]->u[1]), old_face);
			//Vec3 baryv = bary[0] * old_face->v[0]->node->v + bary[1] * old_face->v[1]->node->v + bary[2] * old_face->v[2]->node->v;
			//Vec3 baryV = bary[0] * old_face->v[0]->v + bary[1] * old_face->v[1]->v + bary[2] * old_face->v[2]->v;
			Vector3d vwA;
			Vector3d vwB;
			Vector3d vwC;
			if (old_face->v[0]->node->EoL) {
				vwA = Vector3d(old_face->v[0]->node->v[0], old_face->v[0]->node->v[1], old_face->v[0]->node->v[2]) + -deform_grad(old_face) * Vector2d(old_face->v[0]->v[0], old_face->v[0]->v[1]);
			}
			else {
				vwA = Vector3d(old_face->v[0]->node->v[0], old_face->v[0]->node->v[1], old_face->v[0]->node->v[2]);
			}
			if (old_face->v[1]->node->EoL) {
				vwB = Vector3d(old_face->v[1]->node->v[0], old_face->v[1]->node->v[1], old_face->v[1]->node->v[2]) + -deform_grad(old_face) * Vector2d(old_face->v[1]->v[0], old_face->v[1]->v[1]);
			}
			else {
				vwB = Vector3d(old_face->v[1]->node->v[0], old_face->v[1]->node->v[1], old_face->v[1]->node->v[2]);
			}
			if (old_face->v[2]->node->EoL) {
				vwC = Vector3d(old_face->v[2]->node->v[0], old_face->v[2]->node->v[1], old_face->v[2]->node->v[2]) + -deform_grad(old_face) * Vector2d(old_face->v[2]->v[0], old_face->v[2]->v[1]);
			}
			else {
				vwC = Vector3d(old_face->v[2]->node->v[0], old_face->v[2]->node->v[1], old_face->v[2]->node->v[2]);
			}
			//Vector3d baryvv = Vector3d(baryv[0], baryv[1], baryv[2]);
			//Vector2d baryVV = Vector2d(baryV[0], baryV[1]);
			//Vector3d v_new_world = baryvv + deform_grad(old_face) * baryVV;
			Vector3d v_new_world = bary[0] * vwA + bary[1] * vwB + bary[2] * vwC;
			mesh.nodes[i]->v = Vec3(v_new_world(0), v_new_world(1), v_new_world(2));
			v_old(3 * i) = v_new_world(0);
			v_old(3 * i + 1) = v_new_world(1);
			v_old(3 * i + 2) = v_new_world(2);
		}

		// Fill velocity vector for x
		//v(3 * i) = mesh.nodes[i]->v[0];
		//v(3 * i + 1) = mesh.nodes[i]->v[1];
		//v(3 * i + 2) = mesh.nodes[i]->v[2];
		sparse_size += 3;

		if (mesh.nodes[i]->EoL) {
			regulizeV.block(i * 3, i * 3, 3, 3) = Matrix3d::Identity();
			//MDK_.push_back(T(i * 3, i * 3, 1e-6));
			//MDK_.push_back(T(i * 3 + 1, i * 3 + 1, 1e-6));
			//MDK_.push_back(T(i * 3 + 2, i * 3 + 2, 1e-6));
			//MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2, mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2, 1e-6));
			//MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + 1, mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + 1, 1e-6));
			isEoL(i) = 1;
			//v(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) = mesh.nodes[i]->verts[0]->v[0];
			//v(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) = mesh.nodes[i]->verts[0]->v[1];
			// EoL constraint
			if (is_seam_or_boundary(mesh.nodes[i])) continue;
			Vector2d tan1;
			Vector2d tan2;
			Node* n1;
			Node* n2;
			bool first = true;
			bool second = true;
			for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
				if (mesh.nodes[i]->adje[j]->preserve) {
					if (first) {
						if (mesh.nodes[i] == mesh.nodes[i]->adje[j]->n[0]) {
							tan1 = Vector2d(mesh.nodes[i]->adje[j]->n[1]->verts[0]->u[0] - mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->adje[j]->n[1]->verts[0]->u[1] - mesh.nodes[i]->verts[0]->u[1]);
							n1 = mesh.nodes[i]->adje[j]->n[1];
						}
						else {
							tan1 = Vector2d(mesh.nodes[i]->adje[j]->n[0]->verts[0]->u[0] - mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->adje[j]->n[0]->verts[0]->u[1] - mesh.nodes[i]->verts[0]->u[1]);
							n1 = mesh.nodes[i]->adje[j]->n[0];
						}
						first = false;
					}
					else {
						if (mesh.nodes[i] == mesh.nodes[i]->adje[j]->n[0]) {
							tan2 = -Vector2d(mesh.nodes[i]->adje[j]->n[1]->verts[0]->u[0] - mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->adje[j]->n[1]->verts[0]->u[1] - mesh.nodes[i]->verts[0]->u[1]);
							n2 = mesh.nodes[i]->adje[j]->n[1];
						}
						else {
							tan2 = -Vector2d(mesh.nodes[i]->adje[j]->n[0]->verts[0]->u[0] - mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->adje[j]->n[0]->verts[0]->u[1] - mesh.nodes[i]->verts[0]->u[1]);
							n2 = mesh.nodes[i]->adje[j]->n[0];
						}
						second = false; // TODO: More edges than just two
						break;
					}
				}
			}
			// EOL constraints
			if (!mesh.nodes[i]->on_corner) {
				Vector2d ave;
				if(!first && !second) ave = (tan1.normalized() + tan2.normalized()) / 2.0;
				if (!first && second) ave = tan1;
				ave.normalize();
				//N_.push_back(T(collision + LAG_constraints + EOL_constraints, mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2, ave(0)));
				//N_.push_back(T(collision + LAG_constraints + EOL_constraints, mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1, ave(1)));
				if (!first) {
					fill_lc.push_back(0.0);
					fill_uc.push_back(0.0);
					Aeq_.push_back(T(aeqsize + EOL_constraints, mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2, ave(0)));
					Aeq_.push_back(T(aeqsize + EOL_constraints, mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1, ave(1)));
					EOL_constraints++;
				}
			}

			// New velocities from averages of two on each side
			if (!found) {
				if(mesh.nodes[i]->EoL_secondary) {
					//mesh.nodes[i]->v = (n1->v + n2->v) / 2.0;
					//mesh.nodes[i]->verts[0]->v = (n1->verts[0]->v + n2->verts[0]->v) / 2.0;
					mesh.nodes[i]->v = (n1->v + n2->v) / 2.0;
					v_old(3 * i) = (n1->v[0] + n2->v[0]) / 2;
					v_old(3 * i + 1) = (n1->v[1] + n2->v[1]) / 2;
					v_old(3 * i + 2) = (n1->v[2] + n2->v[2]) / 2;
					mesh.nodes[i]->verts[0]->v = (n1->verts[0]->v + n2->verts[0]->v) / 2.0;
					v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) = (n1->verts[0]->v[0] + n2->verts[0]->v[0]) / 2;
					v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) = (n1->verts[0]->v[1] + n2->verts[0]->v[1]) / 2;
					mesh.nodes[i]->EoL_secondary = false;
				}
				else {
					Face* old_face = get_enclosing_face(last_mesh, Vec2(mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->verts[0]->u[1]));
					Matrix2d ftf = deform_grad(old_face).transpose() * deform_grad(old_face);
					Vector2d dtv = -deform_grad(old_face).transpose() * Vector3d(mesh.nodes[i]->v[0], mesh.nodes[i]->v[1], mesh.nodes[i]->v[2]);
					MatrixXd KKTl = MatrixXd::Zero(4, 4);
					VectorXd KKTr(4);
					KKTl.block(0, 0, 2, 2) = ftf;
					KKTl(0, 2) = 1;
					KKTl(1, 3) = 1;
					KKTl(2, 0) = 1;
					KKTl(3, 1) = 1;
					KKTr << dtv, VectorXd::Zero(2);
					ConjugateGradient<MatrixXd, Lower | Upper> cg;
					cg.compute(KKTl);
					VectorXd newv = cg.solve(KKTr);
					mesh.nodes[i]->verts[0]->v = Vec3(newv(0), newv(1), 0.0);
					v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) = newv(0);
					v_old(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) = newv(1);
					Vector3d newvL = Vector3d(mesh.nodes[i]->v[0], mesh.nodes[i]->v[1], mesh.nodes[i]->v[2]) + -deform_grad(old_face) * newv.segment<2>(0);
					mesh.nodes[i]->v = Vec3(newvL(0), newvL(1), newvL(2));
					v_old(3 * i) = newvL(0);
					v_old(3 * i + 1) = newvL(1);
					v_old(3 * i + 2) = newvL(2);
					cout << "newly EOL: " << i << endl;
				}
			}
		}
	}
	ClothTimer[1]->toc();

	ClothTimer[2]->tic();
	for (int i = 0; i < mesh.faces.size(); i++) {
		// Set up corrdaintes
		double xa[3];
		double xb[3];
		double xc[3];
		double Xa[2];
		double Xb[2];
		double Xc[2];

		// Fill temporary vectors 
		// TODO:: Convrsion functions
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
		Map<Vector3d>(xa, 3) = txa;
		Map<Vector3d>(xb, 3) = txb;
		Map<Vector3d>(xc, 3) = txc;
		Map<Vector2d>(Xa, 2) = tXa;
		Map<Vector2d>(Xb, 2) = tXb;
		Map<Vector2d>(Xc, 2) = tXc;
		Map<Vector3d>(g, grav.rows(), grav.cols()) = grav;

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

		int aindex = mesh.faces[i]->v[0]->node->index * 3;
		int bindex = mesh.faces[i]->v[1]->node->index * 3;
		int cindex = mesh.faces[i]->v[2]->node->index * 3;
		int aindexX = mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2;
		int bindexX = mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2;
		int cindexX = mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2;

		// EoL faces
		if ((mesh.faces[i]->v[0]->node->EoL || mesh.faces[i]->v[1]->node->EoL || mesh.faces[i]->v[2]->node->EoL)) {
			MatrixXd F = deform_grad(mesh.faces[i]);
			double fi[15]; // Gravity force vector
			double Mi[225]; //Inertial matrix
			double fm[15];
			double Km[225];
			ComputeMembraneEOL(xa, xb, xc, Xa, Xb, Xc, e, nu, PP, QQ, Wm, fm, Km);
			VectorXd fme = Map<VectorXd>(fm, 15);
			MatrixXd Kme = Map<MatrixXd>(Km, 15, 15);
			ComputeInertiaEOL(xa, xb, xc, Xa, Xb, Xc, g, material.density, Wi, fi, Mi);
			VectorXd fie = Map<VectorXd>(fi, 15);
			MatrixXd Mie = Map<MatrixXd>(Mi, 15, 15);

			// Recompute EOL portions of Kme
			//MatrixXd Kmaa = Kme.block(0, 0, 3, 3); MatrixXd Kmab = Kme.block(0, 5, 3, 3); MatrixXd Kmac = Kme.block(0, 10, 3, 3);
			//Kme.block(5, 0, 3, 3); Kme.block(5, 5, 3, 3); Kme.block(5, 10, 3, 3);
			//Kme.block(10, 0, 3, 3); Kme.block(10, 5, 3, 3); Kme.block(10, 10, 3, 3);

			//Kme.block(0, 3, 3, 2) = -Kme.block(0, 0, 3, 3) * F; Kme.block(0, 8, 3, 2) = -Kme.block(0, 5, 3, 3) * F; Kme.block(0, 13, 3, 2) = -Kme.block(0, 10, 3, 3) * F;
			//Kme.block(5, 3, 3, 2) = -Kme.block(5, 0, 3, 3) * F; Kme.block(5, 8, 3, 2) = -Kme.block(5, 5, 3, 3) * F; Kme.block(5, 13, 3, 2) = -Kme.block(5, 10, 3, 3) * F;
			//Kme.block(10, 3, 3, 2) = -Kme.block(10, 0, 3, 3) * F; Kme.block(10, 8, 3, 2) = -Kme.block(10, 5, 3, 3) * F; Kme.block(10, 13, 3, 2) = -Kme.block(10, 10, 3, 3) * F;
			//Kme.block(3, 0, 2, 3) = -F.transpose() * Kme.block(0, 0, 3, 3); Kme.block(3, 0, 2, 3) = -F.transpose() * Kme.block(0, 5, 3, 3); Kme.block(3, 0, 2, 3) = -F.transpose() * Kme.block(0, 10, 3, 3);
			//Kme.block(8, 0, 2, 3) = -F.transpose() * Kme.block(5, 0, 3, 3); Kme.block(8, 0, 2, 3) = -F.transpose() * Kme.block(5, 5, 3, 3); Kme.block(8, 0, 2, 3) = -F.transpose() * Kme.block(5, 10, 3, 3);
			//Kme.block(13, 0, 2, 3) = -F.transpose() * Kme.block(10, 0, 3, 3); Kme.block(13, 0, 2, 3) = -F.transpose() * Kme.block(10, 5, 3, 3); Kme.block(13, 0, 2, 3) = -F.transpose() * Kme.block(10, 10, 3, 3);
			//Kme.block(3, 3, 2, 2) = F.transpose() * Kme.block(0, 0, 3, 3) * F; Kme.block(3, 8, 2, 2) = F.transpose() * Kme.block(0, 5, 3, 3) * F; Kme.block(3, 13, 2, 2) = F.transpose() * Kme.block(0, 10, 3, 3) * F;
			//Kme.block(8, 3, 2, 2) = F.transpose() * Kme.block(5, 0, 3, 3) * F; Kme.block(8, 8, 2, 2) = F.transpose() * Kme.block(5, 5, 3, 3) * F; Kme.block(8, 13, 2, 2) = F.transpose() * Kme.block(5, 10, 3, 3) * F;
			//Kme.block(13, 3, 2, 2) = F.transpose() * Kme.block(10, 0, 3, 3) * F; Kme.block(13, 8, 2, 2) = F.transpose() * Kme.block(10, 5, 3, 3) * F; Kme.block(13, 13, 2, 2) = F.transpose() * Kme.block(10, 10, 3, 3) * F;

			int ja = 0; int jA = 3; int jb = 5; int jB = 8; int jc = 10; int jC = 13;
			Matrix3d Kmaa = Kme.block(ja, ja, 3, 3); Matrix3d Kmab = Kme.block(ja, jb, 3, 3); Matrix3d Kmac = Kme.block(ja, jc, 3, 3);
			Matrix3d Kmba = Kme.block(jb, ja, 3, 3); Matrix3d Kmbb = Kme.block(jb, jb, 3, 3); Matrix3d Kmbc = Kme.block(jb, jc, 3, 3);
			Matrix3d Kmca = Kme.block(jc, ja, 3, 3); Matrix3d Kmcb = Kme.block(jc, jb, 3, 3); Matrix3d Kmcc = Kme.block(jc, jc, 3, 3);
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
			f.segment<3>(mesh.faces[i]->v[0]->node->index * 3) += (fie.segment<3>(0) + fme.segment<3>(0));
			if (mesh.faces[i]->v[0]->node->eol_index != -1) f.segment<2>(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2) += -F.transpose() * (fie.segment<3>(0) + fme.segment<3>(0));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[0]->node->index * 3 + k, Maa(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[0]->node->index * 3 + k, MDKaa(j, k)));
				}
			}

			Matrix3d Mbb = Mie.block(5, 5, 3, 3);
			Matrix3d Kbb = Kme.block(5, 5, 3, 3);
			Matrix3d MDKbb = Mbb + damping(0) * h * Mbb + damping(1) * h * h * Kbb;
			f.segment<3>(mesh.faces[i]->v[1]->node->index * 3) += (fie.segment<3>(5) + fme.segment<3>(5));
			if (mesh.faces[i]->v[1]->node->eol_index != -1)f.segment<2>(mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2) += -F.transpose() * (fie.segment<3>(5) + fme.segment<3>(5));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + j, mesh.faces[i]->v[1]->node->index * 3 + k, Mbb(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + j, mesh.faces[i]->v[1]->node->index * 3 + k, MDKbb(j, k)));
				}
			}

			Matrix3d Mcc = Mie.block(10, 10, 3, 3);
			Matrix3d Kcc = Kme.block(10, 10, 3, 3);
			Matrix3d MDKcc = Mcc + damping(0) * h * Mcc + damping(1) * h * h * Kcc;
			f.segment<3>(mesh.faces[i]->v[2]->node->index * 3) += (fie.segment<3>(10) + fme.segment<3>(10));
			if (mesh.faces[i]->v[2]->node->eol_index != -1) f.segment<2>(mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2) += -F.transpose() * (fie.segment<3>(10) + fme.segment<3>(10));
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, Mcc(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, MDKcc(j, k)));
				}
			}

			Matrix3d Mab = Mie.block(0, 5, 3, 3);
			Matrix3d Kab = Kme.block(0, 5, 3, 3);
			Matrix3d MDKab = Mab + damping(0) * h * Mab + damping(1) * h * h * Kab;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[1]->node->index * 3 + k, Mab(j, k)));
					M_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + k, mesh.faces[i]->v[0]->node->index * 3 + j, Mab(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[1]->node->index * 3 + k, MDKab(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + k, mesh.faces[i]->v[0]->node->index * 3 + j, MDKab(j, k)));
				}
			}

			Matrix3d Mac = Mie.block(0, 10, 3, 3);
			Matrix3d Kac = Kme.block(0, 10, 3, 3);
			Matrix3d MDKac = Mac + damping(0) * h * Mac + damping(1) * h * h * Kac;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, Mac(j, k)));
					M_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + k, mesh.faces[i]->v[0]->node->index * 3 + j, Mac(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[0]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, MDKac(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + k, mesh.faces[i]->v[0]->node->index * 3 + j, MDKac(j, k)));
				}
			}

			Matrix3d Mbc = Mie.block(5, 10, 3, 3);
			Matrix3d Kbc = Kme.block(5, 10, 3, 3);
			Matrix3d MDKbc = Mbc + damping(0) * h * Mbc + damping(1) * h * h * Kbc;
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					M_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, Mbc(j, k)));
					M_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + k, mesh.faces[i]->v[1]->node->index * 3 + j, Mbc(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[1]->node->index * 3 + j, mesh.faces[i]->v[2]->node->index * 3 + k, MDKbc(j, k)));
					MDK_.push_back(T(mesh.faces[i]->v[2]->node->index * 3 + k, mesh.faces[i]->v[1]->node->index * 3 + j, MDKbc(j, k)));
				}
			}

			// X values
			if (mesh.faces[i]->v[0]->node->eol_index != -1) {
				Matrix2d MAA = Mie.block(3, 3, 2, 2);
				Matrix2d KAA = Kme.block(3, 3, 2, 2);
				Matrix2d MDKAA = MAA + damping(0) * h * MAA + damping(1) * h * h * KAA;
				//f.segment<2>(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2) += (fie.segment<2>(3) + fme.segment<2>(3));
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + k, MAA(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + k, MDKAA(j, k)));
					}
				}
			}

			if (mesh.faces[i]->v[1]->node->eol_index != -1) {
				Matrix2d MBB = Mie.block(8, 8, 2, 2);
				Matrix2d KBB = Kme.block(8, 8, 2, 2);
				Matrix2d MDKBB = MBB + damping(0) * h * MBB + damping(1) * h * h * KBB;
				//f.segment<2>(mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2) += (fie.segment<2>(8) + fme.segment<2>(8));
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + k, MBB(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + k, MDKBB(j, k)));
					}
				}
			}

			if (mesh.faces[i]->v[2]->node->eol_index != -1) {
				Matrix2d MCC = Mie.block(13, 13, 2, 2);
				Matrix2d KCC = Kme.block(13, 13, 2, 2);
				Matrix2d MDKCC = MCC + damping(0) * h * MCC + damping(1) * h * h * KCC;
				//f.segment<2>(mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2) += (fie.segment<2>(13) + fme.segment<2>(13));
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, MCC(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, MDKCC(j, k)));
					}
				}
			}

			if (mesh.faces[i]->v[0]->node->eol_index != -1 && mesh.faces[i]->v[1]->node->eol_index != -1) {
				Matrix2d MAB = Mie.block(3, 8, 2, 2);
				Matrix2d KAB = Kme.block(3, 8, 2, 2);
				Matrix2d MDKAB = MAB + damping(0) * h * MAB + damping(1) * h * h * KAB;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + k, MAB(j, k)));
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + k, mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, MAB(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + k, MDKAB(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + k, mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, MDKAB(j, k)));
					}
				}
			}

			if (mesh.faces[i]->v[0]->node->eol_index != -1 && mesh.faces[i]->v[2]->node->eol_index != -1) {
				Matrix2d MAC = Mie.block(3, 13, 2, 2);
				Matrix2d KAC = Kme.block(3, 13, 2, 2);
				Matrix2d MDKAC = MAC + damping(0) * h * MAC + damping(1) * h * h * KAC;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, MAC(j, k)));
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, MAC(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, MDKAC(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, mesh.nodes.size() * 3 + mesh.faces[i]->v[0]->node->eol_index * 2 + j, MDKAC(j, k)));
					}
				}
			}

			if (mesh.faces[i]->v[1]->node->eol_index != -1 && mesh.faces[i]->v[2]->node->eol_index != -1) {
				Matrix2d MBC = Mie.block(8, 13, 2, 2);
				Matrix2d KBC = Kme.block(8, 13, 2, 2);
				Matrix2d MDKBC = MBC + damping(0) * h * MBC + damping(1) * h * h * KBC;
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, MBC(j, k)));
						M_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + j, MBC(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + j, mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, MDKBC(j, k)));
						MDK_.push_back(T(mesh.nodes.size() * 3 + mesh.faces[i]->v[2]->node->eol_index * 2 + k, mesh.nodes.size() * 3 + mesh.faces[i]->v[1]->node->eol_index * 2 + j, MDKBC(j, k)));
					}
				}
			}

			// x-X values
			if (mesh.faces[i]->v[0]->node->eol_index != -1) {
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

			if (mesh.faces[i]->v[1]->node->eol_index != -1) {
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

			if (mesh.faces[i]->v[2]->node->eol_index != -1) {
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
		// Lagrangian faces
		else {
			double fi[9]; // Gravity force vector
			double Mi[81]; //Inertial matrix
			double fm[9];
			double Km[81];
			ComputeMembraneLAG(xa, xb, xc, Xa, Xb, Xc, e, nu, PP, QQ, Wm, fm, Km);
			VectorXd fme = Map<VectorXd>(fm, 9);
			MatrixXd Kme = Map<MatrixXd>(Km, 9, 9);
			ComputeInertiaLAG(xa, xb, xc, Xa, Xb, Xc, g, material.density, Wi, fi, Mi);
			VectorXd fie = Map<VectorXd>(fi, 9);
			MatrixXd Mie = Map<MatrixXd>(Mi, 9, 9);

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
	}

	for (int i = 0; i < mesh.edges.size(); i++) {

		if (mesh.edges[i]->preserve) eol_edges++;

		//if (is_seam_or_boundary(mesh.edges[i])) {
		//	double added_weight = (material.edge_density * edge_length(mesh.edges[i])) / 2;
		//	for (int j = 0; j < 3; j++) {
		//		M_.push_back(T(mesh.edges[i]->n[0]->index * 3 + j, mesh.edges[i]->n[0]->index * 3 + j, added_weight));
		//		M_.push_back(T(mesh.edges[i]->n[1]->index * 3 + j, mesh.edges[i]->n[1]->index * 3 + j, added_weight));
		//	}
		//}

		if (mesh.edges[i]->adjf[0] == NULL || mesh.edges[i]->adjf[1] == NULL) {
			continue;
		}

		// A check for EoL calc
		bool to_eolA = mesh.edges[i]->n[0]->EoL;
		bool to_eolB = mesh.edges[i]->n[1]->EoL;
		bool to_eolC = false;
		bool to_eolD = false;

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
		int aindexX;
		int bindexX;
		int cindexX;
		int dindexX;

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
		aindexX = mesh.nodes.size() * 3 + mesh.edges[i]->n[0]->eol_index * 2;
		txb << mesh.edges[i]->n[1]->x[0], mesh.edges[i]->n[1]->x[1], mesh.edges[i]->n[1]->x[2];
		bindex = mesh.edges[i]->n[1]->index * 3;
		bindexX = mesh.nodes.size() * 3 + mesh.edges[i]->n[1]->eol_index * 2;
		for (int j = 0; j < 3; j++) {
			if (mesh.edges[i]->adjf[0]->v[j]->node->x != mesh.edges[i]->n[0]->x && mesh.edges[i]->adjf[0]->v[j]->node->x != mesh.edges[i]->n[1]->x) {
				txc << mesh.edges[i]->adjf[0]->v[j]->node->x[0], mesh.edges[i]->adjf[0]->v[j]->node->x[1], mesh.edges[i]->adjf[0]->v[j]->node->x[2];
				cindex = mesh.edges[i]->adjf[0]->v[j]->node->index * 3;
				cindexX = mesh.nodes.size() * 3 + mesh.edges[i]->adjf[0]->v[j]->node->eol_index * 2;
				if (mesh.edges[i]->adjf[0]->v[j]->node->EoL) to_eolC = true;
			}
		}
		for (int j = 0; j < 3; j++) {
			if (mesh.edges[i]->adjf[1]->v[j]->node->x != mesh.edges[i]->n[0]->x && mesh.edges[i]->adjf[1]->v[j]->node->x != mesh.edges[i]->n[1]->x) {
				txd << mesh.edges[i]->adjf[1]->v[j]->node->x[0], mesh.edges[i]->adjf[1]->v[j]->node->x[1], mesh.edges[i]->adjf[1]->v[j]->node->x[2];
				dindex = mesh.edges[i]->adjf[1]->v[j]->node->index * 3;
				dindexX = mesh.nodes.size() * 3 + mesh.edges[i]->adjf[1]->v[j]->node->eol_index * 2;
				if (mesh.edges[i]->adjf[1]->v[j]->node->EoL) to_eolD = true;
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
		Map<Vector3d>(xa, 3) = txa;
		Map<Vector3d>(xb, 3) = txb;
		Map<Vector3d>(xc, 3) = txc;
		Map<Vector3d>(xd, 3) = txd;
		Map<Vector2d>(Xa, 2) = tXa;
		Map<Vector2d>(Xb, 2) = tXb;
		Map<Vector2d>(Xc, 2) = tXc;
		Map<Vector2d>(Xd, 2) = tXd;

		if (to_eolA || to_eolB || to_eolC || to_eolD) {
			MatrixXd F1 = deform_grad(mesh.edges[i]->adjf[0]);
			MatrixXd F2 = deform_grad(mesh.edges[i]->adjf[1]);
			MatrixXd F = (F1 + F2) / 2;
			double fb[20]; // Bending force vector
			double Kb[400]; //Bending stiffness matrix
			ComputeBendingEOL(xa, xb, xc, xd, Xa, Xb, Xc, Xd, beta, Wb, fb, Kb);
			VectorXd fbe = Map<VectorXd>(fb, 20);
			MatrixXd Kbe = Map<MatrixXd>(Kb, 20, 20);

			int ja = 0; int jA = 3; int jb = 5; int jB = 8; int jc = 10; int jC = 13; int jd = 15; int jD = 18;
			Matrix3d Kbaa = Kbe.block(ja, ja, 3, 3); Matrix3d Kbab = Kbe.block(ja, jb, 3, 3); Matrix3d Kbac = Kbe.block(ja, jc, 3, 3); Matrix3d Kbad = Kbe.block(ja, jd, 3, 3);
			Matrix3d Kbba = Kbe.block(jb, ja, 3, 3); Matrix3d Kbbb = Kbe.block(jb, jb, 3, 3); Matrix3d Kbbc = Kbe.block(jb, jc, 3, 3); Matrix3d Kbbd = Kbe.block(jb, jd, 3, 3);
			Matrix3d Kbca = Kbe.block(jc, ja, 3, 3); Matrix3d Kbcb = Kbe.block(jc, jb, 3, 3); Matrix3d Kbcc = Kbe.block(jc, jc, 3, 3); Matrix3d Kbcd = Kbe.block(jc, jd, 3, 3);
			Matrix3d Kbda = Kbe.block(jd, ja, 3, 3); Matrix3d Kbdb = Kbe.block(jd, jb, 3, 3); Matrix3d Kbdc = Kbe.block(jd, jc, 3, 3); Matrix3d Kbdd = Kbe.block(jd, jd, 3, 3);
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
				//f.segment<2>(aindexX) += fbe.segment<2>(3);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(aindexX + j, aindexX + k, K00X(j, k)));
					}
				}
			}

			if (to_eolB) {
				Matrix2d K11X = damping(1) * h * h * Kbe.block(8, 8, 2, 2);
				//f.segment<2>(bindexX) += fbe.segment<2>(8);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(bindexX + j, bindexX + k, K11X(j, k)));
					}
				}
			}

			if (to_eolC) {
				Matrix2d K22X = damping(1) * h * h * Kbe.block(13, 13, 2, 2);
				//f.segment<2>(cindexX) += fbe.segment<2>(13);
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						MDK_.push_back(T(cindexX + j, cindexX + k, K22X(j, k)));
					}
				}
			}

			if (to_eolD) {
				Matrix2d K33X = damping(1) * h * h * Kbe.block(18, 18, 2, 2);
				//f.segment<2>(dindexX) += fbe.segment<2>(18);
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
			double fb[12]; // Bending force vector
			double Kb[144]; //Bending stiffness matrix
			ComputeBendingLAG(xa, xb, xc, xd, Xa, Xb, Xc, Xd, beta, Wb, fb, Kb);
			VectorXd fbe = Map<VectorXd>(fb, 12);
			MatrixXd Kbe = Map<MatrixXd>(Kb, 12, 12);

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
	ClothTimer[2]->toc();

	ClothTimer[3]->tic();
	if (collision == 0 && false) {
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

			// Fix all X points
			for (int i = mesh.nodes.size() * 3; i < mesh.nodes.size() * 3 + eols * 2; i++) {
				MDK_.push_back(T(mesh.nodes.size() * 3 + i, mesh.nodes.size() * 3 + eols * 2 + i, 1));
				MDK_.push_back(T(mesh.nodes.size() * 3 + i + 1, mesh.nodes.size() * 3 + eols * 2 + i + 1, 1));
				MDK_.push_back(T(mesh.nodes.size() * 3 + eols * 2 + i, mesh.nodes.size() * 3 + i));
				MDK_.push_back(T(mesh.nodes.size() * 3 + eols * 2 + i + 1, mesh.nodes.size() * 3 + i + 1));
			}

			Eigen::SparseMatrix<double> M(sparse_size, sparse_size);
			Eigen::SparseMatrix<double> MDK(sparse_size + (holdcount * 3), sparse_size + (holdcount * 3));
			VectorXd b;
			VectorXd KKT_b;

			M.setFromTriplets(M_.begin(), M_.end());
			MDK.setFromTriplets(MDK_.begin(), MDK_.end());
			b = (M*v + h*f);

			VectorXd append_zeros = VectorXd::Zero(eols * 2 + holdcount * 3);
			int append_count = 0;
			if (step > pull_step) {
				for (int j = 0; j < fixed.rows(); j++) {
					if (fixed(j, 0) != -1) {
						if (j >= 4) {
							append_zeros(eols * 2 + append_count * 3 + 0) = fixed(j, 3);
							append_zeros(eols * 2 + append_count * 3 + 1) = fixed(j, 4);
							append_zeros(eols * 2 + append_count * 3 + 2) = fixed(j, 5);
							append_count++;
						}
					}
				}
			}

			// Append zeros for X
			//VectorXd append_zeros = VectorXd::Zero(eols * 2);

			//KKT_b.resize(b.size() + (holdcount * 3));
			KKT_b.resize(b.size() + (eols * 2) + (holdcount * 3));
			KKT_b << b, append_zeros;

			LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
			lscg.compute(MDK);
			v = lscg.solve(KKT_b);
		}
	}
	else {
		//cout << "Quadratic Solver" << endl;
		Eigen::SparseMatrix<double> M(sparse_size + eols * 2, sparse_size + eols * 2);
		Eigen::SparseMatrix<double> MDK(sparse_size + eols * 2, sparse_size + eols * 2);
		VectorXd b;

		M.setFromTriplets(M_.begin(), M_.end());
		MDK.setFromTriplets(MDK_.begin(), MDK_.end());
		//b = -(M*v + h*f);

		//Eigen::SparseMatrix<double> N(collision + LAG_constraints + EOL_constraints, sparse_size + eols * 2);
		//N.setFromTriplets(N_.begin(), N_.end());
		Eigen::SparseMatrix<double> Aineq(collision + aineqsize, sparse_size + eols * 2);
		//Eigen::SparseMatrix<double> Aeq(aeqsize + EOL_constraints, sparse_size + eols * 2);
		Aineq.setFromTriplets(Aineq_.begin(), Aineq_.end());
		//Aeq.setFromTriplets(Aeq_.begin(), Aeq_.end());
		//igl::mosek::MosekData mosek_data;
		//auto lc = Eigen::VectorXd(collision + LAG_constraints + EOL_constraints);
		//auto uc = Eigen::VectorXd(collision + LAG_constraints + EOL_constraints);
		//auto lx = Eigen::VectorXd(sparse_size + eols * 2);
		//auto ux = Eigen::VectorXd(sparse_size + eols * 2);

		// Constrain EOL boundaries
		//for (int i = sparse_size; i < (sparse_size + eols * 2); i++) {
		//	lx(i) = -MSK_INFINITY;
		//	ux(i) = +MSK_INFINITY;
		//}
		for (int i = 0; i < mesh.nodes.size(); i++) {
			if (mesh.nodes[i]->EoL && is_seam_or_boundary(mesh.nodes[i])) {
				if (mesh.nodes[i]->verts[0]->u[0] == Xmin || mesh.nodes[i]->verts[0]->u[0] == Xmax) {
					//lx(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) = 0.0;
					//ux(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) = 0.0;
					Aeq_.push_back(T(aeqsize + EOL_constraints, mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2, 1));
					aeqsize++;
				}
				else {
					//lx(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) = 0.0;
					//ux(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) = 0.0;
					Aeq_.push_back(T(aeqsize + EOL_constraints, mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1, 1));
					aeqsize++;
				}
			}
			//lx(i * 3) = -MSK_INFINITY;
			//ux(i * 3) = +MSK_INFINITY;
			//lx(i * 3 + 1) = -MSK_INFINITY;
			//ux(i * 3 + 1) = +MSK_INFINITY;
			//lx(i * 3 + 2) = -MSK_INFINITY;
			//ux(i * 3 + 2) = +MSK_INFINITY;
		}

		// TODO:: Better way
		int append_count = 0;
		vector<double> beq_;
		double expofilA = 0.01;
		if (pull_free) {
			for (int j = 0; j < fixed.rows(); j++) {
				if (fixed(j, 0) != -1) {
					if (j == 4) {
						//lx(0) = fixed(j, 3);
						//ux(0) = fixed(j, 3);
						//lx(1) = fixed(j, 4);
						//ux(1) = fixed(j, 4);
						//lx(2) = fixed(j, 5);
						//ux(2) = fixed(j, 5);
						if (step > pull_step) {
							if (fixed(j, 0) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[0]->v[0] + expofilA * fixed(j, 3));
							if (fixed(j, 1) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[0]->v[1] + expofilA * fixed(j, 4));
							if (fixed(j, 2) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[0]->v[2] + expofilA * fixed(j, 5));
						}
						if (fixed(j, 0) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, 0, fixed(j, 0)));
							aeqsize++;
							append_count++;
						}
						if (fixed(j, 1) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, 1, fixed(j, 1)));
							aeqsize++;
							append_count++;
						}
						if (fixed(j, 2) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, 2, fixed(j, 2)));
							aeqsize++;
							append_count++;
						}
						//append_count++;
					}
					else if (j == 5) {
						//lx((rows * (cols - 1)) * 3) = fixed(j, 3);
						//ux((rows * (cols - 1)) * 3) = fixed(j, 3);
						//lx((rows * (cols - 1)) * 3 + 1) = fixed(j, 4);
						//ux((rows * (cols - 1)) * 3 + 1) = fixed(j, 4);
						//lx((rows * (cols - 1)) * 3 + 2) = fixed(j, 5);
						//ux((rows * (cols - 1)) * 3 + 2) = fixed(j, 5);
						if (step > pull_step) {
							if (fixed(j, 0) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows * (cols - 1)]->v[0] + expofilA * fixed(j, 3));
							if (fixed(j, 1) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows * (cols - 1)]->v[1] + expofilA * fixed(j, 4));
							if (fixed(j, 2) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows * (cols - 1)]->v[2] + expofilA * fixed(j, 5));
						}
						if (fixed(j, 0) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows * (cols - 1)) * 3, fixed(j, 0)));
							aeqsize++;
							append_count++;
						}
						if (fixed(j, 1) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows * (cols - 1)) * 3 + 1, fixed(j, 1)));
							aeqsize++;
							append_count++;
						}
						if (fixed(j, 2) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows * (cols - 1)) * 3 + 2, fixed(j, 2)));
							aeqsize++;
							append_count++;
						}
						//append_count++;
					}
					else if (j == 6) {
						//lx((rows * cols - 1) * 3) = fixed(j, 3);
						//ux((rows * cols - 1) * 3) = fixed(j, 3);
						//lx((rows * cols - 1) * 3 + 1) = fixed(j, 4);
						//ux((rows * cols - 1) * 3 + 1) = fixed(j, 4);
						//lx((rows * cols - 1) * 3 + 2) = fixed(j, 5);
						//ux((rows * cols - 1) * 3 + 2) = fixed(j, 5);
						if (step > pull_step) {
							if (fixed(j, 0) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows * cols - 1]->v[0] + expofilA * fixed(j, 3));
							if (fixed(j, 1) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows * cols - 1]->v[1] + expofilA * fixed(j, 4));
							if (fixed(j, 2) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows * cols - 1]->v[2] + expofilA * fixed(j, 5));
						}
						if (fixed(j, 0) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows * cols - 1) * 3, fixed(j, 0)));
							aeqsize++;
							append_count++;
						}
						if (fixed(j, 1) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows * cols - 1) * 3 + 1, fixed(j, 1)));
							aeqsize++;
							append_count++;
						}
						if (fixed(j, 2) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows * cols - 1) * 3 + 2, fixed(j, 2)));
							aeqsize++;
							append_count++;
						}
						//append_count++;
					}
					else if (j == 7) {
						//lx((rows - 1) * 3) = fixed(j, 3);
						//ux((rows - 1) * 3) = fixed(j, 3);
						//lx((rows - 1) * 3 + 1) = fixed(j, 4);
						//ux((rows - 1) * 3 + 1) = fixed(j, 4);
						//lx((rows - 1) * 3 + 2) = fixed(j, 5);
						//ux((rows - 1) * 3 + 2) = fixed(j, 5);
						if (step > pull_step) {
							if (fixed(j, 0) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows - 1]->v[0] + expofilA * fixed(j, 3));
							if (fixed(j, 1) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows - 1]->v[1] + expofilA * fixed(j, 4));
							if (fixed(j, 2) == 1.0) beq_.push_back((1 - expofilA) * mesh.nodes[rows - 1]->v[2] + expofilA * fixed(j, 5));
						}
						if (fixed(j, 0) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows - 1) * 3, fixed(j, 0)));
							aeqsize++;
							append_count++;
						}
						if (fixed(j, 1) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows - 1) * 3 + 1, fixed(j, 1)));
							aeqsize++;
							append_count++;
						}
						if (fixed(j, 2) == 1.0) {
							Aeq_.push_back(T(aeqsize + EOL_constraints, (rows - 1) * 3 + 2, fixed(j, 2)));
							aeqsize++;
							append_count++;
						}
						//append_count++;
					}
				}
			}
		}

		//if (pull_free) {
		//	// Pull sides
		//	for (int i = 0; i < mesh.nodes.size(); i++) {
		//		if (is_seam_or_boundary(mesh.nodes[i])) {
		//			if (mesh.nodes[i]->preserve) continue;
		//			if (mesh.nodes[i]->verts[0]->u[0] == 0.0) {
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3, 1));
		//				aeqsize++;
		//				append_count++;
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 1, 1));
		//				aeqsize++;
		//				append_count++;
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 2, 1));
		//				aeqsize++;
		//				append_count++;
		//				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * fixed(4, 3));
		//				if (step < side_pull_change_step) {
		//					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * 0);
		//				}
		//				else {
		//					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * fixed(4,4));
		//				}
		//				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[2] + expofilA * fixed(4, 5));
		//			}
		//			if (mesh.nodes[i]->verts[0]->u[0] == 1.0) {
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3, 1));
		//				aeqsize++;
		//				append_count++;
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 1, 1));
		//				aeqsize++;
		//				append_count++;
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 2, 1));
		//				aeqsize++;
		//				append_count++;
		//				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * fixed(6, 3));
		//				if (step < side_pull_change_step) {
		//					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * 0);
		//				}
		//				else {
		//					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * fixed(4, 4));
		//				}
		//				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[2] + expofilA * fixed(4, 5));
		//			}
		//			if (mesh.nodes[i]->verts[0]->u[1] == 0.0) {
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3, 1));
		//				aeqsize++;
		//				append_count++;
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 1, 1));
		//				aeqsize++;
		//				append_count++;
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 2, 1));
		//				aeqsize++;
		//				append_count++;
		//				if (step < side_pull_change_step) {
		//					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * 0);
		//				}
		//				else {
		//					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * fixed(4, 3));
		//				}
		//				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * fixed(4, 4));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[2] + expofilA * fixed(4, 5));
		//			}
		//			if (mesh.nodes[i]->verts[0]->u[1] == 1.0) {
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3, 1));
		//				aeqsize++;
		//				append_count++;
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 1, 1));
		//				aeqsize++;
		//				append_count++;
		//				Aeq_.push_back(T(aeqsize + EOL_constraints, i * 3 + 2, 1));
		//				aeqsize++;
		//				append_count++;
		//				if (step < side_pull_change_step) {
		//					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * 0);
		//				}
		//				else {
		//					beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[0] + expofilA * fixed(5, 3));
		//				}
		//				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[1] + expofilA * fixed(5, 4));
		//				beq_.push_back((1 - expofilA) * mesh.nodes[i]->v[2] + expofilA * fixed(5, 5));
		//			}
		//		}
		//	}
		//}

		Eigen::SparseMatrix<double> Aeq(aeqsize + EOL_constraints, sparse_size + eols * 2);
		Aeq.setFromTriplets(Aeq_.begin(), Aeq_.end());
		VectorXd bineq(collision + aineqsize);
		VectorXd beq(aeqsize + EOL_constraints);
		bineq.setZero();
		beq.setZero();

		for (int i = 0; i < ineq_move_index.size(); i++) {
			bineq(ineq_move_index[i]) = ineq_move[i];
		}

		for (int i = 0; i < eq_move_index.size(); i++) {
			beq(eq_move_index[i]) = eq_move[i];
		}

		if (step > pull_step) {
			for (int i = 0; i < beq_.size(); i++) {
				beq((aeqsize + EOL_constraints) - (append_count) + i) = beq_[i];
			}
		}
		//cout << beq << endl;

		// Base collision inequality constraint 
		//for (int i = 0; i < collision; i++) {
		//	//lc(i) = -(0.1/h)* 0.001; // BAUMGARTE
		//	lc(i) = 0.0;
		//	uc(i) = +MSK_INFINITY;
		//}
		// Eqaulity constraint for further lagrangian contraining
		//for (int i = collision; i < collision + eols; i++) {
		//	//lc(i) = -(0.1/h)* 0.001; // BAUMGARTE
		//	lc(i) = 0.0;
		//	uc(i) = 0.0;
		//}
		//for (int i = 0; i < collision + LAG_constraints + EOL_constraints; i++) {
		//	lc(i) = fill_lc[i];
		//	uc(i) = fill_uc[i];
		//}

		// Velocity transfer
		// Regularizing
		//MatrixXd regulize = MatrixXd::Zero(M.rows(), M.cols());
		//regulize.block(mesh.nodes.size() * 3, mesh.nodes.size() * 3, eols * 2, eols * 2) = MatrixXd::Identity(eols * 2, eols * 2);
		//igl::mosek::mosek_quadprog(M+(1e-8)*regulize, -M*v_old, 0, N, lc, uc, lx, ux, mosek_data, v);

		b = -(M*v_old + h*f);

		if (matlab_debug_physics) {
			//sparse_to_file_as_dense(MDK, "MDK");
			ofstream ofs;
			ofs.open("test.m", ofstream::out | ofstream::trunc);
			ofs.close();
			double_to_file(h, "h");
			double_to_file(grav(2), "grav");
			double_to_file(material.density, "rho");
			double_to_file(e, "e");
			double_to_file(nu, "nu");
			mat_to_file(x_X, "x_X");
			vec_to_file(isEoL, "isEol");
			VectorXi vvv(3);
			vvv << 1, 1, 1;
			mat_to_file(faces2.colwise() += vvv, "faces");
			sparse_to_file_as_sparse_m(MDK, "MDK");
			vec_to_file(b, "b");
			//sparse_to_file_as_dense(N, "N");
			//sparse_to_file_as_dense(M, "M");
			//sparse_to_file_as_sparse_m(N, "N");
			sparse_to_file_as_sparse_m(M, "M");
			vec_to_file(f, "f");
			//vec_to_file(lc, "lc");
			//vec_to_file(uc, "uc");
			//vec_to_file(lx, "lx");
			//vec_to_file(ux, "ux");
			vec_to_file(v_old, "v_last_step");
			//vec_to_file(v, "v_transfer");
			sparse_to_file_as_sparse_m(Aineq, "Aineq");
			sparse_to_file_as_sparse_m(Aeq, "Aeq");
			vec_to_file(bineq, "bineq");
			vec_to_file(beq, "beq");
		}

		igl::mosek::mosek_quadprog(MDK, b, 0, N, lc, uc, lx, ux, mosek_data, v);

		QuadProgMosek *program = new QuadProgMosek();
		double inf = std::numeric_limits<double>::infinity();
		VectorXd xl;
		VectorXd xu;
		xl.setConstant(sparse_size + eols * 2, -inf);
		xu.setConstant(sparse_size + eols * 2, inf);
		program->setNumberOfVariables(sparse_size + eols * 2);
		program->setNumberOfInequalities(collision + aineqsize);
		program->setNumberOfEqualities(aeqsize + EOL_constraints);
		program->setObjectiveMatrix(MDK);
		program->setObjectiveVector(b);
		program->setInequalityMatrix(Aineq);
		program->setInequalityVector(bineq);
		program->setEqualityMatrix(Aeq);
		program->setEqualityVector(beq);
		//program->setLowerVariableBound(xl);
		//program->setUpperVariableBound(xu);
		program->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
		program->setParamInt(MSK_IPAR_LOG, 200);
		program->setParamInt(MSK_IPAR_LOG_INTPNT, 200);
		program->setParamInt(MSK_IPAR_LOG_MIO, 200);
		program->setParamInt(MSK_IPAR_LOG_CUT_SECOND_OPT, 200);
		program->setParamInt(MSK_IPAR_LOG_SIM, 200);
		program->setParamInt(MSK_IPAR_LOG_RESPONSE, 200);
		// https://docs.mosek.com/8.1/capi/parameters.html#mosek.dparam.intpnt_qo_tol_dfeas
		program->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-11);
		program->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-11);
		program->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-11);
		program->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
		program->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-11);
		program->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-11);

		bool success = program->solve();

		if ((success) {
			v = program->getPrimalSolution();
			//cout << "Ineq: " << program->getDualInequality() << endl;
			//cout << "Equal: " << program->getDualEquality() << endl;
		}


		if (matlab_debug_physics) {
		//	sparse_to_file_as_dense(MDK, "MDK");
		//	sparse_to_file_as_dense(b, "b");
		//	sparse_to_file_as_dense(N, "N");
		//	sparse_to_file_as_dense(lc, "lc");
		//	sparse_to_file_as_dense(uc, "uc");
		//	sparse_to_file_as_dense(lx, "lx");
		//	sparse_to_file_as_dense(ux, "ux");
		//	sparse_to_file_as_dense(v, "v");
			vec_to_file(v, "v_solved");
		}
	}
	ClothTimer[3]->toc();

	// Hack alpha
	if (friction) {
		if (alpha > 0.0) {
			alpha -= alpha_change;
		}
		if (alpha < 0.0) {
			alpha = 0.0;
		}
	}

	//ofstream ofs;
	//ofs.open("Z.m", ofstream::out | ofstream::app);
	for (int i = 0; i < mesh.nodes.size(); i++) {
		//if(i == 0) ofs << mesh.nodes[0]->x[2] << " ";
		// Friction
		if (mesh.nodes[i]->coll_case > 0) {
			if (mesh.nodes[i]->coll_case == 1) { // Face
				Vec3 vpre = Vec3(v(i * 3), v(i * 3 + 1), v(i * 3 + 2));
				//Vec3 vnpre = dot(mesh.nodes[i]->nor_ave, vpre) * mesh.nodes[i]->nor_ave;
				//Vec3 vtpre = vpre - vnpre;
				//Vec3 vt = alpha * vtpre;
				mesh.nodes[i]->v = (alpha * (vpre - (dot(mesh.nodes[i]->nor_ave, vpre) * mesh.nodes[i]->nor_ave)));
				mesh.nodes[i]->x = mesh.nodes[i]->x + h * mesh.nodes[i]->v;
			}
			else if (mesh.nodes[i]->coll_case == 2) { // Edge
				Vec3 vlpre = Vec3(v(i * 3), v(i * 3 + 1), v(i * 3 + 2));
				Vec3 vtepre = Vec3(v(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2), v(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1), 0.0);
				Vec3 vnlpre = dot(mesh.nodes[i]->nor_ave, vlpre) * mesh.nodes[i]->nor_ave;
				Vec3 vtlpre = vlpre - vnlpre;
				Vec3 vtl = alpha * vtlpre;
				Vec3 vte = alpha * vtepre;
				mesh.nodes[i]->v = vnlpre + vtl;
				mesh.nodes[i]->verts[0]->v = vte;
				mesh.nodes[i]->x = mesh.nodes[i]->x + h * mesh.nodes[i]->v;
				if (mesh.nodes[i]->verts[0]->u[0] != Xmin && mesh.nodes[i]->verts[0]->u[0] != Xmax) mesh.nodes[i]->verts[0]->u[0] = mesh.nodes[i]->verts[0]->u[0] + h * mesh.nodes[i]->verts[0]->v[0];
				if (mesh.nodes[i]->verts[0]->u[1] != Ymin && mesh.nodes[i]->verts[0]->u[1] != Ymax) mesh.nodes[i]->verts[0]->u[1] = mesh.nodes[i]->verts[0]->u[1] + h * mesh.nodes[i]->verts[0]->v[1];
			}
			else if (mesh.nodes[i]->coll_case == 3) { // Corner
				mesh.nodes[i]->v[0] = v(i * 3);
				mesh.nodes[i]->v[1] = v(i * 3 + 1);
				mesh.nodes[i]->v[2] = v(i * 3 + 2);
				mesh.nodes[i]->x = mesh.nodes[i]->x + h * mesh.nodes[i]->v;
				mesh.nodes[i]->verts[0]->v[0] = alpha * v(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2);
				mesh.nodes[i]->verts[0]->v[1] = alpha * v(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1);
				if (mesh.nodes[i]->verts[0]->u[0] != Xmin && mesh.nodes[i]->verts[0]->u[0] != Xmax) mesh.nodes[i]->verts[0]->u[0] = mesh.nodes[i]->verts[0]->u[0] + h * mesh.nodes[i]->verts[0]->v[0];
				if (mesh.nodes[i]->verts[0]->u[1] != Ymin && mesh.nodes[i]->verts[0]->u[1] != Ymax) mesh.nodes[i]->verts[0]->u[1] = mesh.nodes[i]->verts[0]->u[1] + h * mesh.nodes[i]->verts[0]->v[1];

			}
		}
		else {
			mesh.nodes[i]->v[0] = v(i * 3);
			mesh.nodes[i]->v[1] = v(i * 3 + 1);
			mesh.nodes[i]->v[2] = v(i * 3 + 2);
			mesh.nodes[i]->x = mesh.nodes[i]->x + h * mesh.nodes[i]->v;
			if (mesh.nodes[i]->EoL) {
				//cout << "v: " << mesh.nodes[i]->v << endl;
				//cout << "V: " << mesh.nodes[i]->verts[0]->v << endl;
				//cout << i << " x: " << mesh.nodes[i]->x << endl;
				//cout << i << " X: " << mesh.nodes[i]->verts[0]->u << endl;
				mesh.nodes[i]->verts[0]->v[0] = v(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2) * energyadd;
				mesh.nodes[i]->verts[0]->v[1] = v(mesh.nodes.size() * 3 + mesh.nodes[i]->eol_index * 2 + 1) * energyadd;
				if (mesh.nodes[i]->verts[0]->u[0] != Xmin && mesh.nodes[i]->verts[0]->u[0] != Xmax) mesh.nodes[i]->verts[0]->u[0] = mesh.nodes[i]->verts[0]->u[0] + h * mesh.nodes[i]->verts[0]->v[0];
				if (mesh.nodes[i]->verts[0]->u[1] != Ymin && mesh.nodes[i]->verts[0]->u[1] != Ymax) mesh.nodes[i]->verts[0]->u[1] = mesh.nodes[i]->verts[0]->u[1] + h * mesh.nodes[i]->verts[0]->v[1];
			}
		}
	}

	if (export_timings) {
		ClothTimer[0]->export_csv();
		ClothTimer[1]->export_csv();
		ClothTimer[2]->export_csv();
		ClothTimer[3]->export_csv();
		EOL_outputter(eols, eol_edges,mesh.nodes.size(), mesh.faces.size());
	}

	return true;
}

void Cloth::init()
{
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &posBuf2DID);
	glBindBuffer(GL_ARRAY_BUFFER, posBuf2DID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf2D[0], GL_DYNAMIC_DRAW);
	
	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	
	glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_DYNAMIC_DRAW);
	
	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size()*sizeof(unsigned int), &eleBuf[0], GL_DYNAMIC_DRAW);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	
	assert(glGetError() == GL_NO_ERROR);
}

void Cloth::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p, bool on2D) const
{
	// Draw mesh
	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(p->getUniform("kdBack"),  1, Vector3f(1.0, 1.0, 0.0).data());
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	//for(int i = 0; i < rows; ++i) {
	//	glDrawElements(GL_TRIANGLE_STRIP, 2*cols, GL_UNSIGNED_INT, (const void *)(2*cols*i*sizeof(unsigned int)));
	//}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_DYNAMIC_DRAW);
	glDrawElements(GL_TRIANGLES, eleBuf.size(), GL_UNSIGNED_INT, (void*)0);
	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();

}

// Export
int Cloth::getBrenderCount() const
{
	return 2;
}

vector<string> Cloth::getBrenderNames() const
{
	vector<string> names;
	names.push_back("Cloth2D");
	names.push_back("Cloth3D");
	return names;
}

void Cloth::exportBrender(vector< shared_ptr< ofstream > > outfiles) const
{
	ofstream &outfile = *outfiles[0];

	////vertex positions
	//for (int i = 0; i < posBuf.size(); i = i + 3) {
	//	char vert[50];
	//	//sprintf(vert, "v %f %f %f\n", posBuf[i], posBuf[i + 1], posBuf[i + 2]);
	//	outfile << vert;
	//}
	for (int i = 0; i < mesh.nodes.size(); i++) {
		char vert[50];
		sprintf(vert, "v %f %f %f\n", mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->verts[0]->u[1], 0);
		outfile << vert;
	}
	//texture coordinates
	for (int i = 0; i < texBuf.size(); i = i + 2) {
		char vtex[50];
		sprintf(vtex, "vt %f %f\n", texBuf[i], texBuf[i + 1]);
		outfile << vtex;
	}
	//normal vectors
	for (int i = 0; i < norBuf.size(); i = i + 3) {
		char norm[50];
		sprintf(norm, "vn %f %f %f\n", 0.0, 0.0, 1.0);
		outfile << norm;
	}
	//faces
	for (int i = 0; i < eleBuf.size(); i = i + 3) {
		char face[50];
		sprintf(face, "f %i/%i/%i %i/%i/%i %i/%i/%i\n", eleBuf[i]+1, eleBuf[i]+1, eleBuf[i]+1, eleBuf[i + 1]+1, eleBuf[i + 1]+1, eleBuf[i + 1]+1, eleBuf[i + 2]+1, eleBuf[i + 2]+1, eleBuf[i + 2]+1);
		outfile << face;
	}

	// Vectex position
	//for (int i = 0; i < mesh.nodes.size(); i++) {
	//	char vert[50];
	//	sprintf(vert, "v %f %f %f\n", mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
	//	outfile << vert;
	//}
	////texture coordinates
	//for (int i = 0; i < mesh.verts.size(); i++) {
	//	char vtex[50];
	//	sprintf(vtex, "vt %f %f\n", mesh.verts[i]->u[0], mesh.verts[i]->u[1]);
	//	outfile << vtex;
	//}
	////normal vectors
	//for (int i = 0; i < mesh.nodes.size(); i++) {
	//	char norm[50];
	//	sprintf(norm, "vn %f %f %f\n", mesh.nodes[i]->n[0], mesh.nodes[i]->n[1], mesh.nodes[i]->n[2]);
	//	outfile << norm;
	//}
	////faces
	//for (int i = 0; i < mesh.faces.size(); i++) {
	//	char face[50];
	//	sprintf(face, "f %i/%i/%i %i/%i/%i %i/%i/%i\n", mesh.faces[i]->v[0]->index+1, mesh.faces[i]->v[0]->index + 1, mesh.faces[i]->v[0]->index + 1, mesh.faces[i]->v[1]->index + 1, mesh.faces[i]->v[1]->index + 1, mesh.faces[i]->v[1]->index + 1, mesh.faces[i]->v[2]->index, mesh.faces[i]->v[2]->index + 1, mesh.faces[i]->v[2]->index + 1);
	//	outfile << face;
	//}

	ofstream &outfile2 = *outfiles[1];

	////vertex positions
	for (int i = 0; i < posBuf.size(); i = i + 3) {
		char vert[50];
		sprintf(vert, "v %f %f %f\n", posBuf[i], posBuf[i + 1], posBuf[i + 2]);
		outfile2 << vert;
	}
	//texture coordinates
	for (int i = 0; i < texBuf.size(); i = i + 2) {
		char vtex[50];
		sprintf(vtex, "vt %f %f\n", texBuf[i], texBuf[i + 1]);
		outfile2 << vtex;
	}
	//normal vectors
	for (int i = 0; i < norBuf.size(); i = i + 3) {
		char norm[50];
		sprintf(norm, "vn %f %f %f\n", norBuf[i], norBuf[i + 1], norBuf[i + 2]);
		outfile2 << norm;
	}
	//faces
	for (int i = 0; i < eleBuf.size(); i = i + 3) {
		char face[50];
		sprintf(face, "f %i/%i/%i %i/%i/%i %i/%i/%i\n", eleBuf[i] + 1, eleBuf[i] + 1, eleBuf[i] + 1, eleBuf[i + 1] + 1, eleBuf[i + 1] + 1, eleBuf[i + 1] + 1, eleBuf[i + 2] + 1, eleBuf[i + 2] + 1, eleBuf[i + 2] + 1);
		outfile2 << face;
	}
}