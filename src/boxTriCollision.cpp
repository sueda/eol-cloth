#define _USE_MATH_DEFINES
#include <cmath>
#include <string.h>
#include <iostream>
#include <random>

#include "boxTriCollision.h"

using namespace std;
using namespace Eigen;

#ifdef MATLAB_MEX_FILE
#include "mex.h"
void mexPrintMat(const MatrixXd &A, const char *name = NULL)
{
	if (name) {
		mexPrintf("%s = [\n", name);
	}
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			mexPrintf("% 0.6f ", A(i, j));
		}
		mexPrintf("\n");
	}
	if (name) {
		mexPrintf("];");
	}
	mexPrintf("\n");
}

void mexPrintMati(const MatrixXi &A, const char *name = NULL)
{
	if (name) {
		mexPrintf("%s = [\n", name);
	}
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			mexPrintf("%d ", A(i, j));
		}
		mexPrintf("\n");
	}
	if (name) {
		mexPrintf("];");
	}
	mexPrintf("\n");
}
#endif

void printMat(const MatrixXd &A, const char *name, const char *filename)
{
	FILE *fp = fopen(filename, "a");
	fprintf(fp, "%s = [\n", name);
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			fprintf(fp, "% 0.6f ", A(i, j));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "];");
	fprintf(fp, "\n");
	fclose(fp);
}

void printMati(const MatrixXi &A, const char *name, const char *filename)
{
	FILE *fp = fopen(filename, "a");
	fprintf(fp, "%s = [\n", name);
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			fprintf(fp, "%d ", A(i, j));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "];");
	fprintf(fp, "\n");
	fclose(fp);
}

// From raytri.c by Thomas Akenine-Moller
int intersect_triangle3_inc(
	const double *orig, const double *dir,
	const double *vert0, const double *vert1, const double *vert2,
	double *t, double *u, double *v);

namespace btc
{

	Edge::Edge()
	{
		verts = Vector4i::Zero();
		faces = Vector2i::Zero();
		internal = false;
		angle = 0;
		normals[0] = Vector3d::Zero();
		normals[1] = Vector3d::Zero();
	}

	Edge::~Edge()
	{

	}

	Collision::Collision()
	{
		dist = 0;
		nor1 = Vector3d::Zero();
		nor2 = Vector3d::Zero();
		pos1 = Vector3d::Zero();
		pos2 = Vector3d::Zero();
		count1 = 0;
		count2 = 0;
		verts1 = Vector3i::Zero();
		verts2 = Vector3i::Zero();
		weights1 = Vector3d::Zero();
		weights2 = Vector3d::Zero();
		tri1 = -1;
		tri2 = -1;
		edge2 = -1;
		edgeDir = Vector3d::Zero();
	}

	Collision::~Collision()
	{

	}

	/**
	* Output: an array of structures corresponds to an edge. The vertex ordering
	* is as follows:
	*
	*      x2
	*     /  \
	*    / t0 \
	*   /      \
	* x0--edge--x1
	*   \      /
	*    \ t1 /
	*     \  /
	*      x3
	*/
	void createEdges(
		vector<shared_ptr<Edge> > &edges, // output
		const MatrixXi &faces, // input
		const MatrixXd &verts  // input
	)
	{
		// First, create a list of edges from the triangle list. The third row 
		// stores the "other" vertex. The fourth row stores the triangle index.
		int nf = faces.cols();
		int nv = verts.cols();
		struct FaceEdge
		{
			Vector3i verts;
			int face;
			int hash;
		};
		int n = 3 * nf;
		vector<shared_ptr<FaceEdge> > tmp;
		for (int k = 0; k < nf; ++k) {
			for (int i = 0; i < 3; ++i) {
				auto faceEdge = make_shared<FaceEdge>();
				tmp.push_back(faceEdge);
				faceEdge->verts << (i + 0) % 3, (i + 1) % 3, (i + 2) % 3;
				faceEdge->face = k;
				int v0 = faceEdge->verts(0);
				int v1 = faceEdge->verts(1);
				int kmin = min(faces(v0, k), faces(v1, k)) + 1;
				int kmax = max(faces(v0, k), faces(v1, k)) + 1;
				faceEdge->hash = kmin + (n + 1)*kmax;
			}
		}

		// Sort by hash number
		struct FaceEdgeComp
		{
			bool operator() (const shared_ptr<FaceEdge> &e0, const shared_ptr<FaceEdge> &e1) const
			{
				return (e0->hash < e1->hash);
			}
		};
		FaceEdgeComp comp;
		// We don't need a stable sort, but this will make the ordering the same as
		// Matlab, making it easier to debug.
		stable_sort(tmp.begin(), tmp.end(), comp);

		// Now the twin edges are neighbors in the list
		int k = 0;
		while (k < n) {
			int f0 = tmp[k]->face;
			Vector2i e;
			e(0) = faces(tmp[k]->verts(0), f0);
			e(1) = faces(tmp[k]->verts(1), f0);
			Vector3d xa0 = verts.block<3, 1>(0, faces(0, f0));
			Vector3d xb0 = verts.block<3, 1>(0, faces(1, f0));
			Vector3d xc0 = verts.block<3, 1>(0, faces(2, f0));
			Vector3d n0 = (xb0 - xa0).cross(xc0 - xa0).normalized();
			if (k < n - 1 && tmp[k]->hash == tmp[k + 1]->hash) {
				// k and k+1 are twins
				int f1 = tmp[k + 1]->face;
				Vector3d xa1 = verts.block<3, 1>(0, faces(0, f1));
				Vector3d xb1 = verts.block<3, 1>(0, faces(1, f1));
				Vector3d xc1 = verts.block<3, 1>(0, faces(2, f1));
				Vector3d n1 = (xb1 - xa1).cross(xc1 - xa1).normalized();
				double angle = acos(n0.dot(n1));
				auto edge = make_shared<Edge>();
				edges.push_back(edge);
				edge->verts.segment<2>(0) = e; // edge vertex indices
				edge->verts(2) = faces(tmp[k]->verts(2), f0); // x2 index in figure
				edge->verts(3) = faces(tmp[k + 1]->verts(2), f1); // x3 index in figure
				edge->faces << f0, f1; // t0 and t1 in figure
				edge->internal = true;
				edge->angle = angle; // (TODO: concave edges)
				edge->normals[0] = n0;
				edge->normals[1] = n1;
				k += 2;
			}
			else {
				// k is an external edge (singleton)
				auto edge = make_shared<Edge>();
				edges.push_back(edge);
				edge->verts.segment<2>(0) = e; // edge vertex indices
				edge->verts(2) = faces(tmp[k]->verts(2), f0); // x2 index in figure
				edge->verts(3) = -1; // x3 index in figure
				edge->faces << f0, -1; // t0 and t1 in figure
				edge->internal = false;
				edge->angle = 1e9; // infinite dihedral angle
				edge->normals[0] = n0;
				k += 1;
			}
		}
	}

	MatrixXd createFaceNormals(const MatrixXi &faces, const MatrixXd &verts)
	{
		int nf = faces.cols();
		MatrixXd normals = MatrixXd::Zero(3, nf);
		for (int k = 0; k < nf; ++k) {
			const Vector3i &f = faces.col(k);
			const Vector3d &xa = verts.block<3, 1>(0, f(0));
			const Vector3d &xb = verts.block<3, 1>(0, f(1));
			const Vector3d &xc = verts.block<3, 1>(0, f(2));
			Vector3d dba = xb - xa;
			Vector3d dac = xa - xc;
			normals.col(k) = dba.cross(-dac).normalized();
		}
		return normals;
	}

	MatrixXd createVertNormals(const MatrixXi &faces, const MatrixXd &verts)
	{
		int nv = verts.cols();
		int nf = faces.cols();
		MatrixXd normals = MatrixXd::Zero(3, nv);
		vector<double> angles(nv, 0); // Running sum of angles at the vertex
		for (int k = 0; k < nf; ++k) {
			const Vector3i &f = faces.col(k);
			const Vector3d &xa = verts.block<3, 1>(0, f(0));
			const Vector3d &xb = verts.block<3, 1>(0, f(1));
			const Vector3d &xc = verts.block<3, 1>(0, f(2));
			Vector3d dba = xb - xa;
			Vector3d dcb = xc - xb;
			Vector3d dac = xa - xc;
			Vector3d nor = dba.cross(-dac);
			nor.normalize();
			dba.normalize();
			dcb.normalize();
			dac.normalize();
			double angle1 = acos(dba.dot(-dac));
			double angle2 = acos(dcb.dot(-dba));
			double angle3 = acos(dac.dot(-dcb));
			normals.col(f(0)) += angle1*nor;
			normals.col(f(1)) += angle2*nor;
			normals.col(f(2)) += angle3*nor;
			angles[f(0)] += angle1;
			angles[f(1)] += angle2;
			angles[f(2)] += angle3;
		}
		for (int k = 0; k < nv; ++k) {
			Vector3d nor = normals.col(k) / angles[k];
			normals.col(k) = nor.normalized();
		}
		return normals;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Box geometry data
	///////////////////////////////////////////////////////////////////////////////

	double verts1_data[] = {
		-1, -1, -1,  1,
		-1, -1,  1,  1,
		-1,  1, -1,  1,
		-1,  1,  1,  1,
		1, -1, -1,  1,
		1, -1,  1,  1,
		1,  1, -1,  1,
		1,  1,  1,  1,
		-1,  0,  0,  1,
		1,  0,  0,  1,
		0, -1,  0,  1,
		0,  1,  0,  1,
		0,  0, -1,  1,
		0,  0,  1,  1,
	};

	int faces1_data[] = {
		0,     8,     2,
		1,     8,     0,
		3,     8,     1,
		2,     8,     3,
		4,    10,     0,
		5,    10,     4,
		1,    10,     5,
		0,    10,     1,
		6,     9,     4,
		7,     9,     6,
		5,     9,     7,
		4,     9,     5,
		2,    11,     6,
		3,    11,     2,
		7,    11,     3,
		6,    11,     7,
		1,    13,     3,
		5,    13,     1,
		7,    13,     5,
		3,    13,     7,
		0,    12,     4,
		2,    12,     0,
		6,    12,     2,
		4,    12,     6,
	};


	int edgeVerts1_data[] = {
		0, 1, 8,10,
		2, 0, 8,12,
		1, 3, 8,13,
		3, 2, 8,11,
		0, 4,10,12,
		5, 1,10,13,
		4, 5,10, 9,
		6, 2,11,12,
		4, 6, 9,12,
		3, 7,11,13,
		7, 5, 9,13,
		6, 7, 9,11,
	};

	// ADDED BY NICK
	// This is for keeping track of which edges each corner point is made up of
	int vertsEdge1_data[] = {
		0, 1, 4,
		0, 2, 5,
		1, 3, 7,
		2, 3, 9,
		4, 6, 8,
		5, 6, 10,
		7, 8, 11,
		9, 10, 11
	};

	double vertEdgeWeights1_data[] = {
		1.0,     0.0,     1.0,
		0.0,     1.0,     0.0,
		1.0,     0.0,     0.0,
		0.0,     1.0,     1.0,
		0.0,     1.0,     1.0,
		1.0,     0.0,     0.0,
		1.0,     0.0,     1.0,
		0.0,     1.0,     0.0
	};

	int edgeFaces1_data[] = {
		1,7,
		0,21,
		2,16,
		3,13,
		4,20,
		6,17,
		5,11,
		12,22,
		8,23,
		14,19,
		10,18,
		9,15
	};

	Map<Matrix<int, 3, 24, ColMajor> > faces1(faces1_data);
	Map<Matrix<int, 4, 12, ColMajor> > edgeVerts1(edgeVerts1_data);
	Map<Matrix<int, 3, 8, ColMajor> > vertEdges1(vertsEdge1_data); // ADDED BY NICK
	Map<Matrix<double, 3, 8, ColMajor> > vertEdgeWeights1(vertEdgeWeights1_data);
	Map<Matrix<int, 2, 12, ColMajor> > edgeFaces1(edgeFaces1_data);
	Map<Matrix<double, 4, 14, ColMajor> > verts1_(verts1_data);
	Matrix<double, 4, 14> verts1;
	Matrix<double, 3, 12> edgeNors1c;
	Matrix<double, 3, 12> edgeNors1d;
	MatrixXd faceNors1;
	MatrixXd vertNors1;

	void createBox(vector<shared_ptr<Edge> > &edges1, const Vector3d &whd1, const Matrix4d &E1)
	{
		Matrix4d S = Matrix4d::Identity();
		S(0, 0) = 0.5*whd1(0);
		S(1, 1) = 0.5*whd1(1);
		S(2, 2) = 0.5*whd1(2);
		Matrix4d E = E1 * S;
		verts1 = E * verts1_;
		faceNors1 = createFaceNormals(faces1, verts1);
		vertNors1 = createVertNormals(faces1, verts1);
		for (int k = 0; k < 12; ++k) {
			auto edge = make_shared<Edge>();
			edges1.push_back(edge);
			edge->verts = edgeVerts1.col(k);
			edge->faces = edgeFaces1.col(k);
			edge->internal = false; // for boxes only
			edge->normals[0] = faceNors1.col(edge->faces(0));
			edge->normals[1] = faceNors1.col(edge->faces(1));
			edge->angle = acos(edge->normals[0].dot(edge->normals[1]));
		}
	}

	///////////////////////////////////////////////////////////////////////////////

	// AABB Body
	void build_AABB_B(Matrix<double, 6, 1> &aabb, const MatrixXd &V)
	{
		aabb(0) = V.row(0).minCoeff();
		aabb(1) = V.row(1).minCoeff();
		aabb(2) = V.row(2).minCoeff();
		aabb(3) = V.row(0).maxCoeff();
		aabb(4) = V.row(1).maxCoeff();
		aabb(5) = V.row(2).maxCoeff();
	}

	// AABB face
	void build_AABB_F(MatrixXd &aabb, const MatrixXd &V, const MatrixXi &F)
	{
		for (int k = 0; k < F.cols(); ++k) {
			Vector3i f = F.col(k);
			Vector3d xa = V.block<3, 1>(0, f(0));
			Vector3d xb = V.block<3, 1>(0, f(1));
			Vector3d xc = V.block<3, 1>(0, f(2));
			aabb.block<3, 1>(0, k) = xa;
			aabb.block<3, 1>(3, k) = xa;
			for (int i = 0; i < 3; ++i) {
				aabb(i, k) = min(xb(i), aabb(i, k));
				aabb(i, k) = min(xc(i), aabb(i, k));
				aabb(i + 3, k) = max(xb(i), aabb(i + 3, k));
				aabb(i + 3, k) = max(xc(i), aabb(i + 3, k));
			}
		}
	}

	// AABB edge
	void build_AABB_E(MatrixXd &aabb, const MatrixXd &V, const vector<shared_ptr<Edge> > &edges)
	{
		for (int k = 0; k < edges.size(); ++k) {
			Vector4i e = edges[k]->verts;
			Vector3d x0 = V.block<3, 1>(0, e(0));
			Vector3d x1 = V.block<3, 1>(0, e(1));
			aabb.block<3, 1>(0, k) = x0;
			aabb.block<3, 1>(3, k) = x0;
			for (int i = 0; i < 3; ++i) {
				aabb(i, k) = min(x1(i), aabb(i, k));
				aabb(i + 3, k) = max(x1(i), aabb(i + 3, k));
			}
		}
	}

	// Check between two AABBs
	bool check_AABB(const Matrix<double, 6, 1> &aabb1, const Matrix<double, 6, 1> &aabb2)
	{
		double thresh = 1e-3;
		Vector3d threshVec(thresh, thresh, thresh);
		Vector3d min1 = aabb1.segment<3>(0) - threshVec;
		Vector3d max1 = aabb1.segment<3>(3) + threshVec;
		Vector3d min2 = aabb2.segment<3>(0) - threshVec;
		Vector3d max2 = aabb2.segment<3>(3) + threshVec;
		return
			max1(0) >= min2(0) &&
			min1(0) <= max2(0) &&
			max1(1) >= min2(1) &&
			min1(1) <= max2(1) &&
			max1(2) >= min2(2) &&
			min1(2) <= max2(2);
	}

	void barycentric(
		double &alpha, double &beta,
		const Vector3d &a, const Vector3d &b, const Vector3d &c,
		const Vector3d &p)
	{
		Vector3d v0 = b - a;
		Vector3d v1 = c - a;
		Vector3d v2 = p - a;
		double d00 = v0.dot(v0);
		double d01 = v0.dot(v1);
		double d11 = v1.dot(v1);
		double d20 = v2.dot(v0);
		double d21 = v2.dot(v1);
		double denom = d00 * d11 - d01 * d01;
		beta = (d11 * d20 - d01 * d21) / denom;
		double gamma = (d00 * d21 - d01 * d20) / denom;
		alpha = 1.0 - beta - gamma;
	}

	void lineline(
		double &a, double &b,
		const Vector3d &A1, const Vector3d &A2,
		const Vector3d &B1, const Vector3d &B2)
	{
		// https://www.mathworks.com/matlabcentral/newsreader/view_thread/246420
		// Closest points are
		// A0 = (1-a)*A1 + a*A2;
		// B0 = (1-b)*B1 + b*B2;
		Vector3d B2B1 = B2 - B1;
		Vector3d A1B1 = A1 - B1;
		Vector3d A2A1 = A2 - A1;
		Vector3d A2A1xB2B1 = A2A1.cross(B2B1);
		double nA = B2B1.cross(A1B1).dot(A2A1xB2B1);
		double nB = A2A1.cross(A1B1).dot(A2A1xB2B1);
		double d = A2A1.cross(B2B1).dot(A2A1xB2B1);
		a = nA / d;
		b = nB / d;
	}

	double linepoint(const Vector3d &A, const Vector3d &B, const Vector3d &P)
	{
		// Line: A -> B
		// Point: P
		Vector3d AP = P - A;
		Vector3d AB = B - A;
		double ab2 = AB.dot(AB);
		double apab = AP.dot(AB);
		return apab / ab2;
	}

	Vector3d linepointMinDist(const Vector3d &A, const Vector3d &B, const Vector3d &P)
	{
		// Line: A -> B
		// Point: P
		Vector3d AP = P - A;
		Vector3d AB = B - A;
		double ab2 = AB.dot(AB);
		double apab = AP.dot(AB);
		double t = apab / ab2;
		return (1.0 - t)*A + t*B;
	}

	int intersect_square(
		const Vector3d &x0, const Vector3d &dx,
		const Vector3d &xa, const Vector3d &xb, const Vector3d &xc,
		double &t)
	{
		// Intersects a ray with a square defined by three points
		//   x0: Ray origin
		//   dx: Ray direction
		//   xa, xb, xc: Triangle forming a quarter of the square
		//
		//  xb--------------xa
		//  | \            / |
		//  |   \        /   |
		//  |     \    /     |
		//  |       xc       |
		//  |     /    \     |
		//  |   /        \   |
		//  | /            \ |
		//  xd--------------xe
		//
		// Assume that xc is in the center of the quad

		int intersect = 0;
		double unused1, unused2;

		// Compute xd and xe
		Vector3d xd = xa + 2.0*(xc - xa);
		Vector3d xe = xb + 2.0*(xc - xb);

		// Slow: check all four triangles
		intersect = intersect_triangle3_inc(x0.data(), dx.data(), xa.data(), xb.data(), xc.data(), &t, &unused1, &unused2);
		if (intersect) {
			return 1;
		}
		intersect = intersect_triangle3_inc(x0.data(), dx.data(), xb.data(), xd.data(), xc.data(), &t, &unused1, &unused2);
		if (intersect) {
			return 1;
		}
		intersect = intersect_triangle3_inc(x0.data(), dx.data(), xd.data(), xe.data(), xc.data(), &t, &unused1, &unused2);
		if (intersect) {
			return 1;
		}
		intersect = intersect_triangle3_inc(x0.data(), dx.data(), xe.data(), xa.data(), xc.data(), &t, &unused1, &unused2);
		if (intersect) {
			return 1;
		}

		// No intersection
		t = -1.0;
		return 0;
	}

	void boxTriCollision(
		vector<shared_ptr<Collision> > &collisions,
		double threshold,
		const Vector3d &whd1,
		const Matrix4d &E1,
		const MatrixXd &verts2,
		const MatrixXi &faces2,
		const VectorXi &isEOL2,
		bool EOL)
	{
		vector<shared_ptr<Edge> > edges2;
		createEdges(edges2, faces2, verts2);
		boxTriCollision(collisions, threshold, whd1, E1, verts2, faces2, isEOL2, EOL, edges2);
	}

	void boxTriCollision(
		vector<shared_ptr<Collision> > &collisions,
		double threshold,
		const Vector3d &whd1,
		const Matrix4d &E1,
		const MatrixXd &verts2_,
		const MatrixXi &faces2,
		const VectorXi &isEOL2_,
		bool EOL,
		const vector<shared_ptr<Edge> > &edges2)
	{
		// Make a copy first so we can perturb
		MatrixXd verts2 = verts2_;

		// Assume isEOL=false by default
		VectorXi isEOL2(verts2.cols());
		if (isEOL2_.size() == verts2.cols()) {
			for (int i2 = 0; i2 < verts2.cols(); ++i2) {
				isEOL2(i2) = isEOL2_(i2);
			}
		}
		else {
			for (int i2 = 0; i2 < verts2.cols(); ++i2) {
				isEOL2(i2) = 0;
			}
		}

		// The first body is always the box
		vector<shared_ptr<Edge> > edges1;
		createBox(edges1, whd1, E1);

		// Perturb cloth verts
		std::random_device rd;
		std::mt19937 gen;
		std::uniform_real_distribution<> dis(-1.0, 1.0);
		gen.seed(1);
		for (int i2 = 0; i2 < verts2.cols(); ++i2) {
			Vector3d r;
			r(0) = dis(gen)*threshold*1e-3;
			r(1) = dis(gen)*threshold*1e-3;
			r(2) = dis(gen)*threshold*1e-3;
			verts2.block<3, 1>(0, i2) += r;
		}

		// Precompute the face normals for the cloth
		MatrixXd faceNors2 = createFaceNormals(faces2, verts2);

		// Build AABBs
		Matrix<double, 6, 1> aabbB1;
		Matrix<double, 6, 1> aabbB2;
		build_AABB_B(aabbB1, verts1);
		build_AABB_B(aabbB2, verts2);
		MatrixXd aabbF1(6, faces1.cols());
		build_AABB_F(aabbF1, verts1, faces1);
		MatrixXd aabbE2(6, edges2.size());
		build_AABB_E(aabbE2, verts2, edges2);

		// We don't need any vert2-face1 collisions when creating conformal geometry in EOL
		if (!EOL) {
			// Vertex2-Triangle1
			for (int i2 = 0; i2 < verts2.cols(); ++i2) {
				Vector3d x2 = verts2.block<3, 1>(0, i2);
				// AABB test: check Vertex2 against Body1
				Matrix<double, 6, 1> aabbV2;
				aabbV2.segment<3>(0) = x2;
				aabbV2.segment<3>(3) = x2;
				if (!check_AABB(aabbV2, aabbB1)) {
					continue;
				}
				shared_ptr<Collision> cmin = NULL;
				// Is this vertex inside Body2? If Body2 is convex, we can verify this
				// by doing a half-space test on all the triangles from Body2.
				bool inside = true;
				for (int j1 = 0; j1 < faces1.cols(); ++j1) {
					Vector3i f1 = faces1.col(j1);
					Vector3d x1a = verts1.block<3, 1>(0, f1(0));
					Vector3d nor = faceNors1.col(j1);
					Vector3d dx = x2 - x1a;
					if (nor.dot(dx) > 0.0) {
						inside = false;
						break;
					}
				}
				if (!inside) {
					continue;
				}
				for (int j1 = 0; j1 < faces1.cols(); ++j1) {
					Vector3i f1 = faces1.col(j1);
					Vector3d x1a = verts1.block<3, 1>(0, f1(0));
					Vector3d x1b = verts1.block<3, 1>(0, f1(1));
					Vector3d x1c = verts1.block<3, 1>(0, f1(2));
					// Project x2 onto the triangle 1
					Vector3d nor1 = faceNors1.col(j1);
					Vector3d dx = x2 - x1a;
					double proj = nor1.dot(dx);
					if (proj > 0.0) {
						// Outside the box
						continue;
					}
					Vector3d x1 = x2 - proj*nor1;
					dx = x2 - x1;
					double dist = dx.norm();
					if (dist > 5.0*threshold) {
						// Too far
						continue;
					}
					double u, v;
					barycentric(u, v, x1a, x1b, x1c, x1);
					double w = 1.0 - u - v;
					if (u < 0.0 || 1.0 < u || v < 0.0 || 1.0 < v || w < 0.0 || 1.0 < w) {
						// Projected point is outside the triangle
						continue;
					}
					// Compute cloth normal
					Vector3d nor2 = faceNors2.col(i2);
					if (nor2.dot(nor1) < 0.0) {
						nor2 = -nor2;
					}
					// Create contact object
					auto c = make_shared<Collision>();
					c->dist = dist;
					c->nor1 = nor1;
					c->nor2 = nor2;
					c->pos1 = x1;
					c->pos2 = x2;
					c->count1 = 3;
					c->count2 = 1;
					c->verts1 = f1;
					c->verts2 << i2, -1, -1;
					c->weights1 << u, v, w;
					c->weights2 << 1.0, 0.0, 0.0;
					c->tri1 = j1;
					c->tri2 = -1;
					// Is this the closest one so far?
					if (cmin == NULL) {
						cmin = c;
					}
					else {
						if (c->dist < cmin->dist) {
							cmin = c;
						}
					}
				}
				if (cmin != NULL) {
					collisions.push_back(cmin);
				}
			}
		}

		int nVertCol = collisions.size();
		int nFaceCol = 0;
		int nEdgeCol = 0;

		// Vertex1-Triangle2
		for (int i1 = 0; i1 < 8; ++i1) { // only up to 8 corner points (ignore face points)
			const Vector3d &x1 = verts1.block<3, 1>(0, i1);
			// AABB test: check Vertex1 against Body2
			Matrix<double, 6, 1> aabbV1;
			aabbV1.segment<3>(0) = x1;
			aabbV1.segment<3>(3) = x1;
			if (!check_AABB(aabbV1, aabbB2)) {
				continue;
			}
			shared_ptr<Collision> cmin = NULL;
			Vector3d nor1 = vertNors1.block<3, 1>(0, i1); // vertex normal
			for (int j2 = 0; j2 < faces2.cols(); ++j2) {
				const Vector3i &f2 = faces2.col(j2);
				const Vector3d &x2a = verts2.block<3, 1>(0, f2(0));
				const Vector3d &x2b = verts2.block<3, 1>(0, f2(1));
				const Vector3d &x2c = verts2.block<3, 1>(0, f2(2));
				Vector3d nor2 = faceNors2.col(j2);
				// Make sure the triangle normal points outward wrt the box.
				if (nor1.dot(nor2) < 0.0) {
					nor2 = -nor2;
				}
				// Is x1 on the correct side?
				Vector3d dx = x1 - x2a;
				double proj = dx.dot(nor2);
				if (proj < 0.0) {
					continue;
				}
				// Project x1 onto the triangle
				Vector3d x2 = x1 - proj*nor2;
				dx = x2 - x1;
				double dist = dx.norm();
				if (dist > 5.0*threshold) {
					// Too far
					continue;
				}
				// Compute barycentric coords of x1 wrt tri2
				double u, v;
				barycentric(u, v, x2a, x2b, x2c, x1);
				double w = 1.0 - u - v;
				if (u < 0.0 || 1.0 < u || v < 0.0 || 1.0 < v || w < 0.0 || 1.0 < w) {
					// Projected point is outside the triangle
					continue;
				}
				// Create contact object
				auto c = make_shared<Collision>();
				c->dist = dist;
				c->nor1 = nor1;
				c->nor2 = nor2;
				c->pos1 = x1;
				c->pos2 = x2;
				c->count1 = 1;
				c->count2 = 3;
				c->verts1 << i1, -1, -1;
				c->verts2 = f2;
				c->weights1 << 1.0, 0.0, 0.0;
				c->weights2 << u, v, w;
				c->edge1.push_back(vertEdges1(0, i1));
				c->edge1.push_back(vertEdges1(1, i1));
				c->edge1.push_back(vertEdges1(2, i1));
				c->tri1 = -1;
				c->tri2 = j2;
				// Is this the closest one so far?
				if (cmin == NULL) {
					cmin = c;
				}
				else {
					if (c->dist < cmin->dist) {
						cmin = c;
					}
				}
			}
			if (cmin != NULL) {
				collisions.push_back(cmin);
			}
		}
		nFaceCol = collisions.size() - nVertCol;

		// Edge2-Edge1
		for (int k2 = 0; k2 < edges2.size(); ++k2) {
			auto e2 = edges2[k2];
			const Vector3d &x2a = verts2.block<3, 1>(0, e2->verts(0));
			const Vector3d &x2b = verts2.block<3, 1>(0, e2->verts(1));
			Vector3d dx2 = x2b - x2a;
			double len2 = dx2.norm();
			Vector3d nor2 = (e2->normals[0] + e2->normals[1]).normalized();
			Matrix<double, 6, 1> aabbE2k = aabbE2.col(k2);
			for (int k1 = 0; k1 < edges1.size(); ++k1) {
				auto e1 = edges1[k1];
				if (e1->angle < M_PI / 6.0) {
					// Soft edge
					continue;
				}
				// AABB test: check the two triangles of Edge1 against Edge2
				bool aabbC = check_AABB(aabbF1.col(e1->faces(0)), aabbE2k);
				if (!aabbC) {
					bool aabbD = check_AABB(aabbF1.col(e1->faces(1)), aabbE2k);
					if (!aabbD) {
						continue;
					}
				}
				// x1a and x1b are the vertices of the edge.
				const Vector3d &x1a = verts1.block<3, 1>(0, e1->verts(0));
				const Vector3d &x1b = verts1.block<3, 1>(0, e1->verts(1));
				Vector3d dx1 = x1b - x1a;
				double len1 = dx1.norm();
				Vector3d tan1 = dx1 / len1;
				// Is the box edge parallel to the cloth normal?
				double threshAng = 2.0*M_PI / 180.0;
				double angle = acos(tan1.dot(nor2));
				if (fabs(angle) < threshAng || fabs(M_PI - angle) < threshAng) {
					continue;
				}
				// Are the two edges parallel?
				angle = acos(tan1.dot(dx2) / len2);
				if (fabs(angle) < threshAng || fabs(M_PI - angle) < threshAng) {
					continue;
				}
				Vector3d nor = dx1.cross(dx2).normalized();
				// The two triangles of the edge are (a,b,c) and (b,a,d), and n1c 
				// and n1d are the two triangle normals
				const Vector3d &x1c = verts1.block<3, 1>(0, e1->verts(2));
				const Vector3d &x1d = verts1.block<3, 1>(0, e1->verts(3));
				const Vector3d &n1c = e1->normals[0];
				const Vector3d &n1d = e1->normals[1];
				// Make the computed normal point along the edge normal.
				Vector3d nor1 = n1c + n1d;
				nor1.normalize();
				if (nor.dot(nor1) < 0.0) {
					nor = -nor;
				}
				// The computed normal must lie between the two edge normals.
				// n1d nor
				//  | /
				//  |/   
				//  +--- n1c
				double angleCD = acos(n1c.dot(n1d)); // should be the same as e1.angle modulo sign
				double angleCN = acos(n1c.dot(nor)); // angle between n1c to the computed normal
				if (angleCD < 0.0) {
					angleCD = -angleCD;
					angleCN = -angleCN;
				}
				// angleCN must be between 0 and angleCD
				if (angleCN < -threshAng || angleCN - angleCD > threshAng) {
					continue;
				}
				// Where does the line intersect the box?
				double u2c, u2d, unused1, unused2;
				//int i2c = intersect_triangle3_inc(x2a.data(), dx2.data(), x1a.data(), x1b.data(), x1c.data(), &u2c, &unused1, &unused2);
				//int i2d = intersect_triangle3_inc(x2a.data(), dx2.data(), x1b.data(), x1a.data(), x1d.data(), &u2d, &unused1, &unused2);
				int i2c = intersect_square(x2a, dx2, x1a, x1b, x1c, u2c);
				int i2d = intersect_square(x2a, dx2, x1b, x1a, x1d, u2d);
				// Are these ray collisions line segment collisions?
				i2c = i2c && (0.0 <= u2c && u2c <= 1.0);
				i2d = i2d && (0.0 <= u2d && u2d <= 1.0);
				// Compute closest points on the two lines.
				double u1, u2;
				lineline(u1, u2, x1a, x1b, x2a, x2b);
				// Check if this is a vertex (rather than edge) collision. It is an
				// edge collision if both u1 and u2 are between 0 and 1.
				// For the threshold, use edgeLen to convert to world units.
				double thresh1 = 1.0*threshold / len1;
				double thresh2 = 1.0*threshold / len2;
				if (u1 < -thresh1 || u1 > 1.0 + thresh1 || u2 < -thresh2 || u2 > 1.0 + thresh2) {
					// This is a vertex collision.
					continue;
				}
				Vector3d x1, x2, dx;
				if (i2c && i2d) {
					// The cloth edge intersects both box triangles
					x1 = (1.0 - u1)*x1a + u1*x1b;
					x2 = (1.0 - u2)*x2a + u2*x2b;
				}
				else if (i2c && !i2d) {
					// Only intersects Triangle C.
					// Is x2a inside the box?
					dx = x2a - x1a;
					if (dx.dot(n1c) < 0.0) {
						// Inside... the valid range for u2 is [0, u2c].
						u2 = max(0.0, min(u2c, u2));
					}
					else {
						// Outside... the valid range for u2 is [u2c, 1.0].
						u2 = max(u2c, min(1.0, u2));
					}
					x2 = (1.0 - u2)*x2a + u2*x2b;
					// Recompute u1 based on the new x2.
					u1 = linepoint(x1a, x1b, x2);
					x1 = (1.0 - u1)*x1a + u1*x1b;
				}
				else if (!i2c && i2d) {
					// Only intersects Triangle D.
					// Is x2a inside the box?
					dx = x2a - x1b;
					if (dx.dot(n1d) < 0.0) {
						// Inside... the valid range for u2 is [0, u2d].
						u2 = max(0.0, min(u2d, u2));
					}
					else {
						// Outside... the valid range for u2 is [u2d, 1.0].
						u2 = max(u2d, min(1.0, u2));
					}
					x2 = (1.0 - u2)*x2a + u2*x2b;
					// Recompute u1 based on the new x2.
					u1 = linepoint(x1a, x1b, x2);
					x1 = (1.0 - u1)*x1a + u1*x1b;
				}
				else {
					// The cloth edge does not intersect either of the box triangles
					continue;
				}
				if (u1 < -thresh1 || u1 > 1.0 + thresh1 || u2 < -thresh2 || u2 > 1.0 + thresh2) {
					// This is a vertex, not edge, collision. We need another check
					// here because we modify u1 and u2 in the if block above.
					continue;
				}
				// At this point, x2 is the collision point on the cloth and x1 is
				// the collision point on the box.
				dx = x2 - x1;
				double thresh = 2.0*threshold;
				if (dx.dot(dx) > thresh*thresh) {
					// Too far
					continue;
				}
				// Final check: are the collided points close to the box edge?
				thresh = 2.0*threshold;
				dx = x2 - x1;
				if (dx.dot(dx) > thresh*thresh) {
					continue;
				}
				auto c = make_shared<Collision>();
				c->dist = dx.norm();
				c->nor1 = nor1;
				c->nor2 = nor;
				c->pos1 = x1;
				c->pos2 = x2;
				c->count1 = 2;
				c->count2 = 2;
				c->verts1 << e1->verts.segment<2>(0), -1;
				c->verts2 << e2->verts.segment<2>(0), -1;
				c->weights1 << 1.0 - u1, u1, 0.0;
				c->weights2 << 1.0 - u2, u2, 0.0;
				c->edge1.push_back(k1);
				c->edge2 = k2;
				c->edgeDir = tan1;
				collisions.push_back(c);
			}
		}
		nEdgeCol = collisions.size() - nFaceCol - nVertCol;

		// Each box corner should have at most 1 collision with the cloth (can't be
		// two or more). We'll keep the closest cloth collision to the box corner.
		for (int i1 = 0; i1 < 8; ++i1) {
			const Vector3d &x1 = verts1.block<3, 1>(0, i1);
			// Find the closest collision to the box corner
			int kmin = -1;
			double dmin = 1e9;
			for (int k = 0; k < collisions.size(); ++k) {
				if (collisions[k]->count1 == 1 && collisions[k]->verts1(0) == i1) {
					const Vector3d &x2 = collisions[k]->pos2;
					Vector3d dx = x2 - x1;
					double d = dx.dot(dx);
					if (d < dmin) {
						kmin = k;
						dmin = d;
					}
				}
			}
			if (kmin != -1) {
				// Make a list of collisions to delete
				vector<int> dlist;
				for (int k = 0; k < collisions.size(); ++k) {
					if (collisions[k]->count1 == 1 && collisions[k]->verts1(0) == i1 && k != kmin) {
						dlist.push_back(k);
					}
				}
				// Delete collisions in reverse order
				for (int kdel : dlist) {
					collisions[kdel] = collisions.back();
					collisions.pop_back();
				}
			}
		}

		// Go through the collisions and create `pos1_`, the position of the
		// collision on the box offset by a small amount so that it is inside the
		// box.
		double snapDepth = 0.1*threshold;
		for (auto collision : collisions) {
			collision->pos1_ = collision->pos1 - snapDepth*collision->nor1;
		}

		return;
	}

	///////////////////////////////////////////////////////////////////////////////

	void pointTriCollision(
		vector<shared_ptr<Collision> > &collisions,
		double threshold,
		const Eigen::MatrixXd &verts1,
		const Eigen::MatrixXd &norms1,
		const Eigen::MatrixXd &verts2_,
		const Eigen::MatrixXi &faces2,
		bool EOL)
	{
		// Make a copy first so we can perturb
		MatrixXd verts2 = verts2_;

		// Perturb cloth verts
		std::random_device rd;
		std::mt19937 gen;
		std::uniform_real_distribution<> dis(-1.0, 1.0);
		gen.seed(1);
		for (int i2 = 0; i2 < verts2.cols(); ++i2) {
			Vector3d r;
			r(0) = dis(gen)*threshold*1e-3;
			r(1) = dis(gen)*threshold*1e-3;
			r(2) = dis(gen)*threshold*1e-3;
			verts2.block<3, 1>(0, i2) += r;
		}

		// Precompute the face normals for the cloth
		MatrixXd faceNors2 = createFaceNormals(faces2, verts2);

		// Build AABBs
		Matrix<double, 6, 1> aabbB2;
		build_AABB_B(aabbB2, verts2);

		if (EOL) {
			for (int i2 = 0; i2 < verts2.cols(); ++i2) {
				Vector3d x2 = verts2.block<3, 1>(0, i2);
				shared_ptr<Collision> cmin = NULL;
				for (int i1 = 0; i1 < verts1.cols(); ++i1) {
					Vector3d x1 = verts1.block<3, 1>(0, i1);
					Vector3d dx = x2 - x1;
					double dist = dx.norm();
					if (dist < threshold) {
						Vector3d nor1 = norms1.block<3, 1>(0, i1); // vertex normal
						// Create contact object
						auto c = make_shared<Collision>();
						c->dist = dist;
						c->nor1 = nor1;
						c->nor2 = nor1; // We don't care so hack
						c->pos1 = x1;
						c->pos2 = x2;
						c->count1 = 3;
						c->count2 = 1;
						c->verts1 << i1, -1, -1;
						c->verts2 << i2, -1, -1;
						c->weights1 << 1.0, 0.0, 0.0;
						c->weights2 << 1.0, 0.0, 0.0;
						c->tri1 = -1;
						c->tri2 = -1;
						// Is this the closest one so far?
						if (cmin == NULL) {
							cmin = c;
						}
						else {
							if (c->dist < cmin->dist) {
								cmin = c;
							}
						}
					}
				}
				if (cmin != NULL) {
					collisions.push_back(cmin);
				}
			}
		}

		// Vertex1-Triangle2
		for (int i1 = 0; i1 < verts1.cols(); ++i1) {
			const Vector3d &x1 = verts1.block<3, 1>(0, i1);
			// AABB test: check Vertex1 against Body2
			Matrix<double, 6, 1> aabbV1;
			aabbV1.segment<3>(0) = x1;
			aabbV1.segment<3>(3) = x1;
			if (!check_AABB(aabbV1, aabbB2)) {
				continue;
			}
			shared_ptr<Collision> cmin = NULL;
			Vector3d nor1 = norms1.block<3, 1>(0, i1); // vertex normal
			for (int j2 = 0; j2 < faces2.cols(); ++j2) {
				const Vector3i &f2 = faces2.col(j2);
				const Vector3d &x2a = verts2.block<3, 1>(0, f2(0));
				const Vector3d &x2b = verts2.block<3, 1>(0, f2(1));
				const Vector3d &x2c = verts2.block<3, 1>(0, f2(2));
				Vector3d nor2 = faceNors2.col(j2);
				// Make sure the triangle normal points outward wrt the box.
				if (nor1.dot(nor2) < 0.0) {
					nor2 = -nor2;
				}
				// Is x1 on the correct side?
				Vector3d dx = x1 - x2a;
				double proj = dx.dot(nor2);
				if (proj < 0.0) {
					continue;
				}
				// Project x1 onto the triangle
				Vector3d x2 = x1 - proj*nor2;
				dx = x2 - x1;
				double dist = dx.norm();
				if (dist > 5.0*threshold) {
					// Too far
					continue;
				}
				// Compute barycentric coords of x1 wrt tri2
				double u, v;
				barycentric(u, v, x2a, x2b, x2c, x1);
				double w = 1.0 - u - v;
				if (u < 0.0 || 1.0 < u || v < 0.0 || 1.0 < v || w < 0.0 || 1.0 < w) {
					// Projected point is outside the triangle
					continue;
				}
				// Create contact object
				auto c = make_shared<Collision>();
				c->dist = dist;
				c->nor1 = nor1;
				c->nor2 = nor2;
				c->pos1 = x1;
				c->pos2 = x2;
				c->count1 = 1;
				c->count2 = 3;
				c->verts1 << i1, -1, -1;
				c->verts2 = f2;
				c->weights1 << 1.0, 0.0, 0.0;
				c->weights2 << u, v, w;
				c->tri1 = -1;
				c->tri2 = j2;
				// Is this the closest one so far?
				if (cmin == NULL) {
					cmin = c;
				}
				else {
					if (c->dist < cmin->dist) {
						cmin = c;
					}
				}
			}
			if (cmin != NULL) {
				collisions.push_back(cmin);
			}
		}

		// Go through the collisions and create `pos1_`, the position of the
		// collision on the box offset by a small amount so that it is inside the
		// box.
		double snapDepth = 0.1*threshold;
		for (auto collision : collisions) {
			collision->pos1_ = collision->pos1 - snapDepth*collision->nor1;
		}

		return;
	}

}	