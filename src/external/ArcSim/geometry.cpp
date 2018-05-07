/*
Copyright ©2013 The Regents of the University of California
(Regents). All Rights Reserved. Permission to use, copy, modify, and
distribute this software and its documentation for educational,
research, and not-for-profit purposes, without fee and without a
signed licensing agreement, is hereby granted, provided that the
above copyright notice, this paragraph and the following two
paragraphs appear in all copies, modifications, and
distributions. Contact The Office of Technology Licensing, UC
Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
(510) 643-7201, for commercial licensing opportunities.

IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include "geometry.hpp"
#include <cstdlib>

using namespace std;

double signed_vf_distance(const Vec3 &x,
	const Vec3 &y0, const Vec3 &y1, const Vec3 &y2,
	Vec3 *n, double *w) {
	Vec3 _n; if (!n) n = &_n;
	double _w[4]; if (!w) w = _w;
	*n = cross(normalize(y1 - y0), normalize(y2 - y0));
	if (norm2(*n) < 1e-6)
		return infinity;
	*n = normalize(*n);
	double h = dot(x - y0, *n);
	double b0 = stp(y1 - x, y2 - x, *n),
		b1 = stp(y2 - x, y0 - x, *n),
		b2 = stp(y0 - x, y1 - x, *n);
	w[0] = 1;
	w[1] = -b0 / (b0 + b1 + b2);
	w[2] = -b1 / (b0 + b1 + b2);
	w[3] = -b2 / (b0 + b1 + b2);
	return h;
}

double signed_ve_distance(const Vec3 &x, const Vec3 &y0, const Vec3 &y1,
	Vec3* n, double *w) {
	Vec3 e = y1 - y0;
	double d = dot(x - y0, e) / norm2(e);
	if (d < 0 || d > 1.0)
		return infinity;
	if (w) {
		w[0] = 1;
		w[1] = -(1.0 - d);
		w[2] = -d;
		w[3] = 0;
	}
	Vec3 dist = x - (y0 + d*e);
	double l = norm(dist);
	if (n && fabs(l) > 1e-16) *n = dist / l;
	return l;
}

double signed_ee_distance(const Vec3 &x0, const Vec3 &x1,
	const Vec3 &y0, const Vec3 &y1,
	Vec3 *n, double *w) {
	Vec3 _n; if (!n) n = &_n;
	double _w[4]; if (!w) w = _w;
	*n = cross(normalize(x1 - x0), normalize(y1 - y0));
	if (norm2(*n) < 1e-8) {
		// special case: parallel lines
		Vec3 e0 = normalize(x1 - x0), e1 = normalize(y1 - y0);

		double p0min = dot(x0, e0), p0max = dot(x1, e0), p1min = dot(y0, e0), p1max = dot(y1, e0);
		if (p1max < p1min) swap(p1max, p1min);

		double a = max(p0min, p1min), b = min(p0max, p1max), c = 0.5*(a + b);
		if (a > b) return infinity;

		Vec3 d = (y0 - x0) - dot(y0 - x0, e0)*e0;

		if (n) *n = normalize(-d);
		if (w) {
			w[1] = (c - dot(x0, e0)) / norm(x1 - x0);
			w[0] = 1.0 - w[1];
			w[3] = -(dot(e0, e1)*c - dot(y0, e1)) / norm(y1 - y0);
			w[2] = -1.0 - w[3];
		}
		return norm(d);
	}
	*n = normalize(*n);
	double h = dot(x0 - y0, *n);
	double a0 = stp(y1 - x1, y0 - x1, *n), a1 = stp(y0 - x0, y1 - x0, *n),
		b0 = stp(x0 - y1, x1 - y1, *n), b1 = stp(x1 - y0, x0 - y0, *n);
	w[0] = a0 / (a0 + a1);
	w[1] = a1 / (a0 + a1);
	w[2] = -b0 / (b0 + b1);
	w[3] = -b1 / (b0 + b1);
	return h;
}

bool set_unsigned_ve_distance(const Vec3 &x, const Vec3 &y0, const Vec3 &y1,
	double *_d, Vec3 *_n,
	double *_wx, double *_wy0, double *_wy1) {
	double t = clamp(dot(x - y0, y1 - y0) / dot(y1 - y0, y1 - y0), 0., 1.);
	Vec3 y = y0 + t*(y1 - y0);
	double d = norm(x - y);
	if (d < *_d) {
		*_d = d;
		*_n = normalize(x - y);
		*_wx = 1;
		*_wy0 = 1 - t;
		*_wy1 = t;
		return true;
	}
	return false;
}

bool set_unsigned_vf_distance(const Vec3 &x,
	const Vec3 &y0, const Vec3 &y1, const Vec3 &y2,
	double *_d, Vec3 *_n,
	double *_wx,
	double *_wy0, double *_wy1, double *_wy2) {
	Vec3 n = normalize(cross(normalize(y1 - y0), normalize(y2 - y0)));
	double d = abs(dot(x - y0, n));
	double b0 = stp(y1 - x, y2 - x, n),
		b1 = stp(y2 - x, y0 - x, n),
		b2 = stp(y0 - x, y1 - x, n);
	if (d < *_d && b0 >= 0 && b1 >= 0 && b2 >= 0) {
		*_d = d;
		*_n = n;
		*_wx = 1;
		*_wy0 = -b0 / (b0 + b1 + b2);
		*_wy1 = -b1 / (b0 + b1 + b2);
		*_wy2 = -b2 / (b0 + b1 + b2);
		return true;
	}
	bool success = false;
	if (b0 < 0
		&& set_unsigned_ve_distance(x, y1, y2, _d, _n, _wx, _wy1, _wy2)) {
		success = true;
		*_wy0 = 0;
	}
	if (b1 < 0
		&& set_unsigned_ve_distance(x, y2, y0, _d, _n, _wx, _wy2, _wy0)) {
		success = true;
		*_wy1 = 0;
	}
	if (b2 < 0
		&& set_unsigned_ve_distance(x, y0, y1, _d, _n, _wx, _wy0, _wy1)) {
		success = true;
		*_wy2 = 0;
	}
	return success;
}

bool set_unsigned_ee_distance(const Vec3 &x0, const Vec3 &x1,
	const Vec3 &y0, const Vec3 &y1,
	double *_d, Vec3 *_n,
	double *_wx0, double *_wx1,
	double *_wy0, double *_wy1) {
	Vec3 n = normalize(cross(normalize(x1 - x0), normalize(y1 - y0)));
	double d = abs(dot(x0 - y0, n));
	double a0 = stp(y1 - x1, y0 - x1, n), a1 = stp(y0 - x0, y1 - x0, n),
		b0 = stp(x0 - y1, x1 - y1, n), b1 = stp(x1 - y0, x0 - y0, n);
	if (d < *_d && a0 >= 0 && a1 >= 0 && b0 >= 0 && b1 >= 0) {
		*_d = d;
		*_n = n;
		*_wx0 = a0 / (a0 + a1);
		*_wx1 = a1 / (a0 + a1);
		*_wy0 = -b0 / (b0 + b1);
		*_wy1 = -b1 / (b0 + b1);
		return true;
	}
	bool success = false;
	if (a0 < 0
		&& set_unsigned_ve_distance(x1, y0, y1, _d, _n, _wx1, _wy0, _wy1)) {
		success = true;
		*_wx0 = 0;
	}
	if (a1 < 0
		&& set_unsigned_ve_distance(x0, y0, y1, _d, _n, _wx0, _wy0, _wy1)) {
		success = true;
		*_wx1 = 0;
	}
	if (b0 < 0
		&& set_unsigned_ve_distance(y1, x0, x1, _d, _n, _wy1, _wx0, _wx1)) {
		success = true;
		*_wy0 = 0;
		*_n = -*_n;
	}
	if (b1 < 0
		&& set_unsigned_ve_distance(y0, x0, x1, _d, _n, _wy0, _wx0, _wx1)) {
		success = true;
		*_wy1 = 0;
		*_n = -*_n;
	}
	return success;
}

double unsigned_vf_distance(const Vec3 &x,
	const Vec3 &y0, const Vec3 &y1, const Vec3 &y2,
	Vec3 *n, double w[4]) {
	Vec3 _n; if (!n) n = &_n;
	double _w[4]; if (!w) w = _w;
	double d = infinity;
	set_unsigned_vf_distance(x, y0, y1, y2, &d, n, &w[0], &w[1], &w[2], &w[3]);
	return d;
}

double unsigned_ee_distance(const Vec3 &x0, const Vec3 &x1,
	const Vec3 &y0, const Vec3 &y1,
	Vec3 *n, double w[4]) {
	Vec3 _n; if (!n) n = &_n;
	double _w[4]; if (!w) w = _w;
	double d = infinity;
	set_unsigned_ee_distance(x0, x1, y0, y1, &d, n, &w[0], &w[1], &w[2], &w[3]);
	return d;
}

// NICK
double unsigned_vv_distance(const Vec3 &x, const Vec3 &y) {
	return sqrt(pow(y[0] - x[0],2) + pow(y[1] - x[1], 2));
}

double incedent_angle(const Vert* v, const Face* f) {
	Vec3 x = v->u;
	Vec3 y;
	Vec3 z;
	if (f->v[0] == v) {
		y = f->v[1]->u;
		z = f->v[2]->u;
	}
	else if (f->v[1] == v) {
		y = f->v[0]->u;
		z = f->v[2]->u;
	}
	else if (f->v[2] == v) {
		y = f->v[0]->u;
		z = f->v[1]->u;
	}
	Vec3 a = y - x;
	Vec3 b = z - x;
	return acos(dot(a, b) / (norm(a) * norm(b)));
}

Vec3 get_barycentric_coords(const Vec2& point, const Face* f) {
	// WHY IS THERE A WARNING HERE?
	//cout << "warning bary" << endl;
	// Compute vectors        
	Vec2 v0 = reduce_xy(f->v[0]->u - f->v[2]->u);
	Vec2 v1 = reduce_xy(f->v[1]->u - f->v[2]->u);
	Vec2 v2 = point - reduce_xy(f->v[2]->u);
	// Compute dot products
	double dot00 = dot(v0, v0);
	double dot01 = dot(v0, v1);
	double dot02 = dot(v0, v2);
	double dot11 = dot(v1, v1);
	double dot12 = dot(v1, v2);
	// Compute barycentric coordinates
	double invDenom = 1.f / (dot00 * dot11 - dot01 * dot01);
	double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
	return Vec3(u, v, 1 - u - v);
}

// Is the point within the face?
// Adapted from http://www.blackpawn.com/texts/pointinpoly/default.html
bool is_inside(const Vec2& point, const Face* f) {
	Vec3 bary = get_barycentric_coords(point, f);
	//printf("UV: %f, %f\n", u, v);
	// Check if point is in triangle
	// 10*epsilon: want to be robust for borders
	return ((bary[0] >= -10 * EPSILON) && (bary[1] >= -10 * EPSILON) && (bary[2] >= -100 * EPSILON));
}

// Gets the face that surrounds point u in material space
Face* get_enclosing_face(const Mesh& mesh, const Vec2& u,
	Face *starting_face_hint) {
	for (int f = 0; f < (int)mesh.faces.size(); f++)
		if (is_inside(u, mesh.faces[f]))
			return mesh.faces[f];
	return NULL;
}

template <Space s> Vec3 normal(const Face *face) {
	const Vec3 &x0 = pos<s>(face->v[0]),
		&x1 = pos<s>(face->v[1]),
		&x2 = pos<s>(face->v[2]);
	return normalize(cross(x1 - x0, x2 - x0));
}
template Vec3 normal<MS>(const Face *face);
template Vec3 normal<PS>(const Face *face);
template Vec3 normal<WS>(const Face *face);

template <Space s> Vec3 normal(const Node* node) {
	Vec3 n(0);
	for (size_t v = 0; v < node->verts.size(); v++) {
		const Vert *vert = node->verts[v];
		const vector<Face*> &adjfs = vert->adjf;
		for (size_t i = 0; i < adjfs.size(); i++) {
			Face const* face = adjfs[i];
			int j = find(vert, face->v), j1 = (j + 1) % 3, j2 = (j + 2) % 3;
			Vec3 e1 = pos<s>(face->v[j1]) - pos<s>(vert),
				e2 = pos<s>(face->v[j2]) - pos<s>(vert);
			n += cross(e1, e2) / (2 * norm2(e1)*norm2(e2));
		}
	}
	return normalize(n);
}
template Vec3 normal<MS>(const Node*);
template Vec3 normal<WS>(const Node*);

double dihedral_angle(const Vec3& p0, const Vec3& p1, const Vec3& n0, const Vec3& n1) {
	Vec3 e = normalize(p1 - p0);
	if (norm2(e) == 0) return 0;
	if (norm2(n0) == 0 || norm2(n1) == 0) return 0;
	double cosine = dot(n0, n1), sine = dot(e, cross(n0, n1));
	double theta = atan2(sine, cosine);
	return theta;
}

template <Space s> double dihedral_angle(const Edge *edge) {
	// if (!hinge.edge[0] || !hinge.edge[1]) return 0;
	// const Edge *edge0 = hinge.edge[0], *edge1 = hinge.edge[1];
	// int s0 = hinge.s[0], s1 = hinge.s[1];
	if (!edge->adjf[0] || !edge->adjf[1])
		return 0;
	if (s == MS && is_seam(edge))
		return 0;
	Vec3 e = normalize(pos<s>(edge_vert(edge, 0, 0)) - pos<s>(edge_vert(edge, 0, 1)));
	if (norm2(e) == 0) return 0;
	Vec3 n0 = normal<s>(edge->adjf[0]), n1 = normal<s>(edge->adjf[1]);
	if (norm2(n0) == 0 || norm2(n1) == 0) return 0;
	double cosine = dot(n0, n1), sine = dot(e, cross(n0, n1));
	double theta = atan2(sine, cosine);
	return theta;
}
template double dihedral_angle<MS>(const Edge *edge);
template double dihedral_angle<PS>(const Edge *edge);
template double dihedral_angle<WS>(const Edge *edge);

template <Space s> Mat2x2 projected_curvature(const Face *face, const Mat2x3& base) {
	Mat2x2 S;
	for (int e = 0; e < 3; e++) {
		Vec2 e_mat = base * (face->v[PREV(e)]->u - face->v[NEXT(e)]->u),
			t_mat = perp(normalize(e_mat));

		// TODO: don't refine over sharp creases
		//if ((edge->v[0]->flag & edge->v[1]->flag & Vert::SharpCrease) != 0)
		//    continue;
		double theta = dihedral_angle<s>(face->adje[e]);
		S -= 1 / 2.*theta*norm(e_mat)*outer(t_mat, t_mat);
	}
	S /= face->a;
	return S;
}
template Mat2x2 projected_curvature<PS>(const Face *face, const Mat2x3&);
template Mat2x2 projected_curvature<WS>(const Face *face, const Mat2x3&);

template <Space s> Mat3x3 curvature(const Face *face) {
	Mat3x3 S;
	Vec3 n = normal<MS>(face);
	for (int e = 0; e < 3; e++) {
		Vec3 e_mat = face->v[PREV(e)]->u - face->v[NEXT(e)]->u,
			t_mat = cross(normalize(e_mat), n);

		// TODO: don't refine over sharp creases
		//if ((edge->v[0]->flag & edge->v[1]->flag & Vert::SharpCrease) != 0)
		//    continue;
		double theta = dihedral_angle<s>(face->adje[e]);
		S -= 1 / 2.*theta*norm(e_mat)*outer(t_mat, t_mat);
	}
	S /= face->a;
	return S;
}
template Mat3x3 curvature<MS>(const Face *face);
template Mat3x3 curvature<PS>(const Face *face);
template Mat3x3 curvature<WS>(const Face *face);

double aspect(const Vec3& u0, const Vec3& u1, const Vec3& u2) {
	double perimeter = norm(u0 - u1) + norm(u1 - u2) + norm(u2 - u0);
	return 12 * sqrt(3)* 0.5 * norm(cross(u1 - u0, u2 - u0)) / sq(perimeter);
}

double get_angle(const Vec3& u, const Vec3& v) {
	return acos(clamp(dot(normalize(u), normalize(v)), -1.0, 1.0));
}

// Planes

template<Space s>
Plane plane_fit(const Mesh& mesh) {
	int num = 0;
	Vec3 center(0);
	for (int i = 0; i<mesh.verts.size(); i++) {
		num++;
		center += pos<s>(mesh.verts[i]);
	}
	center *= 1.0 / num;

	Mat3x3 M(0);
	for (int i = 0; i<mesh.verts.size(); i++) {
		Vec3 p = pos<s>(mesh.verts[i]) - center;
		M += outer(p, p);
	}
	Eig<3> eig = eigen_decomposition(M);
	return Plane(center, eig.Q.col(2));
}

Mat3x3 local_base(const Vec3& normal) {
	//Vec3 u = (dot(normal, Vec3(1, 0, 0)) > dot(normal, Vec3(0, 0, 1))) ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
	Vec3 u = (abs(normal[2]) < 1e-6) ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
	u = normalize(u - dot(u, normal)*normal);
	Vec3 v = cross(normal, u);
	return Mat3x3(u, v, normal);
}

template Plane plane_fit<MS>(const Mesh&);
template Plane plane_fit<WS>(const Mesh&);

// Intersection

// triangle x_I, ray x+t*d
bool triangle_ray_test(const Vec3 &x0, const Vec3& x1, const Vec3& x2,
	const Vec3 &s, const Vec3& d, double& z, Vec3* bary) {
	Vec3 e1 = x1 - x0, e2 = x2 - x0, t = s - x0;
	Vec3 p = cross(d, e2), q = cross(t, e1);
	double div = dot(p, e1);
	if (fabs(div) < 1e-18) return false;

	double idiv = 1.0 / div;
	double u = idiv * dot(p, t), v = idiv * dot(q, d);

	if (u<-1e-6 || v<-1e-6 || (u + v - 1.0) > 1e-6) return false;
	z = idiv * dot(q, e2);
	if (bary)
		*bary = Vec3((1.0 - u - v), u, v);
	return true;
}
