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

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "mesh.hpp"
#include "util.hpp"

double signed_vf_distance(const Vec3 &x,
	const Vec3 &y0, const Vec3 &y1, const Vec3 &y2,
	Vec3 *n, double *w);

double signed_ee_distance(const Vec3 &x0, const Vec3 &x1,
	const Vec3 &y0, const Vec3 &y1,
	Vec3 *n, double *w);

double unsigned_vf_distance(const Vec3 &x,
	const Vec3 &y0, const Vec3 &y1, const Vec3 &y2,
	Vec3 *n, double *w);

double unsigned_ee_distance(const Vec3 &x0, const Vec3 &x1,
	const Vec3 &y0, const Vec3 &y1,
	Vec3 *n, double *w);

double signed_ve_distance(const Vec3 &x, const Vec3 &y0, const Vec3 &y1,
	Vec3* n, double *w);

// NICK
double unsigned_vv_distance(const Vec3 &x, const Vec3 &y);
double incedent_angle(const Vert* v, const Face* f);

Vec3 get_barycentric_coords(const Vec2 &point, const Face *face);

Face* get_enclosing_face(const Mesh& mesh, const Vec2& u,
	Face *starting_face_hint = NULL);

enum Space { MS, PS, WS }; // material space, plastic space, world space

template <Space s> const Vec3 &pos(const Node *node);
template <Space s> inline Vec3 &pos(Node *node);
template <Space s> const Vec3 &pos(const Vert *vert);
template <Space s> inline Vec3 &pos(Vert *vert);
template <Space s> Vec3 normal(const Face *face);
template <Space s> Vec3 normal(const Node* node);
template <Space s> double dihedral_angle(const Edge *edge);
double dihedral_angle(const Vec3& p0, const Vec3& p1, const Vec3& n0, const Vec3& n1);
template <Space s> Mat3x3 curvature(const Face *face);
template <Space s> Mat2x2 projected_curvature(const Face *face, const Mat2x3& base);

inline double area(const Vec3& u0, const Vec3& u1, const Vec3& u2) { return 0.5*norm(cross(u1 - u0, u2 - u0)); }
inline double area(const Face* face) { return area(face->v[0]->u, face->v[1]->u, face->v[2]->u); }
double aspect(const Vec3& u0, const Vec3& u1, const Vec3& u2);
inline double aspect(const Face* face) { return aspect(face->v[0]->u, face->v[1]->u, face->v[2]->u); }
double get_angle(const Vec3& u, const Vec3& v);

// Nick
inline double edge_length(const Edge* edge) { return sqrt(pow(edge->n[0]->verts[0]->u[0] - edge->n[1]->verts[0]->u[0], 2) + pow(edge->n[0]->verts[0]->u[1] - edge->n[1]->verts[0]->u[1], 2)); }
inline double edge_length3D(const Edge* edge) { return sqrt(pow(edge->n[0]->x[0] - edge->n[1]->x[0], 2) + pow(edge->n[0]->x[1] - edge->n[1]->x[1], 2) + pow(edge->n[0]->x[2] - edge->n[1]->x[2], 2)); }

template <Space s> Plane plane_fit(const Mesh& mesh);

Mat3x3 local_base(const Vec3& normal);

bool triangle_ray_test(const Vec3 &x0, const Vec3& x1, const Vec3& x2,
	const Vec3 &s, const Vec3& d, double& z, Vec3* bary = 0);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
//

template <> inline const Vec3 &pos<PS>(const Node *node) { return node->y; }
template <> inline const Vec3 &pos<WS>(const Node *node) { return node->x; }
template <> inline Vec3 &pos<PS>(Node *node) { return node->y; }
template <> inline Vec3 &pos<WS>(Node *node) { return node->x; }
template <> inline const Vec3 &pos<MS>(const Vert *vert) { return vert->u; }
template <> inline const Vec3 &pos<PS>(const Vert *vert) { return vert->node->y; }
template <> inline const Vec3 &pos<WS>(const Vert *vert) { return vert->node->x; }
template <> inline Vec3 &pos<MS>(Vert *vert) { return vert->u; }
template <> inline Vec3 &pos<PS>(Vert *vert) { return vert->node->y; }
template <> inline Vec3 &pos<WS>(Vert *vert) { return vert->node->x; }


#endif
