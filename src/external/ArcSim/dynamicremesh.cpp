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

#include "dynamicremesh.hpp"
//#include "mesh.hpp"
#include "geometry.hpp"
//#include "sepstrength.hpp"
//#include "physics.hpp"
#include "magic.hpp"
#include "remesh.hpp"
//#include "simulation.hpp"
#include "tensormax.hpp"
//#include "timer.hpp"
//#include "cloth.hpp"
#include "util.hpp"
//#include "display.hpp"
#include "subset.hpp"
#include <algorithm>
#include <cstdlib>
#include <map>
using namespace std;

static const bool verbose = false;

void create_vert_sizing(vector<Vert*>& verts);

// The algorithm

void flip_edges(MeshSubset* subset, vector<Face*>& active_faces,
	vector<Edge*>* update_edges, vector<Face*>* update_faces);

bool split_worst_edge(MeshSubset* subset, const vector<Edge*>& edges);

bool improve_some_face(MeshSubset* subset, vector<Face*>& active);

void delete_spaced_out(Mesh& mesh);

void static_remesh(Mesh& mesh) {
	for (size_t i = 0; i<mesh.verts.size(); i++) {
		mesh.verts[i]->sizing = Mat3x3(1.f / sq(mesh.parent->remeshing.size_min));
	}
	while (split_worst_edge(0, mesh.edges));
	vector<Face*> active_faces = mesh.faces;
	while (improve_some_face(0, active_faces));
	compute_ms_data(mesh);
}

void dynamic_remesh(Mesh& mesh) {
	delete_spaced_out(mesh);
	create_vert_sizing(mesh.verts);
	vector<Face*> active_faces = mesh.faces;
	//cout << "before\n"; wait_key();
	flip_edges(0, active_faces, 0, 0);
	//cout << "post flip\n"; wait_key();
	while (split_worst_edge(0, mesh.edges));
	//cout << "post split\n"; wait_key();
	active_faces = mesh.faces;
	while (improve_some_face(0, active_faces));
	//cout << "post collapse\n"; wait_key();
	compute_ms_data(mesh);
}

//void dynamic_remesh(MeshSubset& subset, const map<Node*, Plane> &planes) {
//	vector<Vert*> verts = subset.get_verts();
//	create_vert_sizing(verts, planes);
//	vector<Face*> active_faces = subset.get_faces();
//	flip_edges(&subset, active_faces, 0, 0);
//	while (split_worst_edge(&subset, subset.get_edges()));
//	active_faces = subset.get_faces();
//	while (improve_some_face(&subset, active_faces));
//	active_faces = subset.get_faces();
//	compute_ms_data(subset.active_nodes);
//	compute_ms_data(active_faces);
//}

// Sizing

double angle(const Vec3 &n1, const Vec3 &n2) {
	return acos(clamp(dot(n1, n2), -1., 1.));
}

template <int n> Mat<n, n> sqrt(const Mat<n, n> &A) {
	Eig<n> eig = eigen_decomposition(A);
	for (int i = 0; i < n; i++)
		eig.l[i] = eig.l[i] >= 0 ? sqrt(eig.l[i]) : -sqrt(-eig.l[i]);
	return eig.Q*diag(eig.l)*eig.Q.t();
}

Mat2x2 perp(const Mat2x2 &A) {
	return Mat2x2(Vec2(A(1, 1), -A(1, 0)),
		Vec2(-A(0, 1), A(0, 0)));
}

template <Space s>
Mat3x3 deformation_gradient(const Face *face) {
	return derivative(pos<s>(face->v[0]->node), pos<s>(face->v[1]->node),
		pos<s>(face->v[2]->node), normal<s>(face), face) * face->Sp_str;
}

Mat2x2 compression_metric(const Face* face, const Mat3x3 &S2, const Mat3x2& UV, double c) {
	Mat3x3 F = deformation_gradient<WS>(face);
	Mat3x3 G = F.t() * F - Mat3x3(1);
	Mat2x2 e = UV.t() * G * UV;
	Mat2x2 e2 = UV.t() * G.t()*G * UV;
	Mat2x2 Sw2 = UV.t() * S2 * UV;

	Mat2x2 D = e2 - 4.0*sq(c)*perp(Sw2)*::magic.rib_stiffening;
	return get_positive(-e + sqrt(D)) / (2.0*sq(c));
}

Mat3x3 obstacle_metric(const Face *face, const map<Node*, Plane> &planes) {
	Mat3x3 o(0);
	for (int v = 0; v < 3; v++) {
		map<Node*, Plane>::const_iterator it = planes.find(face->v[v]->node);
		if (it == planes.end())
			continue;
		double h[3];
		for (int v1 = 0; v1 < 3; v1++)
			h[v1] = dot(face->v[v1]->node->x - it->second.x0, it->second.n);
		Vec3 dh = derivative(h[0], h[1], h[2], 0, face);

		o += outer(dh, dh) / sq(h[v]);
	}
	return o / 3.;
}

//Mat2x2 fracture_metric(Remeshing& remeshing, const Face* face) {
//	if (remeshing.refine_fracture == 0 || !sim.enabled[Simulation::Fracture])
//		return Mat2x2(0);
//	double fmax = 0;
//	for (int i = 0; i<3; i++) {
//		Face* f = adj_face(face, i);
//		if (f) {
//			Mat3x3 sig = compute_sigma(face);
//			Vec3 l = eigen_values(sig);
//			fmax = max(l[0], fmax);
//		}
//	}
//
//	double refine = remeshing.refine_fracture * 0.5 * face->material->toughness;
//	double ramp0 = refine * 0.5, ramp1 = refine / 0.5;
//	double v = clamp((fmax - ramp0) / (ramp1 - ramp0), 0.0, 1.0);
//	return Mat2x2(v);
//}

Mat3x3 compute_face_sizing(Remeshing& remeshing, const Face *face) {
	// project to in-plane 2D
	Mat3x3 base = local_base(normal<MS>(face));
	Mat3x2 UV(base.col(0), base.col(1));
	Mat2x3 UVt = UV.t();

	//Mat2x2 Sp = projected_curvature<PS>(face,UVt);
	Mat2x2 Sw1 = projected_curvature<WS>(face, UVt);
	Mat3x3 Sw2 = derivative(face->v[0]->node->n, face->v[1]->node->n,
		face->v[2]->node->n, Vec3(0), face);
	//Mat2x2 Mcurvp = !sim.enabled[Simulation::Plasticity] ? Mat2x2(0)
	//              : (Sp.t()*Sp)/sq(remeshing.refine_angle);
	Mat2x2 Mcurvw1 = (Sw1.t()*Sw1) / sq(remeshing.refine_angle);
	Mat2x2 Mcurvw2 = UVt * (Sw2.t()*Sw2) * UV / sq(remeshing.refine_angle);
	Mat3x3 V = derivative(face->v[0]->node->v, face->v[1]->node->v,
		face->v[2]->node->v, Vec3(0), face);
	Mat2x2 Mvel = UVt * (V.t()*V) * UV / sq(remeshing.refine_velocity);

	Mat2x2 Mcomp = compression_metric(face, Sw2.t()*Sw2, UV, remeshing.refine_compression);
	//Mat2x2 Mobs = (planes.empty()) ? Mat2x2(0) : UVt*obstacle_metric(face, planes)*UV;
	Mat2x2 Mobs = Mat2x2(0);

	//Mat2x2 Mfrac = fracture_metric(remeshing, face) / sq(remeshing.size_min);
	Mat2x2 Mfrac = Mat2x2(0);
	//Mat2x2 Mfrac = Mat2x2(1.f) / sq(remeshing.size_min);

	vector<Mat2x2> Ms(6);
	//Ms[0] = Mcurvp;
	Ms[0] = Mcurvw1;
	Ms[1] = Mcurvw2;
	Ms[2] = Mvel;
	Ms[3] = Mcomp;
	Ms[4] = Mobs;
	Ms[5] = Mfrac;

	Mat2x2 s = ::magic.combine_tensors ? tensor_max(Ms)
		: Ms[0] + Ms[1] + Ms[2] + Ms[3] + Ms[4] + Ms[5] + Ms[6];
	//s = Ms[1];
	//cout << s << endl;

	// specific resolution request ?
	for (int i = 0; i<3; i++) {
		//if (face->v[0]->node->flag & Node::FlagResolveUni)
			//s = Mat2x2(1.f / sq(remeshing.size_uniform));
		if (face->v[0]->node->flag & Node::FlagResolveMax)
			s = Mat2x2(1.f / sq(remeshing.size_min));
	}

	//if (debug) {
	//	cout << "curvw1 " << norm_F(Mcurvw1) << endl;
	//	cout << "curvw2 " << norm_F(Mcurvw2) << endl;
	//	cout << "vel " << norm_F(Mvel) << endl;
	//	cout << "comp " << norm_F(Mcomp) << endl;
	//	cout << "obs " << norm_F(Mobs) << endl;
	//	cout << "frac " << norm_F(Mfrac) << endl;
	//	cout << "total " << norm_F(s) << endl;
	//}

	Eig<2> eig = eigen_decomposition(s);
	for (int i = 0; i < 2; i++) {
		eig.l[i] = clamp(eig.l[i],
			1.f / sq(remeshing.size_max),
			1.f / sq(remeshing.size_min));
	}
	double lmax = max(eig.l[0], eig.l[1]);
	double lmin = lmax*sq(remeshing.aspect_min);
	for (int i = 0; i < 2; i++)
		if (eig.l[i] < lmin)
			eig.l[i] = lmin;
	s = eig.Q*diag(eig.l)*eig.Q.t();
	//s(0, 0) = 1;
	//s(0, 1) = 0;
	//s(1, 0) = 0;
	//s(1, 1) = 1;
	//cout << s << endl;
	return UV * s * UVt; // reproject to 3D
}

void create_vert_sizing(vector<Vert*>& verts) {
	Remeshing& remeshing = verts[0]->node->mesh->parent->remeshing;
	map<Face*, Mat3x3> face_sizing;
	for (size_t i = 0; i<verts.size(); i++) {
		Mat3x3 sizing(0);
		Vert* vert = verts[i];
		double wsum = 0;
		for (size_t f = 0; f<vert->adjf.size(); f++) {
			Face *face = vert->adjf[f];
			if (!face) continue;
			if (face_sizing.find(face) == face_sizing.end())
				face_sizing[face] = compute_face_sizing(remeshing, face);
			sizing += face->a * face_sizing[face];
			wsum += face->a;
		}
		vert->sizing = sizing / wsum;
	}
}

double edge_metric(const Vert *vert0, const Vert *vert1) {
	if (!vert0 || !vert1)
		return 0;
	Vec3 du = vert0->u - vert1->u;
	return sqrt((dot(du, vert0->sizing * du) + dot(du, vert1->sizing * du)) / 2.);
}

double edge_metric(const Edge *edge) {
	return max(edge_metric(edge_vert(edge, 0, 0), edge_vert(edge, 0, 1)),
		edge_metric(edge_vert(edge, 1, 0), edge_vert(edge, 1, 1)));
}

// Fixing-upping
vector<Edge*> find_edges_to_flip(vector<Face*>& active_faces);
vector<Edge*> independent_edges(const vector<Edge*> &edges);

bool flip_some_edges(MeshSubset* subset, vector<Face*>& active_faces,
	vector<Edge*>* update_edges, vector<Face*>* update_faces) {
	static int n_edges_prev = 0;
	vector<Edge*> edges = independent_edges(find_edges_to_flip(active_faces));
	if ((int)edges.size() == n_edges_prev) // probably infinite loop
		return false;

	bool did_flip = false;
	n_edges_prev = edges.size();
	for (size_t e = 0; e < edges.size(); e++) {
		RemeshOp op = flip_edge(edges[e]);
		if (op.empty()) continue;

		did_flip = true;
		if (subset)
			op.update(subset->active_nodes);
		if (update_edges)
			op.set_null(*update_edges);
		if (update_faces)
			op.update(*update_faces);
		op.update(active_faces);
		op.done();
	}
	return did_flip;
}

void flip_edges(MeshSubset* subset, vector<Face*>& active_faces,
	vector<Edge*>* update_edges, vector<Face*>* update_faces) {
	int N = 3 * active_faces.size();
	for (int i = 0; i < N; i++) {// don't loop without bound
		if (!flip_some_edges(subset, active_faces, update_edges, update_faces))
			return;
	}
}


bool should_flip(const Edge *edge);

vector<Edge*> find_edges_to_flip(vector<Face*>& active_faces) {
	vector<Edge*> edges;
	for (size_t i = 0; i < active_faces.size(); i++)
		for (int j = 0; j<3; j++)
			include(active_faces[i]->adje[j], edges);

	vector<Edge*> fedges;
	for (size_t e = 0; e < edges.size(); e++) {
		Edge *edge = edges[e];
		if (is_seam_or_boundary(edge) || edge->preserve != 0 || !should_flip(edge))
			continue;
		fedges.push_back(edge);
	}
	return fedges;
}

bool independent(const Edge *edge, const vector<Edge*> &edges) {
	for (int i = 0; i < (int)edges.size(); i++) {
		Edge *edge1 = edges[i];
		if (edge->n[0] == edge1->n[0] || edge->n[0] == edge1->n[1]
			|| edge->n[1] == edge1->n[0] || edge->n[1] == edge1->n[1])
			return false;
	}
	return true;
}

vector<Edge*> independent_edges(const vector<Edge*> &edges) {
	vector<Edge*> iedges;
	for (int e = 0; e < (int)edges.size(); e++)
		if (independent(edges[e], iedges))
			iedges.push_back(edges[e]);
	return iedges;
}

double cross(const Vec2 &u, const Vec2 &v) { return u[0] * v[1] - u[1] * v[0]; }

// from Bossen and Heckbert 1996
inline bool should_flip2(const Vec2& x, const Vec2& y, const Vec2& z, const Vec2& w,
	const Mat2x2& M) {
	double area0 = fabs(wedge(z - y, x - y)), area1 = fabs(wedge(x - w, z - w));
	return area0*dot(x - w, M*(z - w)) + dot(z - y, M*(x - y))*area1
		< -::magic.edge_flip_threshold*(area0 + area1);
}

bool should_flip(const Edge *edge) {
	const double max_angle = 40.0*M_PI / 180.0;

	const Vert *vert0 = edge_vert(edge, 0, 0), *vert1 = edge_vert(edge, 0, 1),
		*vert2 = edge_opp_vert(edge, 0), *vert3 = edge_opp_vert(edge, 1);
	Vec3 x = vert0->u, z = vert1->u, w = vert2->u, y = vert3->u;

	// don't flip if high angles are involved
	if (fabs(dihedral_angle<WS>(edge)) >= max_angle) return false;

	const Vec3& n0 = normalize(cross(y - x, w - x)), n1 = normalize(cross(w - z, y - z));
	if (fabs(dihedral_angle(w, y, n0, n1)) >= max_angle) return false;

	// don't flip if mean normal is inverted
	if (dot(n0 + n1, normal<MS>(edge->adjf[0]) + normal<MS>(edge->adjf[1])) <= 0) return false;

	// project onto a 2D plane which conserves diagonal lengths 
	Vec3 u = normalize(x - z), v = normalize((w - y) - dot(w - y, u)*u);
	Mat3x2 A(u, v);
	Mat2x3 At = A.t();

	Mat3x3 M = (vert0->sizing + vert1->sizing + vert2->sizing + vert3->sizing) / 4.0;

	Mat2x2 Mr = At*M*A;
	return should_flip2(At*x, At*y, At*z, At*w, Mr);
}

// Splitting

vector<Edge*> find_bad_edges(const vector<Edge*>& edges);

Vert *adjacent_vert(const Node *node, const Vert *vert);

bool split_worst_edge(MeshSubset* subset, const vector<Edge*>& edges) {
	vector<Edge*> bad_edges = find_bad_edges(edges);
	for (size_t e = 0; e < bad_edges.size(); e++) {
		Edge *edge = bad_edges[e];
		if (!edge) continue;
		Node *node0 = edge->n[0], *node1 = edge->n[1];
		RemeshOp op = split_edge(edge, 0.5);
		for (size_t v = 0; v < op.added_verts.size(); v++) {
			Vert *vertnew = op.added_verts[v];
			Vert *v0 = adjacent_vert(node0, vertnew),
				*v1 = adjacent_vert(node1, vertnew);
			vertnew->sizing = 0.5 * (v0->sizing + v1->sizing);
		}
		if (subset)
			op.update(subset->active_nodes);
		op.set_null(bad_edges);
		op.done();
		if (verbose)
			cout << "Split " << node0 << " and " << node1 << endl;
		vector<Face*> active = op.added_faces;
		flip_edges(subset, active, &bad_edges, 0);
	}
	return !bad_edges.empty();
}

// don't use edge pointer as secondary sort key, otherwise not reproducible
struct Deterministic_sort {
	inline bool operator()(const std::pair<double, Edge*> &left, const std::pair<double, Edge*> &right) {
		return left.first < right.first;
	}
} deterministic_sort;

vector<Edge*> find_bad_edges(const vector<Edge*>& edges) {
	vector< pair<double, Edge*> > edgems;
	for (size_t e = 0; e < edges.size(); e++) {
		double m = edge_metric(edges[e]);
		if (m > 1)
			edgems.push_back(make_pair(m, edges[e]));
	}
	sort(edgems.begin(), edgems.end(), deterministic_sort);
	vector<Edge*> bad_edges(edgems.size());
	for (int e = 0; e < (int)edgems.size(); e++)
		bad_edges[e] = edgems[edgems.size() - e - 1].second;
	return bad_edges;
}

Vert *adjacent_vert(const Node *node, const Vert *vert) {
	const Edge *edge = get_edge(node, vert->node);
	for (int i = 0; i < 2; i++)
		for (int s = 0; s < 2; s++)
			if (edge_vert(edge, s, i) == vert)
				return edge_vert(edge, s, 1 - i);
	return NULL;
}

// Collapsing

vector<int> sort_edges_by_length(const Face *face);

RemeshOp try_edge_collapse(Edge *edge, int which);

bool improve_some_face(MeshSubset* subset, vector<Face*>& active) {
	for (size_t f = 0; f<active.size(); f++) {
		Face *face = active[f];
		for (int e = 0; e < 3; e++) {
			Edge *edge = face->adje[e];
			RemeshOp op;
			if (op.empty()) op = try_edge_collapse(edge, 0);
			if (op.empty()) op = try_edge_collapse(edge, 1);
			if (op.empty()) continue;
			op.update(active);
			if (subset)
				op.update(subset->active_nodes);
			op.done();
			vector<Face*> fix_active = op.added_faces;
			flip_edges(subset, fix_active, 0, &active);
			return true;
		}
		remove(f--, active);
	}
	return false;
}

bool has_labeled_edges(const Node *node);
bool can_collapse(Remeshing& rem, const Edge *edge, int which);
bool any_nearly_invalid(const vector<Edge*> edges) {
	for (int i = 0; i < (int)edges.size(); i++)
		if (edge_metric(edges[i]) > 0.9) return true;
	return false;
}

RemeshOp try_edge_collapse(Edge *edge, int which) {
	Node *node0 = edge->n[which], *node1 = edge->n[1 - which];
	if (node0->preserve || node0->label != 0
		|| (is_seam_or_boundary(node0) && !is_seam_or_boundary(edge))
		|| (has_labeled_edges(node0) && !edge->preserve))
		return RemeshOp();
	if (!can_collapse(node0->mesh->parent->remeshing, edge, which))
		return RemeshOp();
	RemeshOp op = collapse_edge(edge, which);
	if (op.empty())
		return op;
	if (verbose)
		cout << "Collapsed " << node0 << " into " << node1 << endl;
	return op;
}

bool has_labeled_edges(const Node *node) {
	for (int e = 0; e < (int)node->adje.size(); e++)
		if (node->adje[e]->preserve)
			return true;
	return false;
}

bool can_collapse(Remeshing& remeshing, const Edge *edge, int i) {
	for (int s = 0; s < 2; s++) {
		const Vert *vert0 = edge_vert(edge, s, i), *vert1 = edge_vert(edge, s, 1 - i);
		if (!vert0 || (s == 1 && vert0 == edge_vert(edge, 0, i)))
			continue;
		for (int f = 0; f < (int)vert0->adjf.size(); f++) {
			const Face *face = vert0->adjf[f];
			if (is_in(vert1, face->v))
				continue;
			const Vert *vs[3] = { face->v[0], face->v[1], face->v[2] };
			double a0 = norm(cross(vs[1]->u - vs[0]->u, vs[2]->u - vs[0]->u)) / 2;
			replace(vert0, vert1, vs);
			double a = norm(cross(vs[1]->u - vs[0]->u, vs[2]->u - vs[0]->u)) / 2;
			double asp = aspect(vs[0]->u, vs[1]->u, vs[2]->u);

			if ((a < a0 && a < 0.1*sqr(remeshing.size_min)) || asp < remeshing.aspect_min)
				return false;
			for (int e = 0; e < 3; e++)
				if (vs[e] != vert1 && edge_metric(vs[NEXT(e)], vs[PREV(e)]) > 0.9) {
					return false;
				}
		}
	}
	return true;
}

void delete_spaced_out(Mesh& mesh) {
	RemeshOp op;
	for (size_t i = 0; i<mesh.verts.size(); i++)
		if (norm2(mesh.verts[i]->node->x) > 1e6)
			op.removed_verts.push_back(mesh.verts[i]);

	if (op.removed_verts.empty())
		return;

	bool grow = true;
	while (grow) {
		grow = false;
		for (size_t i = 0; i < op.removed_verts.size(); i++) {
			Vert* v = op.removed_verts[i];
			for (size_t j = 0; j < v->adjf.size(); j++) {
				Face* f = v->adjf[j];
				if (find(f, op.removed_faces) < 0) {
					grow = true;
					op.removed_faces.push_back(f);
					for (int k = 0; k<3; k++) {
						include(f->adje[k], op.removed_edges);
						include(f->v[k], op.removed_verts);
						include(f->v[k]->node, op.removed_nodes);
					}
				}
			}
		}
	}
	cout << "deleted " << op.removed_nodes.size() << " spaced out nodes" << endl;
	op.apply(mesh);
	op.done();
}