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

#include "remesh.hpp"

#include "../../Cloth.h"

#include "blockvectors.hpp"
#include "geometry.hpp"
#include "magic.hpp"
#include "referenceshape.hpp"
#include "util.hpp"

#include <assert.h>
#include <cstdlib>
#include <cstdio>
using namespace std;

// Helpers
template <class T> static void delete_all(const vector<T>& a) { for (size_t i = 0; i<a.size(); i++) delete a[i]; }
template <class T> static void remove_all(const vector<T>& a, Mesh& m) { for (size_t i = 0; i<a.size(); i++) m.remove(a[i]); }
template <class T> static void add_all(const vector<T>& a, Mesh& m) { for (size_t i = 0; i<a.size(); i++) m.add(a[i]); }
template <class T> static void include_all(const vector<T>& a, vector<T>& b) { for (size_t i = 0; i<a.size(); i++) include(a[i], b); }
template <class T> static void exclude_all(const vector<T>& a, vector<T>& b) { for (size_t i = 0; i<a.size(); i++) exclude(a[i], b); }

RemeshOp RemeshOp::inverse() const {
	RemeshOp iop;
	iop.added_verts = removed_verts;
	iop.removed_verts = added_verts;
	iop.added_nodes = removed_nodes;
	iop.removed_nodes = added_nodes;
	iop.added_edges = removed_edges;
	iop.removed_edges = added_edges;
	iop.added_faces = removed_faces;
	iop.removed_faces = added_faces;
	return iop;
}

void RemeshOp::cancel() {
	delete_all(added_verts);
	delete_all(added_nodes);
	delete_all(added_edges);
	delete_all(added_faces);
	added_edges.clear();
	added_faces.clear();
	added_nodes.clear();
	added_verts.clear();
	removed_edges.clear();
	removed_faces.clear();
	removed_nodes.clear();
	removed_verts.clear();
}

void RemeshOp::apply(Mesh &mesh) const {
	remove_all(removed_faces, mesh);
	remove_all(removed_edges, mesh);
	remove_all(removed_nodes, mesh);
	remove_all(removed_verts, mesh);
	add_all(added_verts, mesh);
	add_all(added_nodes, mesh);
	add_all(added_edges, mesh);
	add_all(added_faces, mesh);
	for (size_t f = 0; f < added_faces.size(); f++)
		compute_ms_data(added_faces[f]);
	for (size_t f = 0; f < added_faces.size(); f++)
		for (int i = 0; i < 3; i++)
			compute_ms_data(added_faces[f]->v[i]->node);
}

void RemeshOp::done() const {
	delete_all(removed_verts);
	delete_all(removed_nodes);
	delete_all(removed_edges);
	delete_all(removed_faces);
}

void RemeshOp::set_null(std::vector<Edge*>& v) {
	for (size_t i = 0; i<v.size(); i++)
		if (is_in(v[i], removed_edges))
			v[i] = 0;
}

void RemeshOp::update(std::vector<Face*>& v) {
	exclude_all(removed_faces, v);
	include_all(added_faces, v);
}

void RemeshOp::update(std::vector<Edge*>& v) {
	exclude_all(removed_edges, v);
	include_all(added_edges, v);
}

void RemeshOp::update(std::vector<Node*>& v) {
	exclude_all(removed_nodes, v);
	include_all(added_nodes, v);
}

ostream &operator<< (ostream &out, const RemeshOp &op) {
	out << "removed " << op.removed_verts << ", " << op.removed_nodes << ", "
		<< op.removed_edges << ", " << op.removed_faces << ", added "
		<< op.added_verts << ", " << op.added_nodes << ", " << op.added_edges
		<< ", " << op.added_faces;
	return out;
}

// Project new vertex onto material-space mesh

template<Space s>
Vec3 safe_normal(Face* face) {
	if (!face) return Vec3(0);
	const Vec3 a = pos<s>(face->v[1]) - pos<s>(face->v[0]);
	const Vec3 b = pos<s>(face->v[2]) - pos<s>(face->v[0]);
	Vec3 c = cross(a, b);
	if (fabs(dot(normalize(a), normalize(b))) > 0.9 && norm(c) < 1e-6)
		return Vec3(0); // unstable normal
	return normalize(c);
}

void project_vertex(Vert *vnew, Edge* edge, int s, double d) {
	Vert *v0 = edge_vert(edge, s, s), *v1 = edge_vert(edge, s, 1 - s);
	if (s != 0) d = 1 - d;
	Vec3 n = safe_normal<MS>(edge->adjf[0]) + safe_normal<MS>(edge->adjf[1]);
	Vec3 u = (1 - d)*v0->u + d*v1->u;
	//Vec3 v = (1 - d)*v0->v + d*v1->v; // NICK
	if (norm(n) > 0.2) {
		if (!edge->n[0]->mesh->ref->raycast(u, normalize(n))) {
			cout << "split:raycast failed" << endl;
			exit(1);
		}
	}
	vnew->u = u;
	//vnew->v = v; // NICK
	vnew->sizing = (1.0 - d)*v0->sizing + d*v1->sizing;
}

// Helpers for localopt

//void embedding_from_plasticity(const vector<Face*> &fs) {
//	vector<Node*> nodes;
//	vector<Face*> faces;
//	vector<Edge*> edges;
//	for (size_t f = 0; f < fs.size(); f++) {
//		const Face *face = fs[f];
//		include((Face*)face, faces);
//		for (int i = 0; i < 3; i++) {
//			include(face->v[i]->node, nodes);
//			const Edge *edge = face->adje[i];
//			include((Edge*)edge, edges);
//			int s = (face == edge->adjf[0] ? 0 : 1);
//			if (edge->adjf[1 - s]) { // include adjacent "flap"
//				include(edge->adjf[1 - s], faces);
//				include(edge_opp_vert(edge, 1 - s)->node, nodes);
//			}
//		}
//	}
//	for (size_t n = 0; n < nodes.size(); n++)
//		nodes[n]->y = nodes[n]->x;
//	vector<Constraint*> no_cons;
//	local_opt<PS>(nodes, faces, edges, no_cons);
//}

//void plasticity_from_embedding(const vector<Face*> &faces) {
//	for (size_t f = 0; f < faces.size(); f++) {
//		Face *face = faces[f];
//		if (face->material->plastic_flow == 0)
//			continue;
//		face->Sp_str = stretch_plasticity_from_embedding(face);
//	}
//	vector<Edge*> edges;
//	for (size_t f = 0; f < faces.size(); f++)
//		for (int i = 0; i < 3; i++)
//			include(faces[f]->adje[i], edges);
//	for (size_t e = 0; e < edges.size(); e++) {
//		Edge *edge = edges[e];
//		if (!edge->adjf[0] || !edge->adjf[1])
//			continue;
//		if (edge->adjf[0]->material->yield_curv < infinity
//			|| edge->adjf[1]->material->yield_curv < infinity)
//			edge->theta_ideal = dihedral_angle<PS>(edge);
//		else
//			edge->theta_ideal = dihedral_angle<MS>(edge);
//	}
//	for (size_t f = 0; f < faces.size(); f++)
//		recompute_Sp_bend(faces[f]);
//}

//struct PlasticityStash {
//	vector<Mat3x3> Sp_str;
//	vector<double> theta_ideal;
//	PlasticityStash(vector<Face*> &faces, vector<Edge*> &edges)
//		: Sp_str(faces.size()), theta_ideal(edges.size()) {
//		for (size_t f = 0; f < faces.size(); f++) {
//			Face *face = faces[f];
//			Sp_str[f] = face->Sp_str;
//			face->Sp_str = Mat3x3(1);
//		}
//		for (size_t e = 0; e < edges.size(); e++) {
//			Edge *edge = edges[e];
//			theta_ideal[e] = edge->theta_ideal;
//			edge->theta_ideal = dihedral_angle<MS>(edge);
//		}
//	}
//	void apply(vector<Face*> &faces, vector<Edge*> &edges) {
//		for (size_t f = 0; f < faces.size(); f++) {
//			faces[f]->Sp_str = Sp_str[f];
//		}
//		for (size_t e = 0; e < edges.size(); e++) {
//			edges[e]->theta_ideal = theta_ideal[e];
//		}
//	}
//};

//void optimize_node(Node *node) {
//	vector<Node*> nodes(1, node);
//	vector<Face*> faces;
//	for (size_t v = 0; v < node->verts.size(); v++)
//		append(faces, node->verts[v]->adjf);
//	vector<Edge*> edges;
//	for (size_t f = 0; f < faces.size(); f++) {
//		const Face *face = faces[f];
//		for (int i = 0; i < 3; i++)
//			include(face->adje[i], edges);
//	}
//	vector<Constraint*> no_cons;
//	//PlasticityStash stash(faces, edges);
//	//local_opt<PS>(nodes, faces, edges, no_cons); /////
//	//stash.apply(faces, edges);
//	activate_nodes(nodes);
//	deactivate_nodes(nodes);
//	compute_ws_data(nodes);
//	compute_ws_data(faces);
//	//local_opt<WS>(nodes, faces, edges, no_cons); /////
//}

//void local_pop_filter(const vector<Face*> &fs) {
//	if (!::magic.enable_localopt)
//		return;
//	vector<Node*> nodes;
//	vector<Face*> faces;
//	vector<Edge*> edges;
//	for (size_t f = 0; f < fs.size(); f++)
//		for (int i = 0; i < 3; i++)
//			include(fs[f]->v[i]->node, nodes);
//	for (size_t n = 0; n < nodes.size(); n++) {
//		for (size_t v = 0; v < nodes[n]->verts.size(); v++) {
//			const Vert *vert = nodes[n]->verts[v];
//			for (size_t f = 0; f < vert->adjf.size(); f++)
//				include(vert->adjf[f], faces);
//		}
//	}
//	for (size_t f = 0; f < faces.size(); f++)
//		for (int i = 0; i < 3; i++)
//			include(faces[f]->adje[i], edges);
//	vector<Constraint*> cons;
//	cons = proximity_constraints(sim.cloth_meshes, sim.obstacle_meshes,
//		sim.friction, sim.obs_friction, false, true);
//	for (int h = 0; h < (int)sim.handles.size(); h++)
//		append(cons, sim.handles[h]->get_constraints(sim.time));
//
//	/*if (sim.frame > 0) {
//	for (int i=0; i<nodes.size(); i++)
//	Annotation::add(nodes[i]);
//	for (int i=0; i<edges.size(); i++)
//	Annotation::add(edges[i], Vec3(0,0,1));
//	wait_key();}*/
//	local_opt<WS>(nodes, faces, edges, cons);
//	//if (sim.frame>0) wait_key();
//}

// This operation is needed when we introduce new nodes, i.e. in split_edge
// If both node are EOL AND are along a shared edge then the property is transfered
void transferEoL(Node *n, Node *n0, Node *n1) {
	if (n0->EoL && n1->EoL) {
		if (n0->cornerID >= 0 && n1->cornerID >= 0) return;
		if (n0->cornerID >= 0) {
			for (int j = 0; j < n0->cdEdges.size(); j++) {
				if (n1->cdEdges[0] == n0->cdEdges[j]) {
					n->EoL = true;
					n->cdEdges = n1->cdEdges;
					return;
				}
			}
		}
		else if (n1->cornerID >= 0) {
			for (int j = 0; j < n1->cdEdges.size(); j++) {
				if (n0->cdEdges[0] == n1->cdEdges[j]) {
					n->EoL = true;
					n->cdEdges = n0->cdEdges;
					return;
				}
			}
		}
		else if (n0->cdEdges[0] == n1->cdEdges[0]) {
			n->EoL = true;
			n->cdEdges = n0->cdEdges;
			return;
		}
	}
	return;
}

// The actual operations

RemeshOp split_edge(Edge* edge, double d) {
	Mesh& mesh = *edge->n[0]->mesh;
	RemeshOp op;
	Node *node0 = edge->n[0],
		*node1 = edge->n[1],
		*node = new Node((1 - d)*node0->y + d*node1->y,
		(1 - d)*node0->x + d*node1->x,
			(1 - d)*node0->v + d*node1->v,
			0, //node0->label & node1->label,
			node0->flag & node1->flag,
			false);
	node->acceleration = (1 - d)*node0->acceleration + d*node1->acceleration;
	transferEoL(node, node0, node1); // NICK
	if (node->EoL) node->cdEdges = node0->cdEdges;
	op.added_nodes.push_back(node);
	op.removed_edges.push_back(edge);
	op.added_edges.push_back(new Edge(node0, node, edge->theta_ideal,
		edge->preserve));
	op.added_edges.push_back(new Edge(node, node1, edge->theta_ideal,
		edge->preserve));
	Vert *vnew[2] = { NULL, NULL };
	for (int s = 0; s < 2; s++) {
		if (!edge->adjf[s])
			continue;
		Vert *v0 = edge_vert(edge, s, s),
			*v1 = edge_vert(edge, s, 1 - s),
			*v2 = edge_opp_vert(edge, s);
		if (s == 0 || is_seam_or_boundary(edge)) {
			vnew[s] = new Vert(Vec3(0), Vec3(0));
			project_vertex(vnew[s], edge, s, d);
			connect(vnew[s], node);
			op.added_verts.push_back(vnew[s]);
		}
		else
			vnew[s] = vnew[0];
		op.added_edges.push_back(new Edge(v2->node, node, 0, 0));
		Face *f = edge->adjf[s];
		op.removed_faces.push_back(f);
		Face* nf0 = new Face(v0, vnew[s], v2, f->Sp_str, f->Sp_bend, f->material, f->damage);
		Face* nf1 = new Face(vnew[s], v1, v2, f->Sp_str, f->Sp_bend, f->material, f->damage);
		if (min(aspect(nf0), aspect(nf1)) < 1e-3) {
			op.cancel();
			return op;
		}
		op.added_faces.push_back(nf0);
		op.added_faces.push_back(nf1);
	}
	//embedding_from_plasticity(op.removed_faces);
	op.apply(mesh);
	node->y = (1 - d)*node0->y + d*node1->y;
	//optimize_node(node);
	//plasticity_from_embedding(op.added_faces);
	//local_pop_filter(op.added_faces);
	return op;
}

RemeshOp collapse_edge(Edge* edge, int i) {
	/*if (is_seam_or_boundary(edge)) {
	Annotation::add(edge);
	cout << "collapse" << endl;
	cout << edge->n[i]->preserve << endl;
	wait_key();
	}*/
	Mesh& mesh = *edge->n[0]->mesh;
	RemeshOp op;
	Node *node0 = edge->n[i], *node1 = edge->n[1 - i];
	op.removed_nodes.push_back(node0);
	for (size_t e = 0; e < node0->adje.size(); e++) {
		Edge *edge1 = node0->adje[e];
		op.removed_edges.push_back(edge1);
		Node *node2 = (edge1->n[0] != node0) ? edge1->n[0] : edge1->n[1];
		if (node2 != node1 && !get_edge(node1, node2))
			op.added_edges.push_back(new Edge(node1, node2, edge1->theta_ideal,
				edge1->preserve));
		// Preserve in weird situations, also and issue with ArcSim?
		//if (node2 != node1 && (get_edge(node1, node2) != NULL && get_edge(node0, node2) != NULL)) {
		//	if (get_edge(node0, node2)->preserve) {
		//		get_edge(node1, node2)->preserve = true;
		//	}
		//}
	}
	for (int s = 0; s < 2; s++) {
		Vert *vert0 = edge_vert(edge, s, i), *vert1 = edge_vert(edge, s, 1 - i);
		if (!vert0 || (s == 1 && vert0 == edge_vert(edge, 0, i)))
			continue;
		op.removed_verts.push_back(vert0);
		for (size_t f = 0; f < vert0->adjf.size(); f++) {
			Face *face = vert0->adjf[f];
			op.removed_faces.push_back(face);
			if (!is_in(vert1, face->v)) {
				Vert *verts[3] = { face->v[0], face->v[1], face->v[2] };
				replace(vert0, vert1, verts);
				Face* new_face = new Face(verts[0], verts[1], verts[2],
					face->Sp_str, face->Sp_bend, face->material, face->damage);
				op.added_faces.push_back(new_face);
				// inversion test
				if (dot(normal<MS>(face), normal<MS>(new_face)) < 0) {
					op.cancel();
					return RemeshOp();
				}
				// degenerate test
				bool enforce = false;
				double asp_old = aspect(face), asp_new = aspect(new_face);
				if (asp_new < mesh.parent->remeshing.aspect_min / 4 &&
					asp_old >= mesh.parent->remeshing.aspect_min / 4 && !enforce) {
					op.cancel();
					return RemeshOp();
				}
			}
		}
	}
	//wait_key();
	//embedding_from_plasticity(op.removed_faces);
	op.apply(mesh);
	//plasticity_from_embedding(op.added_faces);
	//local_pop_filter(op.added_faces);
	return op;
}

RemeshOp flip_edge(Edge* edge) {
	RemeshOp op;
	Vert *vert0 = edge_vert(edge, 0, 0), *vert1 = edge_vert(edge, 1, 1),
		*vert2 = edge_opp_vert(edge, 0), *vert3 = edge_opp_vert(edge, 1);
	Face *face0 = edge->adjf[0], *face1 = edge->adjf[1];
	double A = area(face0), B = area(face1);
	Mat3x3 sp = (A * face0->Sp_str + B * face1->Sp_str) / (A + B);
	Mat3x3 sb = (A * face0->Sp_bend + B * face1->Sp_bend) / (A + B);
	double damage = (A * face0->damage + B * face1->damage) / (A + B);
	op.removed_edges.push_back(edge);
	op.added_edges.push_back(new Edge(vert2->node, vert3->node,
		-edge->theta_ideal, edge->preserve));
	op.removed_faces.push_back(face0);
	op.removed_faces.push_back(face1);
	op.added_faces.push_back(new Face(vert0, vert3, vert2, sp, sb, face0->material, damage));
	op.added_faces.push_back(new Face(vert1, vert2, vert3, sp, sb, face1->material, damage));
	//embedding_from_plasticity(op.removed_faces);
	op.apply(*edge->n[0]->mesh);
	//plasticity_from_embedding(op.added_faces);
	//local_pop_filter(op.added_faces);
	return op;
}

bool try_move_node(Node* node, Edge* edge, double d) {
	// TODO
	if (d<1e-6)
		return true;
	if (node->preserve)
		return false;
	vector<Face*> faces;
	for (size_t v = 0; v < node->verts.size(); v++)
		append(faces, node->verts[v]->adjf);
	//embedding_from_plasticity(faces);
	int v = edge->n[0] == node ? 0 : 1;
	if (v != 0) d = 1 - d;
	for (int i = 0; i<2; i++) {
		if (edge->adjf[i] && (i == 0 || is_seam_or_boundary(edge)))
			project_vertex(edge_vert(edge, i, v), edge, i, d);
	}
	for (size_t f = 0; f < faces.size(); f++)
		compute_ms_data(faces[f]);
	for (size_t f = 0; f < faces.size(); f++)
		for (int i = 0; i < 3; i++)
			compute_ms_data(faces[f]->v[i]->node);
	//optimize_node(node);
	//plasticity_from_embedding(faces);
	//local_pop_filter(faces);
	return true;
}
