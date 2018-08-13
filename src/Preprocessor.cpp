#include "Preprocessor.h"

#include "remeshExtension.h"
#include "matlabOutputs.h"
#include "conversions.h"

#include "external\ArcSim\geometry.hpp"
#include "external\ArcSim\subset.hpp"

#include <stdlib.h>

using namespace std;
//using namespace Eigen;

double thresh = 0.01; // TODO:: Move
double boundary = 0.01; // TODO:: Move;

Vert *adjacent_vert(const Node *node, const Vert *vert);

double linepoint(const Vec3 &A, const Vec3 &B, const Vec3 &P)
{
	// Line: A -> B
	// Point: P
	Vec3 AP = P - A;
	Vec3 AB = B - A;
	double ab2 = dot(AB,AB);
	double apab = dot(AP,AB);
	return apab / ab2;
}

void markWasEOL(Mesh& mesh) {
	for (int n = 0; n < mesh.nodes.size(); n++) {
		if (mesh.nodes[n]->EoL) {
			mesh.nodes[n]->EoL_state = Node::WasEOL;
		}
	}
}

void addGeometry(Mesh& mesh, const vector<shared_ptr<btc::Collision> > cls)
{
	for (int i = 0; i < cls.size(); i++) {
		// EOL nodes will be detected here
		// If they aren't found then the EOL node has been lifted off and we need to take note
		if (cls[i]->count1 == 3 && cls[i]->count2 == 1) {
			Node* node = mesh.nodes[cls[i]->verts2(0)];
			if (node->EoL) {
				node->EoL_state = Node::IsEOL;
			}
		}
		else if (cls[i]->count1 == 1 && cls[i]->count2 == 3) {
			if ((mesh.nodes[cls[i]->verts2(0)]->EoL && mesh.nodes[cls[i]->verts2(0)]->cornerID == cls[i]->verts1(0)) ||
				(mesh.nodes[cls[i]->verts2(1)]->EoL && mesh.nodes[cls[i]->verts2(1)]->cornerID == cls[i]->verts1(0)) ||
				(mesh.nodes[cls[i]->verts2(2)]->EoL && mesh.nodes[cls[i]->verts2(2)]->cornerID == cls[i]->verts1(0))) continue;
				
			// We don't add new points throughout this loop so this is safe
			Vert *v0 = mesh.verts[cls[i]->verts2(0)],
				*v1 = mesh.verts[cls[i]->verts2(1)],
				*v2 = mesh.verts[cls[i]->verts2(2)];
			//Face *f0 = mesh.faces[cls[i]->tri2];
				
			// TODO:: This may be overkill if the CD weights match the mesh barycoords
			double xX = cls[i]->weights2(0) * v0->u[0] +
				cls[i]->weights2(1) * v1->u[0] +
				cls[i]->weights2(2) * v2->u[0];
			double yX = cls[i]->weights2(0) * v0->u[1] +
				cls[i]->weights2(1) * v1->u[1] +
				cls[i]->weights2(2) * v2->u[1];

			// Boundary
			// TODO:: We want to use a generic check instead of a 0 to 1 range
			if (xX - boundary <= 0.0 || xX + boundary >= 1.0 ||
				yX - boundary <= 0.0 || yX + boundary >= 1.0) {
				continue;
			}
				
			// Faces CAN be added and deleted throughout this loop so we can't just use the CD returned tri
			Face *f0 = get_enclosing_face(mesh, Vec2(xX, yX));
			Vec3 bary = get_barycentric_coords(Vec2(xX, yX), f0);

			bool use_edge = false;
			for (int j = 0; j < 3; j++) {
				if (bary[j] < 1e-3) use_edge = true; // TODO:: No magic
			}

			// If this point is on a cloth edge, we should edge split instead
			if (use_edge) {
				double least = 1.0;
				int which = -1;
				for (int j = 0; j < 3; j++) {
					if (bary[j] < least) {
						least = bary[j];
						which = j;
					}
				}
				Edge *e0 = get_opp_edge(f0, f0->v[which]->node);
				double d;
				if (e0->n[0] == f0->v[0]->node) {
					d = 1.0 - bary[0];
				}
				else if (e0->n[0] == f0->v[1]->node) {
					d = 1.0 - bary[1];
				}
				else {
					d = 1.0 - bary[2];
				}
				Node *node0 = e0->n[0], *node1 = e0->n[1];
				RemeshOp op = split_edgeForced(e0, d, -1);
				for (size_t v = 0; v < op.added_verts.size(); v++) {
					Vert *vertnew = op.added_verts[v];
					Vert *v0 = adjacent_vert(node0, vertnew),
						*v1 = adjacent_vert(node1, vertnew);
					vertnew->sizing = 0.5 * (v0->sizing + v1->sizing);
				}
				op.done();
			}
			else {
				RemeshOp op = split_face(f0, bary);
				for (size_t v = 0; v < op.added_verts.size(); v++) {
					Vert *vertnew = op.added_verts[v];
					vertnew->sizing = (v0->sizing + v1->sizing + v2->sizing) / 3.0;
				}
				op.done();
			}

			Node *n = mesh.nodes.back();
			n->EoL = true;
			n->EoL_state = Node::NewEOL;
			n->preserve = true;
			n->cornerID = cls[i]->verts1(0);
			n->cdEdges = cls[i]->edge1;
			n->x = e2v(cls[i]->pos1_); // We want to offset the node slightly inside the obect so it gets detected by the CD until it moves out on it's own
		}
		if (cls[i]->count1 == 2 && cls[i]->count2 == 2) {
			// TODO:: Is this enough of a check?
			if (mesh.nodes[cls[i]->verts2(0)]->EoL ||
				mesh.nodes[cls[i]->verts2(1)]->EoL) continue;

			// The verts won't change since we only ever add verts
			Edge *e0 = get_edge(mesh.nodes[cls[i]->verts2(0)], mesh.nodes[cls[i]->verts2(1)]);
			double d;

			// If a previous collisions has already split this edge this get_edge will return NULL and we have to do a little more work to find it
			if (e0 == NULL) {
				Vert *v0 = mesh.verts[cls[i]->verts2(0)],
					*v1 = mesh.verts[cls[i]->verts2(1)];

				double xX = cls[i]->weights2(0) * v0->u[0] +
					cls[i]->weights2(1) * v1->u[0];
				double yX = cls[i]->weights2(0) * v0->u[1] +
					cls[i]->weights2(1) * v1->u[1];

				// Boundary
				// TODO:: We want to use a generic check instead of a 0 to 1 range
				if (xX - boundary <= 0.0 || xX + boundary >= 1.0 ||
					yX - boundary <= 0.0 || yX + boundary >= 1.0) {
					continue;
				}

				// A barycentric tests largest two values should be the nodes of the new edge we need to split
				Face *f0 = get_enclosing_face(mesh, Vec2(xX, yX));
				Vec3 bary = get_barycentric_coords(Vec2(xX, yX), f0);
				double least = 1.0;
				int which = -1;
				for (int j = 0; j < 3; j++) {
					if (bary[j] < least) {
						least = bary[j];
						which = j;
					}
				}
				e0 = get_opp_edge(f0, f0->v[which]->node);
				if (e0->n[0] == f0->v[0]->node) {
					d = 1.0 - bary[0];
				}
				else if (e0->n[0] == f0->v[1]->node) {
					d = 1.0 - bary[1];
				}
				else {
					d = 1.0 - bary[2];
				}
			}
			else {
				// ArcSim splits an edge using the weight from n[1]
				// Double check in case the CD does not index the same way
				if (e0->n[0]->index == cls[i]->verts2(1)) {
					d = cls[i]->weights2(0);
				}
				else {
					d = cls[i]->weights2(1);
				}
			}

			Node *node0 = e0->n[0], *node1 = e0->n[1];
			RemeshOp op = split_edgeForced(e0, d, -1);
			for (size_t v = 0; v < op.added_verts.size(); v++) {
				Vert *vertnew = op.added_verts[v];
				Vert *v0 = adjacent_vert(node0, vertnew),
					*v1 = adjacent_vert(node1, vertnew);
				vertnew->sizing = 0.5 * (v0->sizing + v1->sizing);
			}
			op.done();

			Node *n = mesh.nodes.back();
			n->EoL = true;
			n->EoL_state = Node::NewEOL;
			n->cdEdges = cls[i]->edge1;
			n->x = e2v(cls[i]->pos1_); // We want to offset the node slightly inside the obect so it gets detected by the CD until it moves out on it's own
		}
	}
}

void revertWasEOL(Mesh& mesh)
{
	// If something is still marked as WasEOL then it has lifted off
	for (int n = 0; n < mesh.nodes.size(); n++) {
		Node* node = mesh.nodes[n];
		if (node->EoL_state == Node::WasEOL) {
			node->EoL = false;
			node->preserve = false;
			node->cornerID = -1;
			node->cdEdges.clear();
		}
		// Boundary
		// TODO:: We want to use a generic check instead of a 0 to 1 range
		if (node->verts[0]->u[0] - boundary <= 0.0 || node->verts[0]->u[0] + boundary >= 1.0 ||
			node->verts[0]->u[1] - boundary <= 0.0 || node->verts[0]->u[1] + boundary >= 1.0) {
			node->EoL = false;
			node->preserve = false;
			node->cornerID = -1;
			node->cdEdges.clear();
		}
	}
}

void markPreserve(Mesh& mesh)
{
	for (int i = 0; i < mesh.edges.size(); i++) {
		Edge *e = mesh.edges[i];
		e->preserve = false;
		Node *n0 = e->n[0],
			*n1 = e->n[1];
		if (n0->EoL && n1->EoL) {
			// We never connect corners together
			// TODO:: Is this correct assumption?
			if (n0->cornerID >= 0 && n1->cornerID >= 0) continue;

			// We check the corner to edge, and edge to edge cases
			// A corner will have a list of cdEdges, while an edge is garaunteed to only have one cdEdge
			bool match = false;
			if (n0->cornerID >= 0) {
				for (int j = 0; j < n0->cdEdges.size(); j++) {
					if (n1->cdEdges[0] == n0->cdEdges[j]) {
						match = true;
						break;
					}
				}
			}
			else if (n1->cornerID >= 0) {
				for (int j = 0; j < n1->cdEdges.size(); j++) {
					if (n0->cdEdges[0] == n1->cdEdges[j]) {
						match = true;
						break;
					}
				}
			}
			else if (n0->cdEdges[0] == n1->cdEdges[0]) {
				match = true;
			}

			// If two EoL nodes share an edge AND they share a cdEdge ID
			// then the connected edge corresponds to object geometry and is preserved
			if (match) e->preserve = true;
		}
	}
}

double edge_metric(const Vert *vert0, const Vert *vert1);
bool can_collapseForced(const Edge *edge, int i) {
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

			double sz = 0.1*sqr(0.5);

			if ((a < a0 && a < 1e-6) || asp < 1e-6) {
				bool get_later = false;
				for (int ne = 0; ne < 3; ne++) {
					int nep;
					ne == 2 ? nep = 0 : nep = ne + 1;
					if (unsigned_vv_distance(vs[ne]->u, vs[nep]->u) < thresh) {
						get_later = true;
						break;
					}
				}
				if (!get_later) return false;
			}
			//for (int e = 0; e < 3; e++)
			//	if (vs[e] != vert1 && edge_metric(vs[NEXT(e)], vs[PREV(e)]) > 0.9) {
			//		return false;
			//	}
		}
	}
	return true;
}

// For very illconditioned geometry we can potentially have two preserved edges forming a triangle
// We want to collapse this triangle and just make the non preserved edge the new single preserved edge
// This triangle was practically a line it was so thin that nothing is really being altered dramatically
void pass_collapse(RemeshOp op, Node *n)
{
	for (int e = 0; e < op.removed_edges.size(); e++) {
		if (op.removed_edges[e]->preserve) {
			if (op.removed_edges[e]->n[0] == n || op.removed_edges[e]->n[1] == n) continue;
			Edge *ep = get_edge(n, op.removed_edges[e]->n[0]);
			if(ep == NULL) ep = get_edge(n, op.removed_edges[e]->n[1]);
			if (ep != NULL) ep->preserve = true;
		}
	}
}

bool collapse_conformal(Mesh &mesh, bool &allclear)
{
	for (int i = 0; i < mesh.edges.size(); i++) {
		Edge *e = mesh.edges[i];
		if (e->preserve) {
			if (edge_length(e) < thresh) {
				allclear = false;
				RemeshOp op;
				Node *n0 = e->n[0],
					*n1 = e->n[1];
				if (n0->index == 607 || n1->index == 607) {
					cout << endl;
				}
				if (is_seam_or_boundary(n1) || n1->cornerID >= 0) {
					if (!can_collapseForced(e, 0)) continue;
					op = collapse_edgeForced(e, 0);
					if (op.empty()) continue;
					pass_collapse(op, n1);
					op.done();
					return true;
				}
				else if (is_seam_or_boundary(n0) || n0->cornerID >= 0) {
					if (!can_collapseForced(e, 1)) continue;
					op = collapse_edgeForced(e, 1);
					if (op.empty()) continue;
					pass_collapse(op, n0);
					op.done();
					return true;
				}
				else {
					if (n0->verts[0]->adjf.size() <= n1->verts[0]->adjf.size()) {
						if (!can_collapseForced(e, 1)) op = collapse_edgeForced(e, 1);
						if (op.empty()) {
							if (!can_collapseForced(e, 0)) continue;
							op = collapse_edgeForced(e, 0);
							if (op.empty()) {
								continue;
							}
							pass_collapse(op, n1);
						}
						else pass_collapse(op, n0);
					}
					else {
						if (!can_collapseForced(e, 0)) op = collapse_edgeForced(e, 0);
						if (op.empty()) {
							if (!can_collapseForced(e, 1)) continue;
							op = collapse_edgeForced(e, 1);
							if (op.empty()) {
								continue;
							}
							pass_collapse(op, n0);
						}
						else pass_collapse(op, n1);
					}
					op.done();
					return true;
				}
			}
		}
	}
	return false;
}

bool collapse_nonconformal(Mesh &mesh, bool &allclear)
{
	for (int i = 0; i < mesh.nodes.size(); i++) {
		Node *n = mesh.nodes[i];
		if (n->EoL) {
			Vert *v = n->verts[0];
			for (int f = 0; f < v->adjf.size(); f++) {
				for (int e = 0; e < 3; e++) {
					Edge *e0 = v->adjf[f]->adje[e];
					if (!e0->preserve && edge_length(e0) < thresh) {
						Node *n0 = e0->n[0],
							*n1 = e0->n[1];
						if (n0->EoL && n1->EoL) continue;
						// Don't deal with edges between boundary and inside
						if (!(
							(is_seam_or_boundary(n0) && is_seam_or_boundary(n1)) ||
							(!is_seam_or_boundary(n0) && !is_seam_or_boundary(n1))
							)) continue;
						allclear = false;
						RemeshOp op;
						if (n0->EoL) {
							if (!can_collapseForced(e0, 1)) continue;
							op = collapse_edgeForced(e0, 1);
							if (op.empty()) continue;
							op.done();
							return true;
						}
						else if (n1->EoL) {
							if (!can_collapseForced(e0, 0)) continue;
							op = collapse_edgeForced(e0, 0);
							if (op.empty()) continue;
							op.done();
							return true;
						}
						else {
							if (!n0->preserve) {
								if (can_collapseForced(e0, 0)) op = collapse_edgeForced(e0, 0);
								if (op.empty()) {
									if (!can_collapseForced(e0, 1)) continue;
									op = collapse_edgeForced(e0, 1);
									if (op.empty()) {
										continue;
									}
								}
							}
							else if (!n1->preserve) {
								if (can_collapseForced(e0, 1)) op = collapse_edgeForced(e0, 1);
								if (op.empty()) {
									if (!can_collapseForced(e0, 0)) continue;
									op = collapse_edgeForced(e0, 0);
									if (op.empty()) {
										continue;
									}
								}
							}
							op.done();
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}

//bool collapse_close(Mesh &mesh)
//{
//	for (int i = 0; i < mesh.nodes.size(); i++) {
//		Node *n0 = mesh.nodes[i];
//		if (n0->EoL) {
//			for (int e = 0; e < n0->adje.size(); e++) {
//				Edge *e0 = n0->adje[e];
//				Node *n1 = other_node(e0, n0);
//				if (!n1->EoL && !is_seam_or_boundary(n1) && edge_length(e0) < thresh) {
//					RemeshOp op;
//					if (n0 == e0->n[0]) op = collapse_edgeForced(e0, 1);
//					else op = collapse_edgeForced(e0, 0);
//					if (op.empty()) continue;
//					op.done();
//					return true;
//				}
//			}
//		}
//	}
//	return false;
//}

int conformalCount(Face *f)
{
	int count = 0;
	for (int e = 0; e < 3; e++) {
		if (f->adje[e]->preserve) count++;
	}
	return count;
}

Node *single_eol_from_face(Face *f)
{
	for (int v = 0; v < 3; v++) {
		if (f->v[v]->node->EoL) return f->v[v]->node;
	}
	return NULL;
}

Edge *single_conformal_edge_from_face(Face *f)
{
	for (int e = 0; e < 3; e++) {
		if (f->adje[e]->preserve) return f->adje[e];
	}
	return NULL;
}

Edge *single_nonconformal_edge_from_face(Face *f)
{
	for (int e = 0; e < 3; e++) {
		if (!f->adje[e]->preserve) return f->adje[e];
	}
	return NULL;
}

double face_altitude(Edge* edge, Face* face) {
	return (2 * area(face)) / edge_length(edge);
}


// TODO:: Better metric than face altitude?
bool split_illconditioned_faces(Mesh &mesh)
{
	vector<Edge*> bad_edges;
	for (int i = 0; i < mesh.faces.size(); i++) {
		Face *f0 = mesh.faces[i];
		int cc = conformalCount(f0);
		if (cc == 1) {
			Node *n0 = single_eol_from_face(f0);
			Edge *e0 = get_opp_edge(f0, n0);
			if (face_altitude(e0, f0) < thresh / 2) {
				cout << face_altitude(e0, f0) << endl;
				bad_edges.push_back(e0);
			}
		}
		else if (cc == 2) {
			Edge *e0 = single_conformal_edge_from_face(f0);
			if (face_altitude(e0, f0) < (thresh / 2)) {
				cout << face_altitude(e0, f0) << endl;
				bad_edges.push_back(e0);
			}
		}
		else if(cc == 3) {
			Edge *e0 = single_nonconformal_edge_from_face(f0);
			if (face_altitude(e0, f0) < thresh / 2) {
				cout << face_altitude(e0, f0) << endl;
				bad_edges.push_back(e0);
			}
		}
	}
	int exclude = 0;
	for (size_t e = 0; e < bad_edges.size(); e++) {
		Edge *edge = bad_edges[e];
		if (!edge) continue;
		Node *node0 = edge->n[0], *node1 = edge->n[1];
		RemeshOp op = split_edgeForced(edge, 0.5, thresh);
		if (op.empty()) {
			exclude++;
			continue;
		}
		for (size_t v = 0; v < op.added_verts.size(); v++) {
			Vert *vertnew = op.added_verts[v];
			Vert *v0 = adjacent_vert(node0, vertnew),
				*v1 = adjacent_vert(node1, vertnew);
			vertnew->sizing = 0.5 * (v0->sizing + v1->sizing);
		}
		// We don't want to try and split these triangles only to make them worse for the next pass
		//bool make_worse = true;
		//for (size_t n = 0; n < op.added_nodes.size(); n++) {
		//	for (int adje = 0; adje < op.added_nodes[n]->adje.size(); adje++) {
		//		if (edge_length(op.added_nodes[n]->adje[adje]) < thresh) make_worse = false;
		//	}
		//}
		//if (make_worse) op.cancel();
		op.set_null(bad_edges);
		op.done();
	}
	return bad_edges.size() == 0 + exclude;
}

void flip_edges(MeshSubset* subset, vector<Face*>& active_faces,
	vector<Edge*>* update_edges, vector<Face*>* update_faces);


// TODO:: Does conformal stalling occur with EOL?
void cleanup(Mesh& mesh)
{
	vector<Face*> active_faces = mesh.faces;
	flip_edges(0, active_faces, 0, 0);
	markPreserve(mesh);
	bool allclear = false;
	while (!allclear) {
		// Iterate until all the bad edges are accounted for
		// If a bad edge exists, but is unsafe to collapse, try in the next iteration where it may become safe, or another operation may take care of it
		while(!allclear) {
			allclear = true;
			while (collapse_nonconformal(mesh, allclear));
			markPreserve(mesh); // Can be removed if the conformal pass doesn't loop over edges and check for preserve status
			while (collapse_conformal(mesh, allclear));
		}
		allclear = split_illconditioned_faces(mesh);
		//allclear = true;
	}
	markPreserve(mesh); // Probably doesn't need to be called so much, but wan't to be safe
}

// TODO:: I think there are problems when a box corner reaches the cloth border but the two box edges still move through the cloth

void preprocess(Mesh& mesh, const vector<shared_ptr<btc::Collision> > cls)
{
	markWasEOL(mesh);
	addGeometry(mesh, cls);
	revertWasEOL(mesh);

	markPreserve(mesh);

	cleanup(mesh);

	compute_ws_data(mesh);
}

void preprocessClean(Mesh& mesh)
{
	cleanup(mesh);
	compute_ws_data(mesh);
}

bool a;
void preprocessPart(Mesh& mesh, const vector<shared_ptr<btc::Collision> > cls, int &part)
{
	if (part == 1) {
		//for (int i = 0; i < mesh.nodes.size(); i++) {
		//	mesh.nodes[i]->EoL = false;
		//	mesh.nodes[i]->cornerID = -1;
		//	mesh.nodes[i]->cdEdges.clear();
		//}
		//for (int i = 0; i < mesh.edges.size(); i++) {
		//	mesh.edges[i]->preserve = false;
		//}
		a = true;
		markWasEOL(mesh);
		addGeometry(mesh, cls);
		revertWasEOL(mesh);
		cout << "Add Geometry" << endl;
	}
	else if (part == 2) {
		markPreserve(mesh);
		cout << "Mark preserve" << endl;
		//part = 6;
	}
	else if (part == 3) {
		vector<Face*> active_faces = mesh.faces;
		flip_edges(0, active_faces, 0, 0);
		markPreserve(mesh);
		cout << "Flipped edges" << endl;
	}
	else if (part == 4) {
		a = true;
		while (collapse_nonconformal(mesh,a));
		markPreserve(mesh);
		cout << "Collapse nonfornformal" << endl;
	}
	else if (part == 5) {
		while (collapse_conformal(mesh,a));
		if (!a) part = 3;
		cout << "Collapse conformal" << endl;
	}
	else if (part == 6) {
		bool allclear = split_illconditioned_faces(mesh);
		if (!allclear) {
			part = 3;
			cout << "Split ill-conditioned, not good" << endl;
		}
		else {
			markPreserve(mesh);
			cout << "Split ill-conditioned, all good" << endl;
		}
	}
	else if (part == 7) {
		compute_ws_data(mesh);
		cout << "Compute ws data" << endl;
	}
}