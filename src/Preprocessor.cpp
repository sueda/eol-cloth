#include "Preprocessor.h"

#include "remeshExtension.h"
#include "matlabOutputs.h"

#include "external\ArcSim\geometry.hpp"
#include "external\ArcSim\subset.hpp"

#include <stdlib.h>

using namespace std;
//using namespace Eigen;

double upp = 0.0;
int wh = 0;

double thresh = 0.01;

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

void addGeometry(Mesh& mesh, vector<shared_ptr<btc::Collision> > cls)
{
	for (int i = 0; i < cls.size(); i++) {
		//if (i > 8) break;
		if (cls[i]->count1 == 1 && cls[i]->count2 == 3) {
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
				
			// Faces CAN be added and deleted throughout this loop so we can't just use the CD returned tri
			Face *f0 = get_enclosing_face(mesh, Vec2(xX, yX));
			Vec3 bary = get_barycentric_coords(Vec2(xX, yX), f0);

			RemeshOp op = split_face(f0, bary);
			for (size_t v = 0; v < op.added_verts.size(); v++) {
				Vert *vertnew = op.added_verts[v];
				vertnew->sizing = (v0->sizing + v1->sizing + v2->sizing) / 3.0;
			}
			op.done();

			Node *n = mesh.nodes.back();
			n->EoL = true;
			n->preserve = true;
			n->cornerID = cls[i]->verts1(0);
			n->cdEdges = cls[i]->edge1;
		}
		if (cls[i]->count1 == 2 && cls[i]->count2 == 2) {
			// TODO:: Is this enough of a check?
			if (mesh.nodes[cls[i]->verts2(0)]->EoL ||
				mesh.nodes[cls[i]->verts2(1)]->EoL) continue;

			// The CD edges don't map to the changing mesh edges so we need to find it
			// The verts won't change since we only ever add verts
			Edge *e0 = get_edge(mesh.nodes[cls[i]->verts2(0)], mesh.nodes[cls[i]->verts2(1)]);
			double d;
			if (e0 == NULL) {
				Vert *v0 = mesh.verts[cls[i]->verts2(0)],
					*v1 = mesh.verts[cls[i]->verts2(1)];

				double xX = cls[i]->weights2(0) * v0->u[0] +
					cls[i]->weights2(1) * v1->u[0];
				double yX = cls[i]->weights2(0) * v0->u[1] +
					cls[i]->weights2(1) * v1->u[1];

				// Faces CAN be added and deleted throughout this loop so we can't just use the CD returned tri
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
					d = bary[0];
				}
				else if (e0->n[0] == f0->v[1]->node) {
					d = bary[1];
				}
				else {
					d = bary[2];
				}
				cout << endl;
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
			RemeshOp op = split_edgeForced(e0, d);
			for (size_t v = 0; v < op.added_verts.size(); v++) {
				Vert *vertnew = op.added_verts[v];
				Vert *v0 = adjacent_vert(node0, vertnew),
					*v1 = adjacent_vert(node1, vertnew);
				vertnew->sizing = 0.5 * (v0->sizing + v1->sizing);
			}
			op.done();

			Node *n = mesh.nodes.back();
			n->EoL = true;
			n->cdEdges = cls[i]->edge1;
		}
	}

	//mesh2m(mesh, "mesh.m", true);
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

bool collapse_conformal(Mesh &mesh)
{
	for (int i = 0; i < mesh.edges.size(); i++) {
		Edge *e = mesh.edges[i];
		if (e->preserve) {
			if (edge_length(e) < thresh) {
				RemeshOp op;
				Node *n0 = e->n[0],
					*n1 = e->n[1];
				if (is_seam_or_boundary(n1) || n1->cornerID >= 0) {
					op = collapse_edgeForced(e, 0);
					if (op.empty()) continue;
					op.done();
					return true;
				}
				else if (is_seam_or_boundary(n0) || n0->cornerID >= 0) {
					op = collapse_edgeForced(e, 1);
					if (op.empty()) continue;
					op.done();
					return true;
				}
				else {
					if (n0->verts[0]->adjf.size() <= n1->verts[0]->adjf.size()) {
						op = collapse_edgeForced(e, 1);
						if (op.empty()) {
							op = collapse_edgeForced(e, 0);
							if (op.empty()) {
								continue;
							}
						}
					}
					else {
						op = collapse_edgeForced(e, 0);
						if (op.empty()) {
							op = collapse_edgeForced(e, 1);
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
	return false;
}

bool collapse_nonconformal(Mesh &mesh)
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
						// Don't deal with edges between boundary and inside
						if (!(
							(is_seam_or_boundary(n0) && is_seam_or_boundary(n1)) ||
							(!is_seam_or_boundary(n0) && !is_seam_or_boundary(n1))
							)) continue;
						RemeshOp op;
						if (n0->EoL) {
							op = collapse_edgeForced(e0, 1);
							if (op.empty()) continue;
							op.done();
							return true;
						}
						else if (n1->EoL) {
							op = collapse_edgeForced(e0, 0);
							if (op.empty()) continue;
							op.done();
							return true;
						}
						else {
							if (!n0->preserve) {
								op = collapse_edgeForced(e0, 0);
								if (op.empty()) {
									op = collapse_edgeForced(e0, 1);
									if (op.empty()) {
										continue;
									}
								}
							}
							else if (!n1->preserve) {
								op = collapse_edgeForced(e0, 1);
								if (op.empty()) {
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

bool collapse_close(Mesh &mesh)
{
	for (int i = 0; i < mesh.nodes.size(); i++) {
		Node *n0 = mesh.nodes[i];
		if (n0->EoL) {
			for (int e = 0; e < n0->adje.size(); e++) {
				Edge *e0 = n0->adje[e];
				Node *n1 = other_node(e0, n0);
				if (!n1->EoL && !is_seam_or_boundary(n1) && edge_length(e0) < thresh) {
					RemeshOp op;
					if (n0 == e0->n[0]) op = collapse_edgeForced(e0, 1);
					else op = collapse_edgeForced(e0, 0);
					if (op.empty()) continue;
					op.done();
					return true;
				}
			}
		}
	}
	return false;
}

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
	for (size_t e = 0; e < bad_edges.size(); e++) {
		Edge *edge = bad_edges[e];
		if (!edge) continue;
		Node *node0 = edge->n[0], *node1 = edge->n[1];
		RemeshOp op = split_edgeForced(edge, 0.5);
		for (size_t v = 0; v < op.added_verts.size(); v++) {
			Vert *vertnew = op.added_verts[v];
			Vert *v0 = adjacent_vert(node0, vertnew),
				*v1 = adjacent_vert(node1, vertnew);
			vertnew->sizing = 0.5 * (v0->sizing + v1->sizing);
		}
		op.set_null(bad_edges);
		op.done();
	}
	return bad_edges.size() == 0;
}

void flip_edges(MeshSubset* subset, vector<Face*>& active_faces,
	vector<Edge*>* update_edges, vector<Face*>* update_faces);

void cleanup(Mesh& mesh)
{
	vector<Face*> active_faces = mesh.faces;
	flip_edges(0, active_faces, 0, 0);
	bool allclear = false;
	while (!allclear) {
		while (collapse_nonconformal(mesh));
		while (collapse_conformal(mesh));
		allclear = split_illconditioned_faces(mesh);
		allclear = true;
	}
}

// TODO::
// Boundary condition
// Don't collapse corners 
// ArcSim sizing on new nodes?

void preprocess(Mesh& mesh, vector<shared_ptr<btc::Collision> > cls)
{
	//double r0 = (double) rand() / RAND_MAX;
	//double r1 = (double) rand() / RAND_MAX;
	////double r0 = 0.6 + upp;
	////double r1 = 0.65 + upp;
	////upp += 0.1;
	//Face *f = get_enclosing_face(mesh, Vec2(r0, r1));
	//Vec3 bary = get_barycentric_coords(Vec2(r0, r1), f);
	////RemeshOp op = split_face(f, bary);
	//////op.update(mesh);
	////op.done();
	////RemeshOp op2 = split_edge(mesh.edges[wh], .3);
	////op2.done();
	////wh++;
	//
	//vector<shared_ptr<btc::Collision> > test;
	//auto c = make_shared<btc::Collision>();
	//c->count1 = 1;
	//c->count2 = 3;
	//c->verts1 << wh, -1, -1;
	//c->verts2 << f->v[0]->index, f->v[1]->index, f->v[2]->index;
	//c->weights2 << bary[0], bary[1], bary[2];
	//c->tri2 = f->index;
	//c->edge1.push_back(0);
	//test.push_back(c);
	//wh++;

	//int r2 = rand() % (mesh.edges.size()-1);
	//double r3 = (double)rand() / RAND_MAX;
	//Edge *e = mesh.edges[r2];
	//auto c1 = make_shared<btc::Collision>();
	//c1->count1 = 2;
	//c1->count2 = 2;
	//c1->verts2 << e->n[0]->index, e->n[1]->index, -1;
	//c1->weights2 << 1.0 - r3, r3, 0.0;
	//c1->edge1.push_back(0);
	//test.push_back(c1);

	for (int i = 0; i < mesh.nodes.size(); i++) {
		mesh.nodes[i]->EoL = false;
		mesh.nodes[i]->cornerID = -1;
		mesh.nodes[i]->cdEdges.clear();
	}
	for (int i = 0; i < mesh.edges.size(); i++) {
		mesh.edges[i]->preserve = false;
	}

	addGeometry(mesh, cls);

	//while (collapse_close(mesh));

	markPreserve(mesh);

	//cleanup(mesh);

	compute_ws_data(mesh);
}

void preprocessPart(Mesh& mesh, vector<shared_ptr<btc::Collision> > cls, int &part)
{
	if (part == 1) {
		addGeometry(mesh, cls);
		cout << "Add Geometry" << endl;
	}
	//else if (part == 2) {
	//	//while (collapse_close(mesh));
	//	cout << "Collapse Close" << endl;
	//}
	else if (part == 2) {
		markPreserve(mesh);
		cout << "Mark preserve" << endl;
		part = 6;
	}
	else if (part == 3) {
		vector<Face*> active_faces = mesh.faces;
		flip_edges(0, active_faces, 0, 0);
		cout << "Flipped edges" << endl;
	}
	else if (part == 4) {
		while (collapse_nonconformal(mesh));
		cout << "Collapse nonfornformal" << endl;
	}
	else if (part == 5) {
		while (collapse_conformal(mesh));
		cout << "Collapse conformal" << endl;
	}
	else if (part == 6) {
		bool allclear = split_illconditioned_faces(mesh);
		if (!allclear) {
			part = 3;
			cout << "Split ill-conditioned, not good" << endl;
		}
		else {
			cout << "Split ill-conditioned, all good" << endl;
		}
	}
	else if (part == 7) {
		compute_ws_data(mesh);
		cout << "Compute ws data" << endl;
	}
}