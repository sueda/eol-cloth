#include "remeshExtension.h"

#include "external\ArcSim\util.hpp"
#include "external\ArcSim\blockvectors.hpp"
#include "external\ArcSim\geometry.hpp"

RemeshOp split_edgeForced(Edge* edge, double d, double thresh) {
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



		op.added_faces.push_back(nf0);
		op.added_faces.push_back(nf1);
	}

	bool make_worse = true;
	if (thresh > 0) {
		for (size_t n = 0; n < op.added_nodes.size(); n++) {
			for (int adje = 0; adje < op.added_nodes[n]->adje.size(); adje++) {
				if (edge_length(op.added_nodes[n]->adje[adje]) < thresh) make_worse = false;
			}
		}
	}
	else {
		make_worse = false;
	}

	if (make_worse) {
		op.cancel();
	}
	else {
		op.apply(mesh);
		node->y = (1 - d)*node0->y + d*node1->y;
	}
	return op;
}

RemeshOp collapse_edgeForced(Edge* edge, int i) {
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
		// Preserve in weird situations
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
				if (area(new_face) == 0) {
					op.cancel();
					return RemeshOp();
				}
			}
		}
	}
	op.apply(mesh);
	return op;
}

RemeshOp split_face(Face* face, Vec3 b) {
	Mesh& mesh = *face->v[0]->node->mesh;
	RemeshOp op;
	Node *node0 = face->v[0]->node,
		*node1 = face->v[1]->node,
		*node2 = face->v[2]->node,
		*node = new Node(b[0] * node0->y + b[1] * node1->y + b[2] * node2->y,
			b[0] * node0->x + b[1] * node1->x + b[2] * node2->x,
			b[0] * node0->v + b[1] * node1->v + b[2] * node2->v,
			0,
			node0->flag & node1->flag & node2->flag,
			false);
	node->acceleration = b[0] * node0->acceleration + b[1] * node1->acceleration + b[2] * node2->acceleration;
	op.added_nodes.push_back(node);
	op.added_edges.push_back(new Edge(node0, node, 0.0,
		false));
	op.added_edges.push_back(new Edge(node1, node, 0.0,
		false));
	op.added_edges.push_back(new Edge(node2, node, 0.0,
		false));
	Vert *v0 = face->v[0],
		*v1 = face->v[1],
		*v2 = face->v[2];
	Vert *v = new Vert(b[0] * v0->u + b[1] * v1->u + b[2] * v2->u,
		b[0] * v0->v + b[1] * v1->v + b[2] * v2->v);
	v->sizing = b[0] * v0->sizing + b[1] * v1->sizing + b[2] * v2->sizing;
	connect(v, node);
	op.added_verts.push_back(v);
	op.removed_faces.push_back(face);
	Face* nf0 = new Face(v0, v1, v, face->Sp_str, face->Sp_bend, face->material, face->damage);
	Face* nf1 = new Face(v1, v2, v, face->Sp_str, face->Sp_bend, face->material, face->damage);
	Face* nf2 = new Face(v2, v0, v, face->Sp_str, face->Sp_bend, face->material, face->damage);

	op.added_faces.push_back(nf0);
	op.added_faces.push_back(nf1);
	op.added_faces.push_back(nf2);

	op.apply(mesh);
	return op;
}