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

#ifndef MESH_HPP
#define MESH_HPP

#include "transformation.hpp"
#include "vectors.hpp"
#include <utility>
#include <vector>
#include <memory>

struct Serialize;

// material space (not fused at seams)
struct Vert;
struct Face;
// world space (fused)
struct Node;
struct Edge;
struct Material;
struct ReferenceShape;
struct Cloth;
struct Mesh;

struct Plane {
	Plane() {}
	Plane(const Vec3& x0, const Vec3& n) : x0(x0), n(n) {}
	Vec3 x0, n;
};

extern int uuid_src;

struct Vert {
	Vec3 u; // material space
	// NICK
	Vec3 v; // material velocity
	Node *node; // world space
				// topological data
	std::vector<Face*> adjf; // adjacent faces
	int index; // position in mesh.verts
			   // derived material-space data that only changes with remeshing
			   // remeshing data
	Mat3x3 sizing;
	// constructors
	Vert() : node(0), index(-1) {}
	explicit Vert(const Vec3 &u, const Vec3 &v) :
		u(u), v(v) {}

	//void serializer(Serialize& s);
};

struct Node {
	enum NodeFlags {
		FlagNone = 0, FlagActive = 1, FlagMayBreak = 2,
		FlagResolveUni = 4, FlagResolveMax = 8
	};

	int uuid;

	Mesh* mesh;

	bool EoL;
	int EoL_index;
	int cornerID;
	std::vector<int> cdEdges;

	int label;
	int flag;
	std::vector<Vert*> verts;
	Vec3 y; // plastic embedding
	Vec3 x, x0, v; // position, old (collision-free) position, velocity
	bool preserve; // don't remove this node
				   // topological data
	int index; // position in mesh.nodes
	std::vector<Edge*> adje; // adjacent edges
							 // derived world-space data that changes every frame
	Vec3 n; // local normal, approximate
			// derived material-space data that only changes with remeshing
	double a, m; // area, mass
	Mat3x3 curvature; // filtered curvature for bending fracture
					  // pop filter data
	Vec3 acceleration;
	Node() : uuid(uuid_src++), label(0), flag(0), preserve(false), index(-1), a(0), m(0) {}
	explicit Node(const Vec3 &y, const Vec3 &x, const Vec3 &v, int label, int flag,
		bool preserve) :
		EoL(false), EoL_index(-1), cornerID(-1), // EoL specifics
		uuid(uuid_src++), mesh(0), label(label), flag(flag), y(y), x(x), x0(x), v(v), preserve(preserve),
		curvature(0) {}

	inline bool active() const { return flag & FlagActive; }

	//void serializer(Serialize& s);
};

struct Edge {
	Node *n[2]; // nodes
	int preserve;
	// topological data
	Face *adjf[2]; // adjacent faces
	int index; // position in mesh.edges
			   // plasticity data
	double theta_ideal, damage; // rest dihedral angle, damage parameter
								// constructors
	Edge() : index(-1), theta_ideal(0), damage(0) { n[0] = n[1] = 0; adjf[0] = adjf[1] = 0; }
	explicit Edge(Node *node0, Node *node1, double theta_ideal, int preserve) :
		preserve(preserve), theta_ideal(theta_ideal), damage(0) {
		n[0] = node0;
		n[1] = node1;
	}

	//void serializer(Serialize& s);
};

struct Face {
	Vert* v[3]; // verts
	Material* material;
	int flag;
	// topological data
	Edge *adje[3]; // adjacent edges
	int index; // position in mesh.faces
			   // derived world-space data that changes every frame
	Vec3 n; // local normal, exact
			// derived material-space data that only changes with remeshing
	double a, m; // area, mass
	Mat3x3 Dm, invDm; // finite element matrix
					  // plasticity data
	Mat3x3 Sp_bend; // plastic bending strain
	Mat3x3 Sp_str; // plastic stretching
	Mat3x3 sigma;
	double damage; // accumulated norm of S_plastic/S_yield
				   // constructors
	Face() : material(0), flag(0), index(-1), a(0), m(0), damage(0) {
		for (int i = 0; i<3; i++) { v[i] = 0; adje[i] = 0; }
	}
	explicit Face(Vert *vert0, Vert *vert1, Vert *vert2, const Mat3x3& ps,
		const Mat3x3& pb, Material* mat, double damage) :
		material(mat), flag(0), a(0), m(0), Sp_bend(pb), Sp_str(ps), sigma(0), damage(damage) {
		v[0] = vert0;
		v[1] = vert1;
		v[2] = vert2;
	}

	//void serializer(Serialize& s);
};

struct Mesh {
	ReferenceShape *ref;
	std::shared_ptr<Cloth> parent;
	//CollisionProxy* proxy;

	int EoL_Count;

	std::vector<Vert*> verts;
	std::vector<Node*> nodes;
	std::vector<Edge*> edges;
	std::vector<Face*> faces;
	// These do *not* assume ownership, so no deletion on removal
	void add(Vert *vert);
	void add(Node *node);
	void add(Edge *edge);
	void add(Face *face);
	void remove(Vert *vert);
	void remove(Node *node);
	void remove(Edge *edge);
	void remove(Face *face);

	Mesh() : ref(0), parent(0), EoL_Count(0) {};

	//void serializer(Serialize& s);
};

template <typename Prim> inline const std::vector<Prim*> &get(const Mesh &mesh);
template <typename Prim> inline int count_elements(const std::vector<Mesh*>& meshes);

void connect(Vert *vert, Node *node); // assign vertex to node

bool check_that_pointers_are_sane(const Mesh &mesh);
bool check_that_contents_are_sane(const Mesh &mesh);

void compute_ms_data(Mesh &mesh); // call after mesh topology changes
void compute_ws_data(Mesh &mesh); // call after vert positions change
void compute_ms_data(std::vector<Face*>& face);
void compute_ms_data(std::vector<Node*>& node);
void compute_ws_data(std::vector<Face*>& face);
void compute_ws_data(std::vector<Node*>& node);
void compute_ms_data(Face* face);
void compute_ms_data(Node* node);

// NICK
//double calc_edge_weight(Node* n);

inline Edge* get_opp_edge(const Face* face, const Node* opp);
inline Vert* get_vert(const Face* face, const Node* node);
inline Edge *get_edge(const Node *node0, const Node *node1);
inline Vert *edge_vert(const Edge *edge, int side, int i);
inline Vert *edge_opp_vert(const Edge *edge, int side);
inline Node *other_node(const Edge* edge, const Node* node0);
inline Face *adj_face(const Face* face0, int num);

inline Edge *next_edge_ccw(const Edge* edge, Node* center);
inline Edge *next_edge_cw(const Edge* edge, Node* center);
inline Face *next_face_ccw(const Edge* edge, Node* center);
inline Face *next_face_cw(const Edge* edge, Node* center);

void set_indices(Mesh &mesh);
void set_indices(std::vector<Mesh*> &meshes);
void mark_nodes_to_preserve(Mesh &mesh);

inline Vec3 derivative(double a0, double a1, double a2, double az, const Face *face);
inline Mat3x3 derivative(const Vec3& w0, const Vec3& w1,
	const Vec3& w2, const Vec3& dz, const Face *face);

void apply_transformation_onto(const Mesh& start_state, Mesh& onto,
	const Transformation& tr);
void apply_transformation(Mesh& mesh, const Transformation& tr);

void update_x0(Mesh &mesh);

Mesh deep_copy(Mesh &mesh);
void delete_mesh(Mesh &mesh);

// ADDED BY NICK
void reindex_nodes(std::vector<Node*>& nodes);

void activate_nodes(std::vector<Node*>& nodes);
void deactivate_nodes(std::vector<Node*>& nodes);

//
// IMPLEMENTATION OF INLINE FUNCTIONS
//

inline Vec3 derivative(double a0, double a1, double a2, double az, const Face *face) {
	return face->invDm.t() * Vec3(a1 - a0, a2 - a0, az);
}

inline Mat3x3 derivative(const Vec3& w0, const Vec3& w1,
	const Vec3& w2, const Vec3& dz, const Face *face) {
	return Mat3x3(w1 - w0, w2 - w0, dz) * face->invDm;
}

inline Vert* get_vert(const Face* face, const Node* node) {
	if (face->v[0]->node == node) return face->v[0];
	return face->v[1]->node == node ? face->v[1] : face->v[2];
}

inline Node *other_node(const Edge* edge, const Node* node0) {
	return edge->n[0] == node0 ? edge->n[1] : edge->n[0];
}

inline Face *adj_face(const Face* face0, int num) {
	Edge* e = face0->adje[num];
	return e->adjf[0] == face0 ? e->adjf[1] : e->adjf[0];
}

inline Edge *next_edge_ccw(const Edge* edge, Node* center) {
	Face *f = next_face_ccw(edge, center);
	Node *os = other_node(edge, center);
	for (int i = 0; i<3; i++) {
		Edge *e = f->adje[i];
		if ((e->n[0] == center || e->n[1] == center) && e->n[0] != os && e->n[1] != os)
			return e;
	}
	return NULL;
}

inline Edge *next_edge_cw(const Edge* edge, Node* center) {
	Face *f = next_face_cw(edge, center);
	Node *os = other_node(edge, center);
	for (int i = 0; i<3; i++) {
		Edge *e = f->adje[i];
		if ((e->n[0] == center || e->n[1] == center) && e->n[0] != os && e->n[1] != os)
			return e;
	}
	return NULL;
}

inline Face *next_face_ccw(const Edge* edge, Node* center) {
	return (edge->n[0] == center) ? edge->adjf[0] : edge->adjf[1];
}

inline Face *next_face_cw(const Edge* edge, Node* center) {
	return (edge->n[1] == center) ? edge->adjf[0] : edge->adjf[1];
}

// NICK
inline Edge* get_opp_edge(const Face* face, const Node* opp) {
	for (int e = 0; e < 3; e++) {
		Edge *edge = face->adje[e];
		if (edge->n[0] == opp || edge->n[1] == opp) continue;
		return edge;
	}
	return NULL;
}

inline Edge *get_edge(const Node *n0, const Node *n1) {
	for (int e = 0; e < (int)n0->adje.size(); e++) {
		Edge *edge = n0->adje[e];
		if (edge->n[0] == n1 || edge->n[1] == n1)
			return edge;
	}
	return NULL;
}

inline Vert *edge_vert(const Edge *edge, int side, int i) {
	Face *face = (Face*)edge->adjf[side];
	if (!face)
		return NULL;
	for (int j = 0; j < 3; j++)
		if (face->v[j]->node == edge->n[i])
			return face->v[j];
	return NULL;
}

inline Vert *edge_opp_vert(const Edge *edge, int side) {
	Face *face = (Face*)edge->adjf[side];
	if (!face)
		return NULL;
	for (int j = 0; j < 3; j++)
		if (face->v[j]->node == edge->n[side])
			return face->v[j>0 ? j - 1 : j + 2];
	return NULL;
}

template <> inline const std::vector<Vert*> &get(const Mesh &mesh) { return mesh.verts; }
template <> inline const std::vector<Node*> &get(const Mesh &mesh) { return mesh.nodes; }
template <> inline const std::vector<Edge*> &get(const Mesh &mesh) { return mesh.edges; }
template <> inline const std::vector<Face*> &get(const Mesh &mesh) { return mesh.faces; }

template <typename Prim> inline int count_elements(const std::vector<Mesh*>& meshes) {
	int num = 0;
	for (size_t i = 0; i< meshes.size(); i++)
		num += get<Prim>(*meshes[i]).size();
	return num;
}

#endif
