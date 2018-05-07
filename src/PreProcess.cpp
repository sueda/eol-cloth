#include "PreProcess.h"
#include "Cloth.h"
#include "Box.h"
#include "boxTriCollision.h"
#include "MatlabDebug.h"
#include "ChronoTimer.h"
#include "CVM.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <stdlib.h>

// Fade2D
#include "Fade_2D.h"

// ArcSim
#include "mesh.hpp"
#include "vectors.hpp"
#include "util.hpp"
#include "io.hpp"
#include "geometry.hpp"
#include "remesh.hpp"

using namespace std;
using namespace Eigen;

//struct CVM
//{
//	shared_ptr<btc::Collision> coll;
//	//int* vert_index;
//	Vert* vert;
//	CVM(shared_ptr<btc::Collision> c, Vert* v) : coll(c), vert(v) {
//
//	}
//	CVM(Vert* v) : vert(v) {}
//};

void mini_reindex_nodes(vector<Node*>& nodes) {
	for (size_t i = 0; i<nodes.size(); i++) {
		nodes[i]->index = i;
		nodes[i]->verts[0]->index = i;
	}
}

// TODO:: Corner already exists case
double compare_edge_weights(shared_ptr<CVM> a, shared_ptr<CVM> b) {
	return (a->vert->egde_weight[a->which_weight] > b->vert->egde_weight[b->which_weight]);
	//if (a->coll != NULL && b->coll != NULL) {
	//	if (a->coll->count1 == 3 && b->coll->count1 == 3) return (a->coll->weights2(0) > b->coll->weights2(0));
	//	if (a->coll->count1 == 3) return (a->coll->weights2(0) > b->coll->weights1(0));
	//	if (b->coll->count1 == 3) return (a->coll->weights1(0) > b->coll->weights2(0));
	//	return (a->coll->weights1(0) > b->coll->weights1(0));
	//}
	//else if (a->coll != NULL && b->coll == NULL) {
	//	return (a->coll->weights1(0) > b->vert->egde_weight[b->which_weight]);
	//}
	//else if (a->coll == NULL && b->coll != NULL) {
	//	return (a->vert->egde_weight[a->which_weight] > b->coll->weights1(0));
	//}
	//else if (a->coll == NULL && b->coll == NULL) {
	//	return (a->vert->egde_weight[a->which_weight] > b->vert->egde_weight[b->which_weight]);
	//}
}

//double calc_edge_length(const Edge* edge)
//{
//	return sqrt(pow(edge->n[0]->verts[0]->u[0] - edge->n[1]->verts[0]->u[0], 2) + pow(edge->n[0]->verts[0]->u[1] - edge->n[1]->verts[0]->u[1], 2));
//}

bool to_flip_corner(const Vert* a, const Vert* b, const Vert* c)
{
	return (sqrt(pow(a->u[0] - b->u[0], 2) + pow(a->u[1] - b->u[1], 2)) > sqrt(pow(a->u[0] - c->u[0], 2) + pow(a->u[1] - c->u[1], 2)));
}

template <typename T>
void moveItemToBack(std::vector<T>& v, size_t itemIndex)
{
	auto it = v.begin() + itemIndex;
	std::rotate(it, it + 1, v.end());
}

double find_split_weight(const Edge* pedge, const Vert* v1, const Vert* v2)
{
	//Vec3 r = pedge->n[1]->verts[0]->u - pedge->n[0]->verts[0]->u;
	//Vec3 s = v2->u - v1->u;
	//Vec3 rxs = cross(r, s);
	//double t = cross((v1->u - pedge->n[0]->verts[0]->u), r) / rxs;

	double s1x = pedge->n[1]->verts[0]->u[0] - pedge->n[0]->verts[0]->u[0];
	double s1y = pedge->n[1]->verts[0]->u[1] - pedge->n[0]->verts[0]->u[1];
	double s2x = v2->u[0] - v1->u[0];
	double s2y = v2->u[1] - v1->u[1];

	double adf = s2x * (pedge->n[0]->verts[0]->u[1] - v1->u[1]);
	double ahkj = s2y * (pedge->n[0]->verts[0]->u[0] - v1->u[0]);
	double see = (-s2x * s1y + s1x * s2y);

	double t = (s2x * (pedge->n[0]->verts[0]->u[1] - v1->u[1]) - s2y * (pedge->n[0]->verts[0]->u[0] - v1->u[0])) / (-s2x * s1y + s1x * s2y);

	Vec3 intp = Vec3((pedge->n[0]->verts[0]->u[0] + (t * s1x)), (pedge->n[0]->verts[0]->u[1] + (t * s1y)), 0.0);

	if ((intp[0] > pedge->n[0]->verts[0]->u[0] && intp[0] > pedge->n[1]->verts[0]->u[0]) || (intp[0] < pedge->n[0]->verts[0]->u[0] && intp[0] < pedge->n[1]->verts[0]->u[0])) {
		return -1.0;
	}

	double new_seg_length = sqrt(pow(intp[0] - pedge->n[0]->verts[0]->u[0], 2) + pow(intp[1] - pedge->n[0]->verts[0]->u[1], 2));
	double full_seg_length = sqrt(pow(pedge->n[1]->verts[0]->u[0] - pedge->n[0]->verts[0]->u[0], 2) + pow(pedge->n[1]->verts[0]->u[1] - pedge->n[0]->verts[0]->u[1], 2));

	return new_seg_length / full_seg_length;
}

double calc_nonexistant_edge_length(const Vec2& a)
{
	//return sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2));
	return sqrt(a[0] * a[0] + a[1] * a[1]);
}

bool edge_not_within_safety_margin(const Vert* v1, const Vert* v2, const double& margin, const Vector2d& bounds)
{
	if (v1->u[0] < margin ) { //&& v2->u[0] < margin) {
		Vec2 vc1 = Vec2(v2->u[0] - v1->u[0], v2->u[1] - v1->u[1]);
		Vec2 vc2 = Vec2(0.0 - 0.0, bounds[1] - 0.0);
		double angle = acos(dot(vc1, vc2) / (calc_nonexistant_edge_length(vc1) * calc_nonexistant_edge_length(vc2)));
		if (angle < M_PI / 5 || angle >(4 * M_PI / 5)) return true;
	}
	if (v1->u[1] < margin) { // && v2->u[1] < margin) {
		Vec2 vc1 = Vec2(v2->u[0] - v1->u[0], v2->u[1] - v1->u[1]);
		Vec2 vc2 = Vec2(bounds[0] - 0.0, 0.0 - 0.0);
		double angle = acos(dot(vc1, vc2) / (calc_nonexistant_edge_length(vc1) * calc_nonexistant_edge_length(vc2)));
		if (angle < M_PI / 5 || angle >(4 * M_PI / 5)) return true;
	}
	if (v1->u[0] > (bounds[0] - margin) ) { //&& v2->u[0] > (1.0 - margin)) {
		Vec2 vc1 = Vec2(v2->u[0] - v1->u[0], v2->u[1] - v1->u[1]);
		Vec2 vc2 = Vec2(0.0 - 0.0, bounds[1] - 0.0);
		double angle = acos(dot(vc1, vc2) / (calc_nonexistant_edge_length(vc1) * calc_nonexistant_edge_length(vc2)));
		if (angle < M_PI / 5 || angle >(4 * M_PI / 5)) return true;
	}
	if (v1->u[1] > (bounds[1] - margin) ) { //&& v2->u[1] > (1.0 - margin)) {
		Vec2 vc1 = Vec2(v2->u[0] - v1->u[0], v2->u[1] - v1->u[1]);
		Vec2 vc2 = Vec2(bounds[1] - 0.0, 0.0 - 0.0);
		double angle = acos(dot(vc1, vc2) / (calc_nonexistant_edge_length(vc1) * calc_nonexistant_edge_length(vc2)));
		if (angle < M_PI / 5 || angle >(4 * M_PI / 5)) return true;
	}

	return false;

	//return v1->u[0] < margin ||
	//	v1->u[0] > (1.0 - margin) ||
	//	v1->u[1] < margin ||
	//	v1->u[1] > (1.0 - margin);
}

bool corner_not_within_safety_margin(const Vec3 v, const double& margin, const Vector2d& bounds)
{
	return (v[0] < margin ||
		v[0] > (bounds[0] - margin) ||
		v[1] < margin ||
		v[1] > (bounds[1] - margin));
	//return (v[0] < margin && v[1] < margin) ||
	//	(v[0] < margin && v[1] > (1.0 - margin)) ||
	//	(v[0] > (1.0 - margin) && v[1] < margin) ||
	//	(v[0] > (1.0 - margin) && v[1] > (1.0 - margin));
}

bool quick_edge_safety_margin_check(const Vec3 v, const double& margin, const Vector2d& bounds)
{
	return (v[0] < margin && v[1] < margin) ||
		(v[0] < margin && v[1] > (bounds[1] - margin)) ||
		(v[0] > (bounds[0] - margin) && v[1] < margin) ||
		(v[0] > (bounds[0] - margin) && v[1] > (bounds[1] - margin));
}

void export2D(const Mesh& mesh, string file_name) {
	vector<GEOM_FADE2D::Point2> vInputPoints;
	for (int i = 0; i < mesh.verts.size(); i++) {
		vInputPoints.push_back(GEOM_FADE2D::Point2(mesh.verts[i]->u[0], mesh.verts[i]->u[1]));
	}
	GEOM_FADE2D::Fade_2D dt;
	vector<GEOM_FADE2D::Triangle2*> vAllTris;

	dt.insert(vInputPoints);

	vector<GEOM_FADE2D::Segment2> vSegment;
	for (int i = 0; i < mesh.edges.size(); i++) {
		if(mesh.edges[i]->preserve) vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.edges[i]->n[0]->index], vInputPoints[mesh.edges[i]->n[1]->index]));
		//vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.edges[i]->n[0]->index], vInputPoints[mesh.edges[i]->n[1]->index]));
	}
	GEOM_FADE2D::ConstraintGraph2* pCG = dt.createConstraint(vSegment, GEOM_FADE2D::CIS_CONSTRAINED_DELAUNAY);
	dt.applyConstraintsAndZones();

	// display remesh for debugging
	GEOM_FADE2D::Visualizer2 vis2(file_name);

	for (int i = 0; i < vInputPoints.size(); ++i)
	{
		std::string text = std::to_string(i);
		vis2.addObject(GEOM_FADE2D::Label(vInputPoints[i], text, false, 2), GEOM_FADE2D::Color(0, 0, 1, 0.01));
	}

	dt.show(&vis2);
	vis2.writeFile();
}

bool collapse_black_edges(Mesh& mesh, double thresh, bool& repeat)
{
	//cout << "starting black edge collapse" << endl;
	for (int i = 0; i < mesh.edges.size(); i++) {
		if (mesh.edges[i]->preserve == true) {
			for (int k = 0; k < 2; k++) {
				//cout << "faces: " << mesh.edges[i]->n[k]->verts[0]->adjf.size() << endl;
				for (int l = 0; l < mesh.edges[i]->n[k]->verts[0]->adjf.size(); l++) {
					for (int j = 0; j < 3; j++) {
						//cout << mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->index << " -> " << mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->index << ": " << calc_edge_length(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << endl;
						if (edge_length(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) < thresh && mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->preserve == false) {
							if (!((is_seam_or_boundary(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]) && is_seam_or_boundary(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1])) || 
								(!is_seam_or_boundary(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]) && !is_seam_or_boundary(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1])))) {
								continue;
							}
							//if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->preserve_once || mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->preserve_once) continue;
							//if (norm(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->v) > 2.0) continue;
							//cout << norm(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->v) << " : " << norm(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->v) << endl;
							RemeshOp op;
							if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->index == 175 || mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->index == 175) {
								cout << endl;
							}
							if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->on_preserved_edge) {
								//if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->verts[0]->u[1] < mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->verts[0]->u[1] && norm(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->v) > 2.0) continue;
								if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->preserve) { // Preserve point check for cloth corners FIX
									Node* n = mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1];
									double old_edge_num = mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->which_edge;
									Vec3 old_weight = mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->verts[0]->egde_weight;
									op = collapse_edgeP(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j], 0);
									if (op.empty()) {
										//cout << "op empty" << endl;
										repeat = true;
										continue;
									}
									n->EoL = true;
									n->on_preserved_edge = true;
									n->which_edge = old_edge_num;
									n->verts[0]->egde_weight = old_weight;

								}
								else {
									op = collapse_edgeP(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j], 1);
									if (op.empty()) {
										//cout << "op empty" << endl;
										repeat = true;
										continue;
									}
								}
								//cout << "BLACK1: X-> " << edge_length(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << " x-> " << edge_length3D(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << endl;
								op.update(mesh.edges);
								op.done();
								//reindex_nodes(mesh.nodes);

								
								return true;
							}
							else if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->on_preserved_edge) {
								//if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->verts[0]->u[1] < mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->verts[0]->u[1] && norm(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->v) > 2.0) continue;
								if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->preserve) { // Preserve point check for cloth corners FIX
									Node* n = mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0];
									double old_edge_num = mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->which_edge;
									Vec3 old_weight = mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->verts[0]->egde_weight;
									op = collapse_edgeP(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j], 1);
									if (op.empty()) {
										//cout << "op empty" << endl;
										repeat = true;
										continue;
									}
									n->EoL = true;
									n->on_preserved_edge = true;
									n->which_edge = old_edge_num;
									n->verts[0]->egde_weight = old_weight;
								}
								else {
									op = collapse_edgeP(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j], 0);
									if (op.empty()) {
										//cout << "op empty" << endl;
										repeat = true;
										continue;
									}
								}
								//cout << "BLACK2: X-> " << edge_length(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << " x-> " << edge_length3D(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << endl;
								op.update(mesh.edges);
								op.done();
								//reindex_nodes(mesh.nodes);
								
								return true;
							}
							else {
								if (!mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[0]->preserve) {
									//cout << calc_edge_length(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << endl;
									op = collapse_edgeP(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j], 0);
									if (op.empty()) {
										if (mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->preserve) continue;
										op = collapse_edgeP(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j], 1);
										if (op.empty()) {
											//cout << "op empty" << endl;
											repeat = true;
											continue;
										}
									}
									//cout << "BLACK3: X-> " << edge_length(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << " x-> " << edge_length3D(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << endl;
									op.update(mesh.edges);
									op.done();
									//reindex_nodes(mesh.nodes);

									
									return true;
								}
								else if (!mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]->n[1]->preserve) {
									op = collapse_edgeP(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j], 1);
									if (op.empty()) {
										//cout << "op empty" << endl;
										repeat = true;
										continue;
									}
									//cout << "BLACK4: X-> " << edge_length(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << " x-> " << edge_length3D(mesh.edges[i]->n[k]->verts[0]->adjf[l]->adje[j]) << endl;
									op.update(mesh.edges);
									op.done();
									//reindex_nodes(mesh.nodes);

									
									return true;
								}
							}
						}
					}
				}
			}
		}
	}
	return false;
}

bool collapse_black_edges_off_points(Mesh& mesh, double thresh, bool& repeat)
{
	for (int i = 0; i < mesh.nodes.size(); i++) {
		if (mesh.nodes[i]->EoL) {
			bool alone = true;
			//cout << "eol: " << mesh.nodes[i]->index << endl;
			for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
				if (mesh.nodes[i]->adje[j]->preserve) {
					alone = false;
					break;
				}
			}
			if (!alone) continue;
			//cout << "ALONE BAD POINT" << endl;
			for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
				if (edge_length(mesh.nodes[i]->adje[j]) < thresh) {
					RemeshOp op;
					if (mesh.nodes[i] == mesh.nodes[i]->adje[j]->n[0]) op = collapse_edgeP(mesh.nodes[i]->adje[j], 1);
					else op = collapse_edgeP(mesh.nodes[i]->adje[j], 0);
					if (op.empty()) continue;
					op.update(mesh.edges);
					op.done();
					//reindex_nodes(mesh.nodes);
					cout << "ALone EOL collapse" << endl;
					return true;
				}
			}
		}
	}
	return false;
}

bool collapse_green_edges(Mesh& mesh, double thresh, bool& repeat)
{
	for (int i = 0; i < mesh.edges.size(); i++) {
		if (mesh.edges[i]->preserve == true) {
			//cout << calc_edge_length(mesh.edges[i]) << endl;
			if (edge_length(mesh.edges[i]) < thresh) {
				RemeshOp op;
				//if (mesh.edges[i]->n[0]->index == 188 || mesh.edges[i]->n[1]->index == 188) {
				//	cout << endl;
				//}
				if (is_seam_or_boundary(mesh.edges[i]->n[1]) || mesh.edges[i]->n[1]->on_corner) {
					op = collapse_edgeP(mesh.edges[i], 0);
					if (op.empty()) continue;
					op.update(mesh.edges);
					op.done();
					//reindex_nodes(mesh.nodes);
					return true;
				}
				else if (is_seam_or_boundary(mesh.edges[i]->n[0]) || mesh.edges[i]->n[0]->on_corner) {
					op = collapse_edgeP(mesh.edges[i], 1);
					if (op.empty()) continue;
					op.update(mesh.edges);
					op.done();
					//reindex_nodes(mesh.nodes);
					return true;
				}
				else {
					// Use less busy node
					if (mesh.edges[i]->n[0]->verts[0]->adjf.size() <= mesh.edges[i]->n[1]->verts[0]->adjf.size()) {
						op = collapse_edgeP(mesh.edges[i], 1);
						if (op.empty()) {
							op = collapse_edgeP(mesh.edges[i], 0);
							if (op.empty()) {
								repeat = true;
								continue;
							}
						}
					}
					else if (mesh.edges[i]->n[0]->verts[0]->adjf.size() > mesh.edges[i]->n[1]->verts[0]->adjf.size()) {
						op = collapse_edgeP(mesh.edges[i], 1);
						if (op.empty()) {
							op = collapse_edgeP(mesh.edges[i], 0);
							if (op.empty()) {
								repeat = true;
								continue;
							}
						}
					}

					//op = collapse_edgeP(mesh.edges[i], 1);
					//if (op.empty()) {
					//	op = collapse_edgeP(mesh.edges[i], 0);
					//	if (op.empty()) {
					//		repeat = true;
					//		continue;
					//	}
					//}
					op.update(mesh.edges);
					op.done();
					//reindex_nodes(mesh.nodes);
					return true;
				}
			}
		}
	}
	return false;
}

bool remove_awkward_faces(Mesh& mesh)
{
	for (int i = 0; i < mesh.faces.size(); i++) {
		//cout << aspect(mesh.faces[i]) << endl;
		if (aspect(mesh.faces[i]) < 0.00001) {
			cout << "still happening" << endl;
			//double largest_e_length = -1;
			//int whichj = -1;
			//Edge *e;
			//for (int j = 0; j < 3; j++) {
			//	//if ((is_seam_or_boundary(mesh.faces[i]->adje[j]->n[0]) && is_seam_or_boundary(mesh.faces[i]->adje[j]->n[1])) || (!is_seam_or_boundary(mesh.faces[i]->adje[j]->n[0]) && !is_seam_or_boundary(mesh.faces[i]->adje[j]->n[1])))
			//	//{
			//	//	largest_e_length = -1;
			//	//	break;
			//	//}
			//	if (edge_length(mesh.faces[i]->adje[j]) > largest_e_length) {
			//		largest_e_length = edge_length(mesh.faces[i]->adje[j]);
			//		//whichj = mesh.faces[i]->adje[j]->index;
			//		e = mesh.faces[i]->adje[j];
			//	}
			//}
			//if (largest_e_length == -1); continue;

			////mesh.remove(mesh.faces[i]->adje[whichj]);
			//mesh.remove(mesh.faces[i]);
			//if (e->adjf[0] || e->adjf[1]) {
			//	//while (remove_awkward_faces());
			//}
			//mesh.remove(e);
			//return true;
		}
	}
	return false;
}

Vert *adjacent_vert2(const Node *node, const Vert *vert) {
	const Edge *edge = get_edge(node, vert->node);
	for (int i = 0; i < 2; i++)
		for (int s = 0; s < 2; s++)
			if (edge_vert(edge, s, i) == vert)
				return edge_vert(edge, s, 1 - i);
	return NULL;
}

double face_altitude(Edge* edge, Face* face) {
	return (2 * area(face)) / edge_length(edge);
}

bool remove_skiiny_triangles(Mesh& mesh, bool& repeat,bool eolcheck)
{
	vector<Edge*> bad_edges;
	//for (int i = 0; i < mesh.faces.size(); i++) {
	//	if (aspect(mesh.faces[i]) < 0.1) {
	//		//cout << "face: " << i << " ";
	//		//for (int j = 0; j < 3; j++) {
	//		//	cout << mesh.faces[i]->adje[j]->preserve << " ";
	//		//}
	//		//cout << ", ";
	//		//cout.precision(17);
	//		//for (int j = 0; j < 3; j++) {
	//		//	cout << mesh.faces[i]->v[j]->index << ": " << fixed << mesh.faces[i]->v[j]->u[1] << " ";
	//		//}
	//		//cout << endl;
	//		for (int j = 0; j < 3; j++) {
	//			if (mesh.faces[i]->adje[j]->preserve) {
	//				bad_edges.push_back(mesh.faces[i]->adje[j]);
	//				break;
	//			}
	//		}
	//	}
	//}

	for (int i = 0; i < mesh.edges.size(); i++) {
		//if (mesh.edges[i]->preserve) {
			if (mesh.edges[i]->adjf[0] != NULL && face_altitude(mesh.edges[i], mesh.edges[i]->adjf[0]) < 0.004) {
				bad_edges.push_back(mesh.edges[i]);
				//cout << mesh.edges[i]->adjf[0]->v[0]->index << " : " << mesh.edges[i]->adjf[0]->v[1]->index << " : " << mesh.edges[i]->adjf[0]->v[2]->index << endl;
			}
			else if (mesh.edges[i]->adjf[1] != NULL && face_altitude(mesh.edges[i], mesh.edges[i]->adjf[1]) < 0.004) {
				bad_edges.push_back(mesh.edges[i]);
			}
		//}
	}

	cout << "bad edges size: " << bad_edges.size() << endl;
	//for (int i = 0; i < bad_edges.size(); i++) {
	//	for (int j = i; j < bad_edges.size(); j++) {
	//		if (bad_edges[i]->n[0] == bad_edges[j]->n[0] || bad_edges[i]->n[0] == bad_edges[j]->n[1] ||
	//			bad_edges[i]->n[1] == bad_edges[j]->n[0] || bad_edges[i]->n[1] == bad_edges[j]->n[1]) {
	//			Face *f0;
	//			Edge *e0;
	//			Vert *v1, *v2;
	//			if (bad_edges[i]->n[0] == bad_edges[j]->n[0]) {
	//				v1 = bad_edges[i]->n[1]->verts[0];
	//				v2 = bad_edges[j]->n[1]->verts[0];
	//			}
	//			else if (bad_edges[i]->n[0] == bad_edges[j]->n[1]) {
	//				v1 = bad_edges[i]->n[1]->verts[0];
	//				v2 = bad_edges[j]->n[0]->verts[0];
	//			}
	//			else if (bad_edges[i]->n[1] == bad_edges[j]->n[0]) {
	//				v1 = bad_edges[i]->n[0]->verts[0];
	//				v2 = bad_edges[j]->n[1]->verts[0];
	//			}
	//			else if (bad_edges[i]->n[1] == bad_edges[j]->n[1]) {
	//				v1 = bad_edges[i]->n[0]->verts[0];
	//				v2 = bad_edges[j]->n[0]->verts[0];
	//			}
	//			if (bad_edges[i]->adjf[0] == bad_edges[j]->adjf[0]) f0 = bad_edges[i]->adjf[0];
	//			else if (bad_edges[i]->adjf[0] == bad_edges[j]->adjf[1]) f0 = bad_edges[i]->adjf[0];
	//			else if (bad_edges[i]->adjf[1] == bad_edges[j]->adjf[0]) f0 = bad_edges[i]->adjf[1];
	//			else if (bad_edges[i]->adjf[1] == bad_edges[j]->adjf[1]) f0 = bad_edges[i]->adjf[1];
	//			if (f0->adje[0] != bad_edges[i] || f0->adje[0] != bad_edges[j]) e0 = f0->adje[0];
	//			else if (f0->adje[1] != bad_edges[i] || f0->adje[1] != bad_edges[j]) e0 = f0->adje[1];
	//			else if (f0->adje[2] != bad_edges[i] || f0->adje[2] != bad_edges[j]) e0 = f0->adje[2];
	//			double weight = find_split_weight(e0, v1, v2);
	//			cout << "weight: " << weight << endl;
	//			if (weight < 0.0 || weight > 1.0) continue;
	//			Edge* to_collapse1;
	//			Edge* to_collapse2;
	//			Node *node0 = e0->n[0], *node1 = e0->n[1];
	//			RemeshOp op = split_edgeP(e0, weight);
	//			for (size_t v = 0; v < op.added_verts.size(); v++) {
	//				Vert *vertnew = op.added_verts[v];
	//				Vert *v00 = adjacent_vert2(node0, vertnew),
	//					*v11 = adjacent_vert2(node1, vertnew);
	//				vertnew->sizing = weight * (v00->sizing + v11->sizing);
	//				if (eolcheck) {
	//					vertnew->node->EoL = true;
	//					vertnew->node->EoL_secondary = true;
	//				}
	//				vertnew->node->on_preserved_edge = true;
	//				vertnew->egde_weight[0] = calc_edge_weight(vertnew->node);
	//				// For immedietely collapsing small edge
	//				for (int j = 0; j < vertnew->node->adje.size(); j++) {
	//					if (vertnew->node->adje[j]->n[0]->verts[0] == v1 || vertnew->node->adje[j]->n[1]->verts[0] == v1) {
	//						to_collapse1 = vertnew->node->adje[j];
	//					}
	//					else if (vertnew->node->adje[j]->n[0]->verts[0] == v2 || vertnew->node->adje[j]->n[1]->verts[0] == v2) {
	//						to_collapse2 = vertnew->node->adje[j];
	//					}
	//				}
	//			}
	//			op.set_null(bad_edges);
	//			op.done();
	//			repeat = true;
	//			return true;
	//		}
	//	}
	//}
	for (int i = 0; i < bad_edges.size(); i++) {
		//cout << bad_edges[i]->n[0]->index << " : " << bad_edges[i]->n[0]->EoL << " -> " << bad_edges[i]->n[1]->index << bad_edges[i]->n[1]->EoL << endl;
		cout << bad_edges[i]->adjf[0]->v[0]->index << " : " << bad_edges[i]->adjf[0]->v[1]->index << " : " << bad_edges[i]->adjf[0]->v[2]->index << " : " << endl;
		cout << bad_edges[i]->adjf[1]->v[0]->index << " : " << bad_edges[i]->adjf[1]->v[1]->index << " : " << bad_edges[i]->adjf[1]->v[2]->index << " : " << endl;
		repeat = true;
		Edge *edge = bad_edges[i];
		Node *node0 = edge->n[0], *node1 = edge->n[1];
		Vert *v1, *v2;
		for (int j = 0; j < 3; j++) {
			if (node0->verts[0] != edge->adjf[0]->v[j] && node1->verts[0] != edge->adjf[0]->v[j]) {
				v1 = edge->adjf[0]->v[j];
				break;
			}
		}
		for (int j = 0; j < 3; j++) {
			if (node0->verts[0] != edge->adjf[1]->v[j] && node1->verts[0] != edge->adjf[1]->v[j]) {
				v2 = edge->adjf[1]->v[j];
				break;
			}
		}
		if (is_seam_or_boundary(v1->node) && is_seam_or_boundary(v2->node)) {
			repeat = false;
			continue;
		}
		if (is_seam_or_boundary(v1) && is_seam_or_boundary(v2)) continue;
		double weight = find_split_weight(edge, v1, v2);
		cout << "weight: " << weight << endl;
		if (weight < 0.0 || weight > 1.0) continue;
		Edge* to_collapse1;
		Edge* to_collapse2;
		RemeshOp op = split_edgeP(edge, weight);
		for (size_t v = 0; v < op.added_verts.size(); v++) {
			Vert *vertnew = op.added_verts[v];
			Vert *v00 = adjacent_vert2(node0, vertnew),
				*v11 = adjacent_vert2(node1, vertnew);
			vertnew->sizing = weight * (v00->sizing + v11->sizing);
			if (eolcheck) {
				vertnew->node->EoL = true;
				vertnew->node->EoL_secondary = true;
			}
			vertnew->node->on_preserved_edge = true;
			vertnew->egde_weight[0] = calc_edge_weight(vertnew->node);
			// For immedietely collapsing small edge
			for (int j = 0; j < vertnew->node->adje.size(); j++) {
				if (vertnew->node->adje[j]->n[0]->verts[0] == v1 || vertnew->node->adje[j]->n[1]->verts[0] == v1) {
					to_collapse1 = vertnew->node->adje[j];
				}
				else if (vertnew->node->adje[j]->n[0]->verts[0] == v2 || vertnew->node->adje[j]->n[1]->verts[0] == v2) {
					to_collapse2 = vertnew->node->adje[j];
				}
			}
		}
		op.set_null(bad_edges);
		op.done();

		//break; // One at a time

		// For immedietely collapsing small edge
		//RemeshOp op2;
		//if (edge_length(to_collapse1) < edge_length(to_collapse2)) {
		//	if (!is_seam_or_boundary(to_collapse1->n[0]) && !is_seam_or_boundary(to_collapse1->n[1])) {
		//		if (to_collapse1->n[0]->on_preserved_edge) {
		//			op = collapse_edgeP(to_collapse1, 1);
		//		}
		//		else {
		//			op = collapse_edgeP(to_collapse1, 0);
		//		}
		//	}
		//}
		//else {
		//	if (!is_seam_or_boundary(to_collapse2->n[0]) && !is_seam_or_boundary(to_collapse2->n[1])) {
		//		if (to_collapse2->n[0]->on_preserved_edge) {
		//			op = collapse_edgeP(to_collapse2, 1);
		//		}
		//		else {
		//			op = collapse_edgeP(to_collapse2, 0);
		//		}
		//	}
		//}
		//op2.update(mesh.edges);
		//op2.done();
	}
	//if (repeat) {
	//	export2D(mesh, "altidude_collapse.ps");
	//}

	return false;
}

void nudge_check(Mesh& mesh, const double &thresh, const Vector2d &bounds)
{
	Node *n = mesh.nodes.back();
	int size = mesh.nodes.size() - 1;
	for (int i = 0; i < size; i++) {
		if (mesh.nodes[i]->on_preserved_edge) continue;
		int which = -1;
		//if (unsigned_vv_distance(n->verts[0]->u, mesh.nodes[i]->verts[0]->u) < thresh) {
		if ((n->verts[0]->u[1] - mesh.nodes[i]->verts[0]->u[1]) < thresh) {
			cout << "NUDGE CHECK: " << n->index << " -> " << mesh.nodes[i]->index << endl;
			Node *n0 = mesh.nodes[i];

			double sign = n0->verts[0]->u[1] - n->verts[0]->u[1] > 0.0 ? 1.0 : -1.0;
			if (sign < 0.0) continue;
			Vec3 newX(n0->verts[0]->u[0], n->verts[0]->u[1] + sign *(thresh + 1e-2), 0);
			bool toclose = false;
			//cout << "newx " << newX << endl;
			for (int j = 0; j < mesh.nodes.size(); j++) {
				if (unsigned_vv_distance(newX, mesh.nodes[j]->verts[0]->u) < 2.0 * thresh) {
					toclose = true;
					break;
				}
			}
			if (toclose) continue;
			Face* f0 = get_enclosing_face(mesh, Vec2(newX[0],newX[1]));
			Vec3 bary = get_barycentric_coords(Vec2(newX[0],newX[1]), f0);
			//cout << "BARY : " << f0->v[0]->index << " " << f0->v[1]->index << " " << f0->v[2]->index << " " << bary << endl;
			Vec3 tx;
			for (int j = 0; j < 3; j++) {
				tx[j] = bary[0] * f0->v[0]->node->x[j] +
					bary[1] * f0->v[1]->node->x[j] + 
					bary[2] * f0->v[2]->node->x[j];
			}
			Vec3 tv(0);
			for (int j = 0; j < 3; j++) {
				tv[j] = bary[0] * f0->v[0]->node->v[j] +
					bary[1] * f0->v[1]->node->v[j] +
					bary[2] * f0->v[2]->node->v[j];
			}
			Vec3 tV(0);
			for (int j = 0; j < 3; j++) {
				tV[j] = bary[0] * f0->v[0]->v[j] +
					bary[1] * f0->v[1]->v[j] +
					bary[2] * f0->v[2]->v[j];
			}
			if (newX[0] < 1e-12) newX[0] = 0.0;
			if (bounds[0] - newX[0] < 1e-12) newX[0] = bounds[0];
			if (newX[1] < 1e-12) newX[1] = 0.0;
			if (bounds[1] - newX[1] < 1e-12) newX[1] = bounds[1];
			mesh.add(new Vert(newX, tV));
			mesh.add(new Node(tx, tx, Vec3(0), 0, 0, false));
			mesh.nodes.back()->v = tv;
			mesh.nodes.back()->preserve = true;
			connect(mesh.verts.back(), mesh.nodes.back());
			cout << "NEW NUDGE: " << mesh.nodes.back()->index << " from " << n0->x << " -> " << tx << endl;
		}
	}
	//for (int i = 0; i < mesh.edges.size(); i++) {
	//	if (mesh.edges[i]->preserve) {
	//		Node *n0;
	//		Vert *v0;
	//		if (face_altitude(mesh.edges[i], mesh.edges[i]->adjf[0]) < 0.01) {
	//			v0 = edge_opp_vert(mesh.edges[i], 0);
	//			n0 = v0->node;
	//		}
	//		else if (face_altitude(mesh.edges[i], mesh.edges[i]->adjf[1]) < 0.01) {
	//			v0 = edge_opp_vert(mesh.edges[i], 1);
	//			n0 = v0->node;
	//		}
	//		else continue;
	//		cout << "NUDGE CHECK FROM EDGE: " << n0->index << "OFF " << mesh.edges[i]->n[0]->index << " " << mesh.edges[i]->n[1]->index << endl;
	//		double middley = (mesh.edges[i]->n[0]->verts[0]->u[1] + mesh.edges[i]->n[1]->verts[0]->u[1]) / 2;
	//		double sign = n0->verts[0]->u[1] - middley > 0.0 ? 1.0 : -1.0;
	//		if (sign < 0.0) continue;
	//		Vec3 newX(n0->verts[0]->u[0], middley + sign *(thresh + 1e-2), 0);
	//		bool toclose = false;
	//		//cout << "newx " << newX << endl;
	//		for (int j = 0; j < mesh.nodes.size(); j++) {
	//			if (unsigned_vv_distance(newX, mesh.nodes[j]->verts[0]->u) < thresh) {
	//				toclose = true;
	//				break;
	//			}
	//		}
	//		if (toclose) continue;
	//		Face* f0 = get_enclosing_face(mesh, Vec2(newX[0], newX[1]));
	//		Vec3 bary = get_barycentric_coords(Vec2(newX[0], newX[1]), f0);
	//		//cout << "BARY : " << f0->v[0]->index << " " << f0->v[1]->index << " " << f0->v[2]->index << " " << bary << endl;
	//		Vec3 tx;
	//		for (int j = 0; j < 3; j++) {
	//			tx[j] = bary[0] * f0->v[0]->node->x[j] +
	//				bary[1] * f0->v[1]->node->x[j] +
	//				bary[2] * f0->v[2]->node->x[j];
	//		}
	//		Vec3 tv(0);
	//		for (int j = 0; j < 3; j++) {
	//			tv[j] = bary[0] * f0->v[0]->node->v[j] +
	//				bary[1] * f0->v[1]->node->v[j] +
	//				bary[2] * f0->v[2]->node->v[j];
	//		}
	//		Vec3 tV(0);
	//		for (int j = 0; j < 3; j++) {
	//			tV[j] = bary[0] * f0->v[0]->v[j] +
	//				bary[1] * f0->v[1]->v[j] +
	//				bary[2] * f0->v[2]->v[j];
	//		}
	//		if (newX[0] < 1e-12) newX[0] = 0.0;
	//		if (bounds[0] - newX[0] < 1e-12) newX[0] = bounds[0];
	//		if (newX[1] < 1e-12) newX[1] = 0.0;
	//		if (bounds[1] - newX[1] < 1e-12) newX[1] = bounds[1];
	//		mesh.add(new Vert(newX, tV));
	//		mesh.add(new Node(tx, tx, Vec3(0), 0, 0, false));
	//		mesh.nodes.back()->v = tv;
	//		mesh.nodes.back()->preserve = true;
	//		connect(mesh.verts.back(), mesh.nodes.back());
	//		cout << "NEW NUDGE: " << mesh.nodes.back()->index << " from " << n0->x << " -> " << tx << endl;
	//	}
	//}
}

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
int edgeVerts1_data[] = {
	0, 1,
	2, 0,
	1, 3,
	3, 2,
	0, 4,
	5, 1,
	4, 5,
	6, 2,
	4, 6,
	3, 7,
	7, 5,
	6, 7,
};
Map<Matrix<double, 4, 14, ColMajor> > verts1_(verts1_data);
Map<Matrix<int, 2, 12, ColMajor> > edgeVerts1(edgeVerts1_data);
Matrix<double, 4, 14> verts1BOX;

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

void updateEdgeWeights(Mesh& mesh, vector<shared_ptr<Box> > box)
{
	for (int i = 0; i < mesh.nodes.size(); i++) {
		if (mesh.nodes[i]->EoL && !mesh.nodes[i]->on_corner) {
			int which_edge = mesh.nodes[i]->which_edge;
			int which_box = which_edge / 12;
			Matrix4d S = Matrix4d::Identity();
			S(0, 0) = 0.5*box[which_box]->dim(0);
			S(1, 1) = 0.5*box[which_box]->dim(1);
			S(2, 2) = 0.5*box[which_box]->dim(2);
			Matrix4d E = box[which_box]->E1 * S;
			verts1BOX = E * verts1_;
			Vector2i indices = edgeVerts1.col(which_edge % 12);
			Vector3d A = verts1BOX.block<3, 1>(0, indices(0));
			Vector3d B = verts1BOX.block<3, 1>(0, indices(1));
			Vector3d P = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
			double weight = linepoint(A, B, P);
			mesh.verts[i]->egde_weight[0] = 1.0-weight;
			//cout << "A: " << A << endl;
			//cout << "B: " << B << endl;
			//cout << "P: " << P << endl;
			//cout << "Before: " << mesh.verts[i]->egde_weight << endl;
			//cout << "After: " << 1.0-weight << endl;
		}
	}
}

bool collisionRemesh(Mesh& mesh, vector<shared_ptr<Box> > box, vector<std::shared_ptr<btc::Collision> > collisions, VectorXi EOLS, shared_ptr<pp_settings> pps)
{
	// Assume function is only called if collision is on

	// Box Setup
	//collisions.clear();
	MatrixXd verts2(3, mesh.nodes.size());
	MatrixXi faces2(3, mesh.faces.size());
	vector<vector<shared_ptr<CVM> > > box_edges_tracker;
	box_edges_tracker.resize(12 * box.size());

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

	for (int i = 0; i < mesh.nodes.size(); i++) {
		verts2.col(i) = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
		if (mesh.nodes[i]->EoL) {
			//cout << "Index: " << mesh.nodes[i]->index << "on preserve: " << mesh.nodes[i]->on_preserved_edge << endl;
			if (mesh.nodes[i]->on_corner) {
				//if (corner_not_within_safety_margin(mesh.nodes[i]->verts[0]->u, pps->safety_margin)) {
				//	cout << "OUT" << endl;
				//	for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
				//		mesh.nodes[i]->adje[j]->preserve = false;
				//	}
				//	mesh.nodes[i]->preserve = false;
				//	mesh.nodes[i]->EoL = false;
				//	mesh.nodes[i]->on_corner = false;
				//	mesh.nodes[i]->on_preserved_edge = false;
				//	continue;
				//}
				box_edges_tracker[mesh.nodes[i]->corner_egdes[0]].push_back(make_shared<CVM>(mesh.verts[i], 0));
				if (!pps->points) {
					box_edges_tracker[mesh.nodes[i]->corner_egdes[1]].push_back(make_shared<CVM>(mesh.verts[i], 1));
					box_edges_tracker[mesh.nodes[i]->corner_egdes[2]].push_back(make_shared<CVM>(mesh.verts[i], 2));
				}
				continue;
			}
			else {
				box_edges_tracker[mesh.nodes[i]->which_edge].push_back(make_shared<CVM>(mesh.verts[i]));
				continue;
			}
		}
		else {
			if (mesh.nodes[i]->preserve_once && false) {
				mesh.nodes[i]->preserve_once = false;
			}
			else {
				for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
					mesh.nodes[i]->adje[j]->preserve = false;
				}
			}
		}
		//mesh.nodes[i]->EoL = false;
	}

	for (int i = 0; i < mesh.faces.size(); i++) {
		faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}

	int old_max_index = mesh.nodes.size() - 1;
	int newly_added = 0;
	Material* material = mesh.faces[0]->material; // Cleanup

	vector<int> saved_edges;

	pps->PPTimer[0]->tic();
	for (int b = 0; b < box.size(); b++) {
		collisions.clear();
		if(pps->wire) boxTriCollisionHack2(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, EOLS, false);
		else if (pps->points) {
			pointTriCollision(collisions, box[b]->thresh, pps->verts1, pps->norms1, verts2, faces2, false);
			if(pps->move_points) pps->verts1 += pps->norms1 * pps->pmove1 * pps->h;
			cout << collisions.size() << endl;
		}
		else {
			boxTriCollision(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, EOLS, false);
		}
		if (collisions.size() == 0) continue;
		if (pps->matlab_debug_collision) {
			double_to_file(box[b]->thresh, "threshold");
			vec_to_file(box[b]->dim, "whd1");
			mat_to_file(box[b]->E1, "E1");
			mat_to_file(verts2, "verts2");
			VectorXi vvv(3);
			vvv << 1, 1, 1;
			vec_to_file(EOLS, "isEOL2");
			mat_to_file(faces2.colwise() += vvv, "faces2");
		}

		//int numedgecoll = 0;
		//for (int i = 0; i < collisions.size(); i++) {
		//	if (collisions[i]->count1 == 2 && collisions[i]->count2 == 2) {
		//		numedgecoll++;
		//	}
		//}
		//cout << "Edge collisions: " << numedgecoll << endl;
		//cout << collisions.size() << endl;

		//vector<shared_ptr<CVM> > box_edges_tracker[12];
		//bool has_corner_case = false;

		for (int i = 0; i < collisions.size(); i++) {
			if (collisions[i]->count1 == 1 && collisions[i]->count2 == 3) {
				//has_corner_case = true;
				//if (pps->once_hack) continue; // Hack
				if (mesh.nodes[collisions[i]->verts2(0)]->EoL) {
					//bool dontadd = false;
					//for (int j = 0; j < box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->which_edge].size(); j++) {
					//	if (box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->which_edge][j]->vert == mesh.verts[collisions[i]->verts2(0)]) dontadd = true;
					//}
					//if (dontadd) continue;
					//box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->which_edge].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)]));
					continue;
				}
				if (mesh.nodes[collisions[i]->verts2(1)]->EoL) {
					//bool dontadd = false;
					//for (int j = 0; j < box_edges_tracker[mesh.nodes[collisions[i]->verts2(1)]->which_edge].size(); j++) {
					//	if (box_edges_tracker[mesh.nodes[collisions[i]->verts2(1)]->which_edge][j]->vert == mesh.verts[collisions[i]->verts2(1)]) dontadd = true;
					//}
					//if (dontadd) continue;
					//box_edges_tracker[mesh.nodes[collisions[i]->verts2(1)]->which_edge].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(1)]));
					continue;
				}
				if (mesh.nodes[collisions[i]->verts2(2)]->EoL) {
					//bool dontadd = false;
					//for (int j = 0; j < box_edges_tracker[mesh.nodes[collisions[i]->verts2(2)]->which_edge].size(); j++) {
					//	if (box_edges_tracker[mesh.nodes[collisions[i]->verts2(2)]->which_edge][j]->vert == mesh.verts[collisions[i]->verts2(2)]) dontadd = true;
					//}
					//if (dontadd) continue;
					//box_edges_tracker[mesh.nodes[collisions[i]->verts2(2)]->which_edge].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(2)]));
					continue;
				}
				double tx = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[0] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[0] +
					collisions[i]->weights2(2) * mesh.verts[collisions[i]->verts2(2)]->u[0];
				if (tx < 1e-12) tx = 0.0;
				if (pps->bounds[0] - tx < 1e-12) tx = pps->bounds[0];
				double ty = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[1] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[1] +
					collisions[i]->weights2(2) * mesh.verts[collisions[i]->verts2(2)]->u[1];
				if (ty < 1e-12) ty = 0.0;
				if (pps->bounds[1] - ty < 1e-12) ty = pps->bounds[1];
				Vec3 tv(0);
				for (int j = 0; j < 3; j++) {
					tv[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->node->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->node->v[j] +
						collisions[i]->weights2(2) * mesh.verts[collisions[i]->verts2(2)]->node->v[j];
				}
				Vec3 tV(0);
				for (int j = 0; j < 3; j++) {
					tV[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->v[j] +
						collisions[i]->weights2(2) * mesh.verts[collisions[i]->verts2(2)]->v[j];
				}
				Vec3 t(tx, ty, 0);
				if (corner_not_within_safety_margin(t, pps->safety_margin, pps->bounds)) continue; // Don't add collisions within the corner safety margin
																						//Vec3 n(collisions[i]->pos2(0), collisions[i]->pos2(1), collisions[i]->pos2(2)); // pos2
				Vec3 n(collisions[i]->pos1_(0), collisions[i]->pos1_(1), collisions[i]->pos1_(2)); // pos1_
				mesh.add(new Vert(t, tV));
				mesh.add(new Node(n, n, Vec3(0), 0, 0, false));
				mesh.nodes.back()->v = tv;
				connect(mesh.verts.back(), mesh.nodes.back());
				box_edges_tracker[collisions[i]->cornerEdges(0) + (12 * b)].push_back(make_shared<CVM>(mesh.verts.back(), 0));
				if (!pps->points) {
					box_edges_tracker[collisions[i]->cornerEdges(1) + (12 * b)].push_back(make_shared<CVM>(mesh.verts.back(), 1));
					box_edges_tracker[collisions[i]->cornerEdges(2) + (12 * b)].push_back(make_shared<CVM>(mesh.verts.back(), 2));
				}
				mesh.nodes.back()->corner_egdes[0] = collisions[i]->cornerEdges(0) + (12 * b);
				mesh.nodes.back()->corner_egdes[1] = collisions[i]->cornerEdges(1) + (12 * b);
				mesh.nodes.back()->corner_egdes[2] = collisions[i]->cornerEdges(2) + (12 * b);
				mesh.verts.back()->egde_weight[0] = collisions[i]->cornerEdgeWeights(0);
				mesh.verts.back()->egde_weight[1] = collisions[i]->cornerEdgeWeights(1);
				mesh.verts.back()->egde_weight[2] = collisions[i]->cornerEdgeWeights(2);

				// Add collision info to node
				mesh.nodes.back()->has_coll_info = true;
				mesh.nodes.back()->on_corner = true;
				mesh.nodes.back()->preserve = true;
				if (pps->EOLon) mesh.nodes.back()->EoL = true;
				mesh.nodes.back()->coll_norm = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(2));  // FIX
				mesh.nodes.back()->nor_ave = Vec3(collisions[i]->cornerNormal(0), collisions[i]->cornerNormal(1), collisions[i]->cornerNormal(2));

				newly_added++;

			}
			else if (collisions[i]->count1 == 2 && collisions[i]->count2 == 2) {
				//if (pps->once_hack) continue; // Hack
				//if (collisions[i]->edge1 == 11 || collisions[i]->edge1 == 9 || collisions[i]->edge1 == 30 || collisions[i]->edge1 == 39) continue; // Ignore Edge Hack
				if (mesh.nodes[collisions[i]->verts2(0)]->EoL) {
					//bool dontadd = false;
					//for (int j = 0; j < box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->which_edge].size(); j++) {
					//	if (box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->which_edge][j]->vert == mesh.verts[collisions[i]->verts2(0)]) dontadd = true;
					//}
					//if (dontadd) continue;
					//box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->which_edge].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)]));
					//cout << "Before: " << mesh.verts[collisions[i]->verts2(0)]->egde_weight[0] << endl;
					//verts1BOX = box[b]->E1 * verts1_;
					//Vector3d A = verts1BOX.block<3, 1>(0, collisions[i]->verts1(0));
					//Vector3d B = verts1BOX.block<3, 1>(0, collisions[i]->verts1(1));
					//Vector3d P = Vector3d(mesh.nodes[collisions[i]->verts2(0)]->x[0], mesh.nodes[collisions[i]->verts2(0)]->x[1], mesh.nodes[collisions[i]->verts2(0)]->x[2]);
					//double weight = linepoint(A, B, P);
					//mesh.verts[collisions[i]->verts2(0)]->egde_weight[0] = 1.0-weight; // Update the edge weight
					//cout << "After: " << mesh.verts[collisions[i]->verts2(0)]->egde_weight[0] << endl;
					continue;
				}
				if (mesh.nodes[collisions[i]->verts2(1)]->EoL) {
					//bool dontadd = false;
					//for (int j = 0; j < box_edges_tracker[mesh.nodes[collisions[i]->verts2(1)]->which_edge].size(); j++) {
					//	if (box_edges_tracker[mesh.nodes[collisions[i]->verts2(1)]->which_edge][j]->vert == mesh.verts[collisions[i]->verts2(1)]) dontadd = true;
					//}
					//if (dontadd) continue;
					//box_edges_tracker[mesh.nodes[collisions[i]->verts2(1)]->which_edge].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(1)]));
					//cout << "Before: " << mesh.verts[collisions[i]->verts2(1)]->egde_weight[0] << endl;
					//verts1BOX = box[b]->E1 * verts1_;
					//Vector3d A = verts1BOX.block<3, 1>(0, collisions[i]->verts1(0));
					//Vector3d B = verts1BOX.block<3, 1>(0, collisions[i]->verts1(1));
					//Vector3d P = Vector3d(mesh.nodes[collisions[i]->verts2(1)]->x[0], mesh.nodes[collisions[i]->verts2(1)]->x[1], mesh.nodes[collisions[i]->verts2(1)]->x[2]);
					//double weight = linepoint(A, B, P);
					//mesh.verts[collisions[i]->verts2(1)]->egde_weight[0] = 1.0 - weight; // Update the edge weight
					//cout << "After: " << mesh.verts[collisions[i]->verts2(1)]->egde_weight[0] << endl;
					continue;
				}
				double tx = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[0] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[0];
				if (tx < 1e-12) tx = 0.0;
				if (pps->bounds[0] - tx < 1e-12) tx = pps->bounds[0];
				double ty = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[1] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[1];
				if (ty < 1e-12) ty = 0.0;
				if (pps->bounds[1] - ty < 1e-12) ty = pps->bounds[1];
				Vec3 tv(0);
				for (int j = 0; j < 3; j++) {
					tv[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->node->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->node->v[j];
				}
				Vec3 tV(0);
				for (int j = 0; j < 3; j++) {
					tV[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->v[j];
				}
				Vec3 t(tx, ty, 0);
				if (quick_edge_safety_margin_check(t, pps->safety_margin, pps->bounds)) continue; // Don't add collisions within the corner safety margin
																						//Vec3 n(collisions[i]->pos2(0), collisions[i]->pos2(1), collisions[i]->pos2(2)); // pos2
				Vec3 n(collisions[i]->pos1_(0), collisions[i]->pos1_(1), collisions[i]->pos1_(2)); // pos1
				mesh.add(new Vert(t, tV));
				mesh.add(new Node(n, n, Vec3(0), 0, 0, false));
				mesh.nodes.back()->v = tv;
				connect(mesh.verts.back(), mesh.nodes.back());
				box_edges_tracker[collisions[i]->edge1 + (12 * b)].push_back(make_shared<CVM>(mesh.verts.back()));

				// Add collision info to node
				mesh.nodes.back()->has_coll_info = true;
				if (pps->EOLon) mesh.nodes.back()->EoL = true;
				mesh.verts.back()->egde_weight[0] = collisions[i]->weights1(0);
				// Old norm
				//mesh.nodes.back()->coll_norm = Vec3(collisions[i]->nor(0), collisions[i]->nor(1), collisions[i]->nor(2));
				// Ave norm
				mesh.nodes.back()->coll_norm = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(2));

				//nudge_check(mesh, pps->collapse_non_preserve_thresh, pps->bounds);

				newly_added++;
			}
			else if (collisions[i]->count1 == 3 && collisions[i]->count2 == 1) {
				// EOLS
				if (mesh.nodes[collisions[i]->verts2(0)]->EoL) {
					//if (mesh.nodes[collisions[i]->verts2(0)]->on_corner) {
					//	box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->corner_egdes[0]].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)]));
					//	box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->corner_egdes[1]].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)]));
					//	box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->corner_egdes[2]].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)]));
					//	continue;
					//}
					//else {
					//	box_edges_tracker[mesh.nodes[collisions[i]->verts2(0)]->which_edge].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)]));
					//	continue;
					//}
					continue;
				}
				if (collisions[i]->edge1 >= 0) {
					box_edges_tracker[collisions[i]->edge1 + (12 * b)].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)]));
					// Add collision info to node
					mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = true;
					mesh.nodes[collisions[i]->verts2(0)]->coll_norm = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(2));
					mesh.nodes[collisions[i]->verts2(0)]->x = Vec3(collisions[i]->pos1_(0), collisions[i]->pos1_(1), collisions[i]->pos1_(2)); // Reset position to avoid sink
					mesh.nodes[collisions[i]->verts2(0)]->verts[0]->egde_weight[0] = collisions[i]->weights2(0);
				}
			}
			else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 1) {
				box_edges_tracker[collisions[i]->cornerEdges(0) + (12 * b)].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)], 0));
				box_edges_tracker[collisions[i]->cornerEdges(1) + (12 * b)].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)], 1));
				box_edges_tracker[collisions[i]->cornerEdges(2) + (12 * b)].push_back(make_shared<CVM>(mesh.verts[collisions[i]->verts2(0)], 2));
				// Add collision info to node
				mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = true;
				mesh.nodes[collisions[i]->verts2(0)]->on_corner = true;
				mesh.nodes[collisions[i]->verts2(0)]->preserve = true;
				mesh.nodes[collisions[i]->verts2(0)]->coll_norm = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(2));
				mesh.nodes[collisions[i]->verts2(0)]->x = Vec3(collisions[i]->pos1_(0), collisions[i]->pos1_(1), collisions[i]->pos1_(2)); // Reset position to avoid sink
				mesh.nodes[collisions[i]->verts2(0)]->verts[0]->egde_weight[0] = collisions[i]->cornerEdgeWeights(0);
				mesh.nodes[collisions[i]->verts2(0)]->verts[0]->egde_weight[1] = collisions[i]->cornerEdgeWeights(1);
				mesh.nodes[collisions[i]->verts2(0)]->verts[0]->egde_weight[2] = collisions[i]->cornerEdgeWeights(2);
			}
			else if (collisions[i]->count1 == 1 && collisions[i]->count2 == 2) {
				double tx = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[0] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[0];
				if (tx < 1e-12) tx = 0.0;
				if (pps->bounds[0] - tx < 1e-12) tx = pps->bounds[0];
				double ty = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[1] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[1];
				if (ty < 1e-12) ty = 0.0;
				if (pps->bounds[1] - ty < 1e-12) ty = pps->bounds[1];
				Vec3 tv(0);
				for (int j = 0; j < 3; j++) {
					tv[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->node->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->node->v[j];
				}
				Vec3 tV(0);
				for (int j = 0; j < 3; j++) {
					tV[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->v[j];
				}
				Vec3 t(tx, ty, 0);
				if (corner_not_within_safety_margin(t, pps->safety_margin, pps->bounds)) continue; // Don't add collisions within the corner safety margin
																						//Vec3 n(collisions[i]->pos2(0), collisions[i]->pos2(1), collisions[i]->pos2(2)); // pos2
				Vec3 n(collisions[i]->pos1_(0), collisions[i]->pos1_(1), collisions[i]->pos1_(2)); // pos1
				mesh.add(new Vert(t, tV));
				mesh.add(new Node(n, n, Vec3(0), 0, 0, false));
				mesh.nodes.back()->v = tv;
				connect(mesh.verts.back(), mesh.nodes.back());
				box_edges_tracker[collisions[i]->cornerEdges(0) + (12 * b)].push_back(make_shared<CVM>(mesh.verts.back(), 0));
				if (!pps->points) {
					box_edges_tracker[collisions[i]->cornerEdges(1) + (12 * b)].push_back(make_shared<CVM>(mesh.verts.back(), 1));
					box_edges_tracker[collisions[i]->cornerEdges(2) + (12 * b)].push_back(make_shared<CVM>(mesh.verts.back(), 2));
				}
				// Add collision info to node
				mesh.nodes[collisions[i]->verts2(0)]->has_coll_info = true;
				mesh.nodes[collisions[i]->verts2(0)]->on_corner = true;
				mesh.nodes[collisions[i]->verts2(0)]->preserve = true;
				mesh.nodes[collisions[i]->verts2(0)]->coll_norm = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(2));
				mesh.nodes[collisions[i]->verts2(0)]->verts[0]->egde_weight[0] = collisions[i]->cornerEdgeWeights(0);
				mesh.nodes[collisions[i]->verts2(0)]->verts[0]->egde_weight[1] = collisions[i]->cornerEdgeWeights(1);
				mesh.nodes[collisions[i]->verts2(0)]->verts[0]->egde_weight[2] = collisions[i]->cornerEdgeWeights(2);
			}
		}
	}
	pps->PPTimer[0]->toc();

	pps->PPTimer[1]->tic();
	vector<Vert *> outside_margin_v;
	vector<Node *> outside_margin_n;
	// Sort edges, check safety margin, and save order
	for (int i = 0; i < box_edges_tracker.size(); i++) {
		if (box_edges_tracker[i].size() == 0) continue;
		if (box_edges_tracker[i].size() == 1) {
			//continue;
			if (box_edges_tracker[i][0]->vert->node->on_corner) continue;
			if (box_edges_tracker[i][0]->vert->index <= old_max_index) {
				//collisions_passed.push_back(box_edges_tracker[i][j]);
				box_edges_tracker[i][0]->vert->node->EoL = false;
				box_edges_tracker[i][0]->vert->node->on_preserved_edge = false;
				//box_edges_tracker[i][0]->vert->node->lookout_for = true;
				box_edges_tracker[i].erase(box_edges_tracker[i].begin());
			}
			else {
				Vert* vert2 = box_edges_tracker[i][0]->vert;
				Node* node2 = box_edges_tracker[i][0]->vert->node;
				outside_margin_v.push_back(vert2);
				outside_margin_n.push_back(node2);
				mesh.remove(vert2);
				mesh.remove(node2);
				box_edges_tracker[i].erase(box_edges_tracker[i].begin());
			}
			continue;
		}
		sort(box_edges_tracker[i].begin(), box_edges_tracker[i].end(), compare_edge_weights);
		//if (to_flip_corner(box_edges_tracker[i][0]->vert, box_edges_tracker[i][1]->vert, box_edges_tracker[i].back()->vert)) {
		//	moveItemToBack(box_edges_tracker[i], 0);
		//}
		// Check if any of the new edges are in a safety margin
		// Should be structured so nodes aren't added to mesh unless they pass this, instead of deleting them
		for (int j = 0; j + 1 < box_edges_tracker[i].size(); j++) {

			// Duplicate problem
			if (box_edges_tracker[i][j]->vert->u == box_edges_tracker[i][j + 1]->vert->u) {
				cout << "DUPE" << endl;
				if (box_edges_tracker[i][j + 1]->vert->index >= old_max_index) {
					Vert* vert = box_edges_tracker[i][j + 1]->vert;
					Node* node = box_edges_tracker[i][j + 1]->vert->node;
					outside_margin_v.push_back(vert);
					outside_margin_n.push_back(node);
					mesh.remove(vert);
					mesh.remove(node);
					box_edges_tracker[i].erase(box_edges_tracker[i].begin() + j + 1);
				} 
				else if (box_edges_tracker[i][j]->vert->index >= old_max_index) {
					Vert* vert = box_edges_tracker[i][j]->vert;
					Node* node = box_edges_tracker[i][j]->vert->node;
					outside_margin_v.push_back(vert);
					outside_margin_n.push_back(node);
					mesh.remove(vert);
					mesh.remove(node);
					box_edges_tracker[i].erase(box_edges_tracker[i].begin() + j);
				}
				j--;
				continue;
			}

			if (edge_not_within_safety_margin(box_edges_tracker[i][j]->vert, box_edges_tracker[i][j + 1]->vert, pps->safety_margin, pps->bounds)) {
				//pps->safety_margin += 0.008;
				cout << "edge safety alert: " << box_edges_tracker[i][j]->vert->index << endl;
				if (box_edges_tracker[i][j]->vert->index <= old_max_index) {
					//collisions_passed.push_back(box_edges_tracker[i][j]);
					cout << box_edges_tracker[i].size() << endl;
					box_edges_tracker[i][j]->vert->node->EoL = false;
					box_edges_tracker[i][j]->vert->node->on_preserved_edge = false;
					//box_edges_tracker[i][j]->vert->node->lookout_for = true;
					box_edges_tracker[i].erase(box_edges_tracker[i].begin() + j);
				}
				else {
					Vert* vert = box_edges_tracker[i][j]->vert;
					Node* node = box_edges_tracker[i][j]->vert->node;
					outside_margin_v.push_back(vert);
					outside_margin_n.push_back(node);
					mesh.remove(vert);
					mesh.remove(node);
					box_edges_tracker[i].erase(box_edges_tracker[i].begin() + j);
				}
				j--;
				if (box_edges_tracker[i].size() < 2) {
					j++;
					//pps->safety_margin += 0.015;
					if (box_edges_tracker[i][j]->vert->index <= old_max_index) {
						//collisions_passed.push_back(box_edges_tracker[i][j]);
						box_edges_tracker[i][j]->vert->node->EoL = false;
						box_edges_tracker[i][j]->vert->node->on_preserved_edge = false;
						//box_edges_tracker[i][j]->vert->node->lookout_for = true;
						for (int k = 0; k < box_edges_tracker[i][j]->vert->node->adje.size(); k++) {
							if (box_edges_tracker[i][j]->vert->node->adje[k]->preserve) box_edges_tracker[i][j]->vert->node->adje[k]->preserve = false;
						}
						box_edges_tracker[i].erase(box_edges_tracker[i].begin() + j);
						break;
					}
					Vert* vert2 = box_edges_tracker[i][j]->vert;
					Node* node2 = box_edges_tracker[i][j]->vert->node;
					outside_margin_v.push_back(vert2);
					outside_margin_n.push_back(node2);
					mesh.remove(vert2);
					mesh.remove(node2);
					box_edges_tracker[i].erase(box_edges_tracker[i].begin() + j);
					break;
				}
			}
		}
		// This may be needed but was craeting an improper deletion
		//if (box_edges_tracker[i].size() > 2) {
		//	int j = box_edges_tracker[i].size() - 1;
		//	if (edge_not_within_safety_margin(box_edges_tracker[i][j]->vert, box_edges_tracker[i][j - 1]->vert, pps->safety_margin)) {
		//		//pps->safety_margin += 0.008;
		//		if (box_edges_tracker[i][j]->vert->index <= old_max_index) {
		//			//collisions_passed.push_back(box_edges_tracker[i][j]);
		//			box_edges_tracker[i][j]->vert->node->EoL = false;
		//			box_edges_tracker[i][j]->vert->node->on_preserved_edge = false;
		//			box_edges_tracker[i][j]->vert->node->lookout_for = true;
		//			box_edges_tracker[i].erase(box_edges_tracker[i].begin() + j);
		//		}
		//		else {
		//			Vert* vert2 = box_edges_tracker[i][j]->vert;
		//			Node* node2 = box_edges_tracker[i][j]->vert->node;
		//			outside_margin_v.push_back(vert2);
		//			outside_margin_n.push_back(node2);
		//			mesh.remove(vert2);
		//			mesh.remove(node2);
		//			box_edges_tracker[i].erase(box_edges_tracker[i].begin() + j);
		//		}
		//	}
		//}
	}

	mini_reindex_nodes(mesh.nodes);

	for (size_t i = 0; i < outside_margin_v.size(); i++) delete outside_margin_v[i];
	for (size_t i = 0; i < outside_margin_n.size(); i++) delete outside_margin_n[i];

	for (int i = 0; i < box_edges_tracker.size(); i++) {
		// Store in the saved edges vector and mark which edge verts are on for ArcSim purposes
		for (int j = 0; j + 1 < box_edges_tracker[i].size(); j++) {
			//if (box_edges_tracker[i][j]->vert->index == box_edges_tracker[i][j+1]->vert->index) continue;
			saved_edges.push_back(box_edges_tracker[i][j]->vert->index);
			saved_edges.push_back(box_edges_tracker[i][j + 1]->vert->index);
			box_edges_tracker[i][j]->vert->node->which_edge = i; // Which edge its on for later
			if (j + 1 < box_edges_tracker[i].size()) box_edges_tracker[i][j + 1]->vert->node->which_edge = i;
		}
		//box_edges_tracker[i].back()->vert->node->which_edge = i; // which edge its on for later
	}

	if (!pps->points) {
		if (saved_edges.size() == 0) {
			for (int i = 0; i < mesh.edges.size(); i++) {
				mesh.edges[i]->preserve = false;
			}
			pps->PPTimer[1]->toc();
			if (pps->export_timings) {
				pps->PPTimer[0]->export_csv();
				pps->PPTimer[1]->export_csv();
			}
			return false;
		}
	}

	if (saved_edges.size() > 0) {
		pps->once_hack = true;
		// Clear mesh face and edge data
		int before_removal_face_size = mesh.faces.size();
		for (int i = 0; i < before_removal_face_size; i++) {
			Face* cleanup = mesh.faces[0];
			mesh.remove(mesh.faces[0]);
			delete cleanup;
		}
		int before_removal_edge_size = mesh.edges.size();
		for (int i = 0; i < before_removal_edge_size; i++) {
			Edge* cleanup = mesh.edges[0];
			mesh.remove(mesh.edges[0]);
			delete cleanup;
		}

		// FADE2D
		vector<GEOM_FADE2D::Point2> vInputPoints;
		for (int i = 0; i < mesh.verts.size(); i++) {
			vInputPoints.push_back(GEOM_FADE2D::Point2(mesh.verts[i]->u[0], mesh.verts[i]->u[1]));
			vInputPoints.back().setCustomIndex(i * 10); // HACKY?
		}
		GEOM_FADE2D::Fade_2D dt;
		vector<GEOM_FADE2D::Triangle2*> vAllTris;

		vector<GEOM_FADE2D::Segment2> vSegment;
		for (int i = 0; i < saved_edges.size(); i += 2) {
			vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[saved_edges[i]], vInputPoints[saved_edges[i + 1]]));
			// HACKY?
			vInputPoints[saved_edges[i]].setCustomIndex(vInputPoints[saved_edges[i]].getCustomIndex() + 1);
			vInputPoints[saved_edges[i + 1]].setCustomIndex(vInputPoints[saved_edges[i + 1]].getCustomIndex() + 1);
		}

		dt.insert(vInputPoints);

		//dt.show("example3_noConstraints.ps");

		GEOM_FADE2D::ConstraintGraph2* pCG = dt.createConstraint(vSegment, GEOM_FADE2D::CIS_CONSTRAINED_DELAUNAY);
		dt.applyConstraintsAndZones();

		// Refill all the face data
		dt.getTrianglePointers(vAllTris);
		for (std::vector<GEOM_FADE2D::Triangle2*>::iterator it = vAllTris.begin(); it != vAllTris.end(); ++it) {
			GEOM_FADE2D::Triangle2* pT(*it);
			GEOM_FADE2D::Point2* pCorner1(pT->getCorner(0));
			GEOM_FADE2D::Point2* pCorner2(pT->getCorner(1));
			GEOM_FADE2D::Point2* pCorner3(pT->getCorner(2));

			// TODO:: Make sure this is safe
			mesh.add(new Face(mesh.verts[pCorner1->getCustomIndex() / 10], mesh.verts[pCorner2->getCustomIndex() / 10], mesh.verts[pCorner3->getCustomIndex() / 10], Mat3x3(1), Mat3x3(0), material, 0));
			if ((pCorner3->getCustomIndex() % 10 != 0) && (pCorner2->getCustomIndex() % 10 != 0)) {
				if ((mesh.faces.back()->adje[0]->n[0]->which_edge != mesh.faces.back()->adje[0]->n[1]->which_edge) && (!mesh.faces.back()->adje[0]->n[0]->on_corner && !mesh.faces.back()->adje[0]->n[1]->on_corner)) continue;
				mesh.faces.back()->adje[0]->preserve = true;
				mesh.faces.back()->adje[0]->n[0]->on_preserved_edge = true;
				mesh.faces.back()->adje[0]->n[1]->on_preserved_edge = true;
			}
			if ((pCorner2->getCustomIndex() % 10 != 0) && (pCorner1->getCustomIndex() % 10 != 0)) {
				if ((mesh.faces.back()->adje[2]->n[0]->which_edge != mesh.faces.back()->adje[2]->n[1]->which_edge) && (!mesh.faces.back()->adje[2]->n[0]->on_corner && !mesh.faces.back()->adje[2]->n[1]->on_corner)) continue;
				mesh.faces.back()->adje[2]->preserve = true;
				mesh.faces.back()->adje[2]->n[0]->on_preserved_edge = true;
				mesh.faces.back()->adje[2]->n[1]->on_preserved_edge = true;
			}
			if ((pCorner1->getCustomIndex() % 10 != 0) && (pCorner3->getCustomIndex() % 10 != 0)) {
				if ((mesh.faces.back()->adje[1]->n[0]->which_edge != mesh.faces.back()->adje[1]->n[1]->which_edge) && (!mesh.faces.back()->adje[1]->n[0]->on_corner && !mesh.faces.back()->adje[1]->n[1]->on_corner)) continue;
				mesh.faces.back()->adje[1]->preserve = true;
				mesh.faces.back()->adje[1]->n[0]->on_preserved_edge = true;
				mesh.faces.back()->adje[1]->n[1]->on_preserved_edge = true;
			}
		}

		// display remesh for debugging
		if (pps->export_postscript) {
			GEOM_FADE2D::Visualizer2 vis2("debug_constraints.ps");
			dt.show(&vis2, false);

			//vector<GEOM_FADE2D::Point2*> vPointsOfConstraintEdge;
			//pCG->getPolygonVertices(vPointsOfConstraintEdge);
			//for (size_t i = 0; i + 1<vPointsOfConstraintEdge.size(); ++i)
			////for (size_t i = 0; i + 1<2; ++i)
			//{
			//	GEOM_FADE2D::Point2* p0(vPointsOfConstraintEdge[i]);
			//	GEOM_FADE2D::Point2* p1(vPointsOfConstraintEdge[i + 1]);
			//	vis2.addObject(GEOM_FADE2D::Segment2(*p0, *p1), GEOM_FADE2D::Color(1, 0, 0, 0.01));
			//	//vis2.addObject(*p0, GEOM_FADE2D::Color(0, 0, 1, 0.1));
			//	//vis2.addObject(*p1, GEOM_FADE2D::Color(0, 0, 1, 0.1));
			//}
			for (int i = 0; i < saved_edges.size(); i += 2) {
				GEOM_FADE2D::Point2 p0(vInputPoints[saved_edges[i]]);
				GEOM_FADE2D::Point2 p1(vInputPoints[saved_edges[i + 1]]);
				vis2.addObject(GEOM_FADE2D::Segment2(p0, p1), GEOM_FADE2D::Color(1, 0, 0, 0.00001));
			}
			for (int i = 0; i < vInputPoints.size(); ++i)
			{
				std::string text = std::to_string(mesh.nodes[i]->index);
				vis2.addObject(GEOM_FADE2D::Label(vInputPoints[i], text, false, 5), GEOM_FADE2D::Color(0, 0, 1, 0.1));
			}
			vis2.writeFile();
		}
	}

	bool repeat = true;
	int safety_net = 0;
	while (repeat && safety_net < 3) {
		repeat = false;
		safety_net++;
		while (collapse_black_edges(mesh, pps->collapse_non_preserve_thresh, repeat));
		while (collapse_green_edges(mesh, pps->collapse_preserve_thresh, repeat));
		while (remove_skiiny_triangles(mesh, repeat, pps->EOLon));
		//cout << "repeat" << endl;
	}
	while (remove_awkward_faces(mesh));

	//for (int i = 0; i < 3; i++) {
	//	while (collapse_black_edges(mesh, pps->collapse_non_preserve_thresh));
	//	while (collapse_green_edges(mesh, pps->collapse_preserve_thresh));
	//}

	//while (collapse_black_edges(mesh, pps->collapse_non_preserve_thresh));
	//while (collapse_green_edges(mesh, pps->collapse_preserve_thresh));
	//while (remove_awkward_faces());
	//reindex_nodes(mesh.nodes);
	//while (remove_skiiny_triangles(mesh));

	compute_ms_data(mesh);
	//cout << "finished edges" << endl;

	pps->PPTimer[1]->tic();
	if (pps->export_timings) {
		pps->PPTimer[0]->export_csv();
		pps->PPTimer[1]->export_csv();
	}

	return true;

}

struct fnpair {
	int type; // true = corner, false = edge
	int faceid1;
	int faceid2;
	int nodeid;

	explicit fnpair(int t, int nid) :
		type(t), faceid1(-1), faceid2(-1), nodeid(nid) {}

	explicit fnpair(int t, int fid, int nid) :
		type(t), faceid1(fid), faceid2(-1), nodeid(nid) {}

	explicit fnpair(int t, int fid1, int fid2, int nid) :
		type(t), faceid1(fid1), faceid2(fid2), nodeid(nid) {}
};

bool collisionRemeshAlt(Mesh& mesh, vector<shared_ptr<Box> > box, vector<std::shared_ptr<btc::Collision> > collisions, VectorXi EOLS, shared_ptr<pp_settings> pps)
{
	// Assume function is only called if collision is on

	// Box Setup
	MatrixXd verts2(3, mesh.nodes.size());
	MatrixXi faces2(3, mesh.faces.size());
	vector< unique_ptr<fnpair> > fnpairs;

	for (int i = 0; i < mesh.nodes.size(); i++) {
		verts2.col(i) = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
		//if (mesh.nodes[i]->index == 392) cout << mesh.nodes[i]->adje.size() << endl;
		if (mesh.nodes[i]->EoL) {
			if (corner_not_within_safety_margin(mesh.nodes[i]->verts[0]->u, pps->safety_margin, pps->bounds)) {
				mesh.nodes[i]->EoL = false;
				mesh.nodes[i]->EoL_secondary = false;
				mesh.nodes[i]->on_preserved_edge = false;
				mesh.nodes[i]->on_corner = false;
				mesh.nodes[i]->coll_case = -1;
				for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
					mesh.nodes[i]->adje[j]->preserve = false;
				}
				mesh.nodes[i]->preserve_once = true;
				mesh.nodes[i]->preserve = false;
				cout << "EOL removed via safety margin " << i << endl;
			}
			else fnpairs.push_back(make_unique<fnpair>(0, mesh.nodes[i]->index));
		}
		//mesh.nodes[i]->EoL = false;
	}

	for (int i = 0; i < mesh.faces.size(); i++) {
		faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}

	int old_max_index = mesh.nodes.size() - 1;
	int newly_added = 0;
	Material* material = mesh.faces[0]->material; // Cleanup

	vector<int> saved_edges;

	pps->PPTimer[0]->tic();
	for (int b = 0; b < box.size(); b++) {
		collisions.clear();
		if (pps->wire) boxTriCollisionHack2(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, EOLS, false);
		else if (pps->points) {
			pointTriCollision(collisions, box[b]->thresh, pps->verts1, pps->norms1, verts2, faces2, false);
			if (pps->move_points) pps->verts1 += pps->norms1 * pps->pmove1 * pps->h;
		}
		else {
			boxTriCollision(collisions, box[b]->thresh, box[b]->dim, box[b]->E1, verts2, faces2, EOLS, false);
		}
		if (collisions.size() == 0) continue;
		if (pps->matlab_debug_collision) {
			double_to_file(box[b]->thresh, "threshold");
			vec_to_file(box[b]->dim, "whd1");
			mat_to_file(box[b]->E1, "E1");
			mat_to_file(verts2, "verts2");
			VectorXi vvv(3);
			vvv << 1, 1, 1;
			vec_to_file(EOLS, "isEOL2");
			mat_to_file(faces2.colwise() += vvv, "faces2");
		}

		for (int i = 0; i < collisions.size(); i++) {
			if (collisions[i]->count1 == 1 && collisions[i]->count2 == 3) {
				if (mesh.nodes[collisions[i]->verts2(0)]->EoL) {
					continue;
				}
				if (mesh.nodes[collisions[i]->verts2(1)]->EoL) {
					continue;
				}
				if (mesh.nodes[collisions[i]->verts2(2)]->EoL) {
					continue;
				}
				double tx = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[0] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[0] +
					collisions[i]->weights2(2) * mesh.verts[collisions[i]->verts2(2)]->u[0];
				if (tx < 1e-12) tx = 0.0;
				if (pps->bounds[0] - tx < 1e-12) tx = pps->bounds[0];
				double ty = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[1] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[1] +
					collisions[i]->weights2(2) * mesh.verts[collisions[i]->verts2(2)]->u[1];
				if (ty < 1e-12) ty = 0.0;
				if (pps->bounds[1] - ty < 1e-12) ty = pps->bounds[1];
				Vec3 tv(0);
				for (int j = 0; j < 3; j++) {
					tv[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->node->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->node->v[j] +
						collisions[i]->weights2(2) * mesh.verts[collisions[i]->verts2(2)]->node->v[j];
				}
				Vec3 tV(0);
				for (int j = 0; j < 3; j++) {
					tV[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->v[j] +
						collisions[i]->weights2(2) * mesh.verts[collisions[i]->verts2(2)]->v[j];
				}
				Vec3 t(tx, ty, 0);
				if (corner_not_within_safety_margin(t, pps->safety_margin, pps->bounds)) continue; // Don't add collisions within the corner safety margin
																					  //Vec3 n(collisions[i]->pos2(0), collisions[i]->pos2(1), collisions[i]->pos2(2)); // pos2
				Vec3 n(collisions[i]->pos1_(0), collisions[i]->pos1_(1), collisions[i]->pos1_(2)); // pos1_
				mesh.add(new Vert(t, tV));
				mesh.add(new Node(n, n, Vec3(0), 0, 0, false));
				mesh.nodes.back()->v = tv;
				connect(mesh.verts.back(), mesh.nodes.back());

				// Pair of faceid and pointid
				//fnpairs.push_back(make_unique<fnpair>(true, get_enclosing_face(mesh, Vec2(mesh.verts.back()->u[0], mesh.verts.back()->u[1]))->index,mesh.nodes.back()->index));
				fnpairs.push_back(make_unique<fnpair>(1, collisions[i]->tri2, mesh.nodes.back()->index));

				// Don't really need anymore
				mesh.nodes.back()->corner_egdes[0] = collisions[i]->cornerEdges(0) + (12 * b);
				mesh.nodes.back()->corner_egdes[1] = collisions[i]->cornerEdges(1) + (12 * b);
				mesh.nodes.back()->corner_egdes[2] = collisions[i]->cornerEdges(2) + (12 * b);
				mesh.verts.back()->egde_weight[0] = collisions[i]->cornerEdgeWeights(0);
				mesh.verts.back()->egde_weight[1] = collisions[i]->cornerEdgeWeights(1);
				mesh.verts.back()->egde_weight[2] = collisions[i]->cornerEdgeWeights(2);

				// Add collision info to node
				mesh.nodes.back()->has_coll_info = true;
				mesh.nodes.back()->on_corner = true;
				mesh.nodes.back()->preserve = true;
				if (pps->EOLon) mesh.nodes.back()->EoL = true;
				mesh.nodes.back()->coll_norm = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(2));  // FIX
				mesh.nodes.back()->nor_ave = Vec3(collisions[i]->cornerNormal(0), collisions[i]->cornerNormal(1), collisions[i]->cornerNormal(2));

				cout << "EOL corner added " << endl;

				newly_added++;

			}
			else if (collisions[i]->count1 == 2 && collisions[i]->count2 == 2) {
				if (mesh.nodes[collisions[i]->verts2(0)]->EoL) {
					continue;
				}
				if (mesh.nodes[collisions[i]->verts2(1)]->EoL) {
					continue;
				}
				double tx = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[0] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[0];
				if (tx < 1e-12) tx = 0.0;
				if (pps->bounds[0] - tx < 1e-12) tx = pps->bounds[0];
				double ty = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->u[1] +
					collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->u[1];
				if (ty < 1e-12) ty = 0.0;
				if (pps->bounds[1] - ty < 1e-12) ty = pps->bounds[1];
				Vec3 tv(0);
				for (int j = 0; j < 3; j++) {
					tv[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->node->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->node->v[j];
				}
				Vec3 tV(0);
				for (int j = 0; j < 3; j++) {
					tV[j] = collisions[i]->weights2(0) * mesh.verts[collisions[i]->verts2(0)]->v[j] +
						collisions[i]->weights2(1) * mesh.verts[collisions[i]->verts2(1)]->v[j];
				}
				Vec3 t(tx, ty, 0);
				//for (int j = 0; j < mesh.verts.size(); j++) {
				//	if (unsigned_vv_distance(mesh.nodes[j]->verts[0]->u, t) < 1e-6) {
				//		cout << "to close" << endl;
				//		continue;
				//	}
				//}

				if (quick_edge_safety_margin_check(t, pps->safety_margin, pps->bounds)) continue; // Don't add collisions within the corner safety margin
																					 //Vec3 n(collisions[i]->pos2(0), collisions[i]->pos2(1), collisions[i]->pos2(2)); // pos2
				//if (corner_not_within_safety_margin(t, pps->safety_margin, pps->bounds)) continue;

				Vec3 n(collisions[i]->pos1_(0), collisions[i]->pos1_(1), collisions[i]->pos1_(2)); // pos1
				mesh.add(new Vert(t, tV));
				mesh.add(new Node(n, n, Vec3(0), 0, 0, false));
				mesh.nodes.back()->v = tv;
				connect(mesh.verts.back(), mesh.nodes.back());

				// Pair of faceids and nodeid
				fnpairs.push_back(make_unique<fnpair>(2, collisions[i]->adjf[0], collisions[i]->adjf[1], mesh.nodes.back()->index));

				// Add collision info to node
				mesh.nodes.back()->has_coll_info = true;
				mesh.nodes.back()->preserve = true;
				if (pps->EOLon) mesh.nodes.back()->EoL = true;
				mesh.verts.back()->egde_weight[0] = collisions[i]->weights1(0);
				mesh.nodes.back()->which_edge = collisions[i]->edge1 + (12 * b);
				// Old norm
				//mesh.nodes.back()->coll_norm = Vec3(collisions[i]->nor(0), collisions[i]->nor(1), collisions[i]->nor(2));
				// Ave norm
				mesh.nodes.back()->coll_norm = Vec3(collisions[i]->nor1(0), collisions[i]->nor1(1), collisions[i]->nor1(2));

				if(pps->nudge) nudge_check(mesh, pps->collapse_non_preserve_thresh, pps->bounds);

				//cout << "NEW EOL edge point" << endl;

				newly_added++;
			}
		}
	}
	pps->PPTimer[0]->toc();

	pps->PPTimer[1]->tic();
	vector<Vert *> outside_margin_v;
	vector<Node *> outside_margin_n;

	vector<Node *> connected;

	for (int i = 0; i < fnpairs.size(); i++) {
		for (int j = i + 1; j < fnpairs.size(); j++) {
			// 0-0
			//if (fnpairs[i]->type == 0 && fnpairs[j]->type == 0) continue;
			if (fnpairs[i]->type == 0 && fnpairs[j]->type == 0) {
				if (mesh.nodes[fnpairs[i]->nodeid]->on_corner && mesh.nodes[fnpairs[j]->nodeid]->on_corner) continue;
				else if (mesh.nodes[fnpairs[i]->nodeid]->on_corner) {
					if (mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[0] != mesh.nodes[fnpairs[j]->nodeid]->which_edge &&
						mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[1] != mesh.nodes[fnpairs[j]->nodeid]->which_edge &&
						mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[2] != mesh.nodes[fnpairs[j]->nodeid]->which_edge) continue;
				}
				else if (mesh.nodes[fnpairs[j]->nodeid]->on_corner) {
					if (mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[0] != mesh.nodes[fnpairs[i]->nodeid]->which_edge &&
						mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[1] != mesh.nodes[fnpairs[i]->nodeid]->which_edge &&
						mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[2] != mesh.nodes[fnpairs[i]->nodeid]->which_edge) continue;
				}
				else {
					if (mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->which_edge) continue;
				}
				bool match = false;
				for (int k = 0; k < mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf.size(); k++) {
					if (mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf[k]->v[0]->node->index == fnpairs[j]->nodeid ||
						mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf[k]->v[1]->node->index == fnpairs[j]->nodeid || 
						mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf[k]->v[2]->node->index == fnpairs[j]->nodeid) {
						match = true;
						break;
					}
				}
				if (match) {
					// Connect
					connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
					connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
					Edge *e0 = get_edge(mesh.nodes[fnpairs[i]->nodeid], mesh.nodes[fnpairs[j]->nodeid]);
					e0->preserve = true;
				}

			}
			// 0-1
			else if (fnpairs[i]->type == 0 && fnpairs[j]->type == 1) {
				if (mesh.nodes[fnpairs[i]->nodeid]->on_corner) continue;
				if (mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[0] &&
					mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[1] &&
					mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[2]) continue;
				bool match = false;
				for (int k = 0; k < mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf.size(); k++) {
					if (mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf[k]->index == fnpairs[j]->faceid1) {
						match = true;
						break;
					}
				}
				if (match) {
					// Connect
					connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
					connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
				}
			}
			// 0-2
			else if (fnpairs[i]->type == 0 && fnpairs[j]->type == 2) {
				if (mesh.nodes[fnpairs[i]->nodeid]->on_corner) {
					if (mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[0] != mesh.nodes[fnpairs[j]->nodeid]->which_edge &&
						mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[1] != mesh.nodes[fnpairs[j]->nodeid]->which_edge &&
						mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[2] != mesh.nodes[fnpairs[j]->nodeid]->which_edge) continue;
				}
				else if (mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->which_edge) continue;
				bool match = false;
				for (int k = 0; k < mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf.size(); k++) {
					if (mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf[k]->index == fnpairs[j]->faceid1 ||
						mesh.nodes[fnpairs[i]->nodeid]->verts[0]->adjf[k]->index == fnpairs[j]->faceid2) {
						match = true;
						break;
					}
				}
				if (match) {
					// Connect
					connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
					connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
				}
			}
			// 1-0
			else if (fnpairs[i]->type == 1 && fnpairs[j]->type == 0) {
				if (mesh.nodes[fnpairs[j]->nodeid]->on_corner) continue;
				if (mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[0] != mesh.nodes[fnpairs[j]->nodeid]->which_edge &&
					mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[1] != mesh.nodes[fnpairs[j]->nodeid]->which_edge &&
					mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[2] != mesh.nodes[fnpairs[j]->nodeid]->which_edge) continue;
				bool match = false;
				for (int k = 0; k < mesh.nodes[fnpairs[j]->nodeid]->verts[0]->adjf.size(); k++) {
					if (mesh.nodes[fnpairs[j]->nodeid]->verts[0]->adjf[k]->index == fnpairs[i]->faceid1) {
						match = true;
						break;
					}
				}
				if (match) {
					// Connect
					connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
					connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
				}
			}
			// 1-1 Don't ever connect two point EOLs
			else if (fnpairs[i]->type == 1 && fnpairs[j]->type == 1) {
				//if (mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->which_edge) continue;
				//if (fnpairs[i]->faceid1 == fnpairs[j]->faceid1) {
				//	// connect
				//	connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
				//	connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
				//}
			}
			// 1-2
			else if (fnpairs[i]->type == 1 && fnpairs[j]->type == 2) {
				if (mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[0] != mesh.nodes[fnpairs[j]->nodeid]->which_edge &&
					mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[1] != mesh.nodes[fnpairs[j]->nodeid]->which_edge &&
					mesh.nodes[fnpairs[i]->nodeid]->corner_egdes[2] != mesh.nodes[fnpairs[j]->nodeid]->which_edge) continue;
				if (fnpairs[i]->faceid1 == fnpairs[j]->faceid1 ||
					fnpairs[i]->faceid1 == fnpairs[j]->faceid2) {
					// connect
					connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
					connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
				}
			}
			//2-0
			else if (fnpairs[i]->type == 2 && fnpairs[j]->type == 0) {
				if (mesh.nodes[fnpairs[j]->nodeid]->on_corner) {
					if (mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[0] &&
						mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[1] &&
						mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[2]) continue;
				}
				else if (mesh.nodes[fnpairs[j]->nodeid]->which_edge != mesh.nodes[fnpairs[i]->nodeid]->which_edge) continue;
				bool match = false;
				for (int k = 0; k < mesh.nodes[fnpairs[j]->nodeid]->verts[0]->adjf.size(); k++) {
					if (mesh.nodes[fnpairs[j]->nodeid]->verts[0]->adjf[k]->index == fnpairs[i]->faceid1 ||
						mesh.nodes[fnpairs[j]->nodeid]->verts[0]->adjf[k]->index == fnpairs[i]->faceid2) {
						match = true;
						break;
					}
				}
				if (match) {
					// Connect
					connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
					connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
				}
			}
			// 2-1
			else if (fnpairs[i]->type == 2 && fnpairs[j]->type == 1) {
				if (mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[0] &&
					mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[1] &&
					mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->corner_egdes[2]) continue;
				if (fnpairs[i]->faceid1 == fnpairs[j]->faceid1 ||
					fnpairs[i]->faceid2 == fnpairs[j]->faceid1) {
					// connect
					connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
					connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
				}
			}
			// 2-2
			else if (fnpairs[i]->type == 2 && fnpairs[j]->type == 2) {
				if (mesh.nodes[fnpairs[i]->nodeid]->which_edge != mesh.nodes[fnpairs[j]->nodeid]->which_edge) continue;
				if (fnpairs[i]->faceid1 == fnpairs[j]->faceid1 ||
					fnpairs[i]->faceid1 == fnpairs[j]->faceid2 ||
					fnpairs[i]->faceid2 == fnpairs[j]->faceid1 ||
					(fnpairs[i]->faceid2 == fnpairs[j]->faceid2 && fnpairs[i]->faceid2 != -1)) {
					// connect
					connected.push_back(mesh.nodes[fnpairs[i]->nodeid]);
					connected.push_back(mesh.nodes[fnpairs[j]->nodeid]);
				}
			}
		}
	}

	mini_reindex_nodes(mesh.nodes);

	//for (size_t i = 0; i < outside_margin_v.size(); i++) delete outside_margin_v[i];
	//for (size_t i = 0; i < outside_margin_n.size(); i++) delete outside_margin_n[i];

	//if (!pps->points) {
	//	if (connected.size() == 0) {
	//		for (int i = 0; i < mesh.edges.size(); i++) {
	//			mesh.edges[i]->preserve = false;
	//		}
	//		pps->PPTimer[1]->toc();
	//		if (pps->export_timings) {
	//			pps->PPTimer[0]->export_csv();
	//			pps->PPTimer[1]->export_csv();
	//		}
	//		return false;
	//	}
	//}

	if (newly_added > 0) {

		for (int i = 0; i < mesh.edges.size(); i++) {
			if (mesh.edges[i]->preserve) {
				connected.push_back(mesh.edges[i]->n[0]);
				connected.push_back(mesh.edges[i]->n[1]);
			}
		}

		pps->once_hack = true;
		//// Clear mesh face and edge data
		//int before_removal_face_size = mesh.faces.size();
		//for (int i = 0; i < before_removal_face_size; i++) {
		//	Face* cleanup = mesh.faces[0];
		//	mesh.remove(mesh.faces[0]);
		//	delete cleanup;
		//}
		//int before_removal_edge_size = mesh.edges.size();
		//for (int i = 0; i < before_removal_edge_size; i++) {
		//	Edge* cleanup = mesh.edges[0];
		//	mesh.remove(mesh.edges[0]);
		//	delete cleanup;
		//}

		// FADE2D
		vector<GEOM_FADE2D::Point2> vInputPoints;
		for (int i = 0; i < mesh.verts.size(); i++) {
			vInputPoints.push_back(GEOM_FADE2D::Point2(mesh.verts[i]->u[0], mesh.verts[i]->u[1]));
			vInputPoints.back().setCustomIndex(i * 10); // HACKY?
		}
		GEOM_FADE2D::Fade_2D dt;
		vector<GEOM_FADE2D::Triangle2*> vAllTris;

		vector<GEOM_FADE2D::Segment2> vSegment;
		for (int i = 0; i < connected.size(); i += 2) {
			vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[connected[i]->index], vInputPoints[connected[i + 1]->index]));
			// HACKY?
			vInputPoints[connected[i]->index].setCustomIndex(vInputPoints[connected[i]->index].getCustomIndex() + 1);
			vInputPoints[connected[i + 1]->index].setCustomIndex(vInputPoints[connected[i + 1]->index].getCustomIndex() + 1);
		}

		// Fade preserve surrounding edges too
		//for (int i = 0; i < fnpairs.size(); i++) {
		//	switch (fnpairs[i]->type) {
		//	case 0:
		//		for (int j = 0; j < mesh.nodes[fnpairs[i]->nodeid]->adje.size(); j++) {
		//			vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.nodes[fnpairs[i]->nodeid]->adje[j]->n[0]->index], vInputPoints[mesh.nodes[fnpairs[i]->nodeid]->adje[j]->n[1]->index]));
		//		}
		//		break;
		//	case 1:
		//		vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[0]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[0]->n[1]->index]));
		//		vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[1]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[1]->n[1]->index]));
		//		vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[2]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[2]->n[1]->index]));
		//		break;
		//	case 2:
		//		vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[0]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[0]->n[1]->index]));
		//		vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[1]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[1]->n[1]->index]));
		//		vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[2]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid1]->adje[2]->n[1]->index]));
		//		if (fnpairs[i]->faceid2 != -1) {
		//			vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid2]->adje[0]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid2]->adje[0]->n[1]->index]));
		//			vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid2]->adje[1]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid2]->adje[1]->n[1]->index]));
		//			vSegment.push_back(GEOM_FADE2D::Segment2(vInputPoints[mesh.faces[fnpairs[i]->faceid2]->adje[2]->n[0]->index], vInputPoints[mesh.faces[fnpairs[i]->faceid2]->adje[2]->n[1]->index]));
		//		}
		//	}
		//}

		dt.insert(vInputPoints);

		//dt.show("example3_noConstraints.ps");

		// Clear mesh face and edge data
		int before_removal_face_size = mesh.faces.size();
		for (int i = 0; i < before_removal_face_size; i++) {
			Face* cleanup = mesh.faces[0];
			mesh.remove(mesh.faces[0]);
			delete cleanup;
		}
		int before_removal_edge_size = mesh.edges.size();
		for (int i = 0; i < before_removal_edge_size; i++) {
			Edge* cleanup = mesh.edges[0];
			mesh.remove(mesh.edges[0]);
			delete cleanup;
		}

		GEOM_FADE2D::ConstraintGraph2* pCG = dt.createConstraint(vSegment, GEOM_FADE2D::CIS_CONSTRAINED_DELAUNAY);
		dt.applyConstraintsAndZones();

		// Refill all the face data
		dt.getTrianglePointers(vAllTris);
		for (std::vector<GEOM_FADE2D::Triangle2*>::iterator it = vAllTris.begin(); it != vAllTris.end(); ++it) {
			GEOM_FADE2D::Triangle2* pT(*it);
			GEOM_FADE2D::Point2* pCorner1(pT->getCorner(0));
			GEOM_FADE2D::Point2* pCorner2(pT->getCorner(1));
			GEOM_FADE2D::Point2* pCorner3(pT->getCorner(2));

			// TODO:: Make sure this is safe
			mesh.add(new Face(mesh.verts[pCorner1->getCustomIndex() / 10], mesh.verts[pCorner2->getCustomIndex() / 10], mesh.verts[pCorner3->getCustomIndex() / 10], Mat3x3(1), Mat3x3(0), material, 0));
			if ((pCorner3->getCustomIndex() % 10 != 0) && (pCorner2->getCustomIndex() % 10 != 0)) {
				if ((mesh.faces.back()->adje[0]->n[0]->which_edge != mesh.faces.back()->adje[0]->n[1]->which_edge) && (!mesh.faces.back()->adje[0]->n[0]->on_corner && !mesh.faces.back()->adje[0]->n[1]->on_corner)) continue;
				mesh.faces.back()->adje[0]->preserve = true;
				mesh.faces.back()->adje[0]->n[0]->on_preserved_edge = true;
				mesh.faces.back()->adje[0]->n[1]->on_preserved_edge = true;
			}
			if ((pCorner2->getCustomIndex() % 10 != 0) && (pCorner1->getCustomIndex() % 10 != 0)) {
				if ((mesh.faces.back()->adje[2]->n[0]->which_edge != mesh.faces.back()->adje[2]->n[1]->which_edge) && (!mesh.faces.back()->adje[2]->n[0]->on_corner && !mesh.faces.back()->adje[2]->n[1]->on_corner)) continue;
				mesh.faces.back()->adje[2]->preserve = true;
				mesh.faces.back()->adje[2]->n[0]->on_preserved_edge = true;
				mesh.faces.back()->adje[2]->n[1]->on_preserved_edge = true;
			}
			if ((pCorner1->getCustomIndex() % 10 != 0) && (pCorner3->getCustomIndex() % 10 != 0)) {
				if ((mesh.faces.back()->adje[1]->n[0]->which_edge != mesh.faces.back()->adje[1]->n[1]->which_edge) && (!mesh.faces.back()->adje[1]->n[0]->on_corner && !mesh.faces.back()->adje[1]->n[1]->on_corner)) continue;
				mesh.faces.back()->adje[1]->preserve = true;
				mesh.faces.back()->adje[1]->n[0]->on_preserved_edge = true;
				mesh.faces.back()->adje[1]->n[1]->on_preserved_edge = true;
			}
		}

		// display remesh for debugging
		if (pps->export_postscript) {
			GEOM_FADE2D::Visualizer2 vis2("debug_constraints.ps");
			dt.show(&vis2, false);

			//vector<GEOM_FADE2D::Point2*> vPointsOfConstraintEdge;
			//pCG->getPolygonVertices(vPointsOfConstraintEdge);
			//for (size_t i = 0; i + 1<vPointsOfConstraintEdge.size(); ++i)
			////for (size_t i = 0; i + 1<2; ++i)
			//{
			//	GEOM_FADE2D::Point2* p0(vPointsOfConstraintEdge[i]);
			//	GEOM_FADE2D::Point2* p1(vPointsOfConstraintEdge[i + 1]);
			//	vis2.addObject(GEOM_FADE2D::Segment2(*p0, *p1), GEOM_FADE2D::Color(1, 0, 0, 0.01));
			//	//vis2.addObject(*p0, GEOM_FADE2D::Color(0, 0, 1, 0.1));
			//	//vis2.addObject(*p1, GEOM_FADE2D::Color(0, 0, 1, 0.1));
			//}
			for (int i = 0; i < connected.size(); i += 2) {
				GEOM_FADE2D::Point2 p0(vInputPoints[connected[i]->index]);
				GEOM_FADE2D::Point2 p1(vInputPoints[connected[i + 1]->index]);
				vis2.addObject(GEOM_FADE2D::Segment2(p0, p1), GEOM_FADE2D::Color(1, 0, 0, 0.00001));
			}
			for (int i = 0; i < vInputPoints.size(); ++i)
			{
				std::string text = std::to_string(mesh.nodes[i]->index);
				vis2.addObject(GEOM_FADE2D::Label(vInputPoints[i], text, false, 5), GEOM_FADE2D::Color(0, 0, 1, 0.1));
			}
			vis2.writeFile();
		}
	}

	bool repeat = true;
	int safety_net = 0;
	while (repeat && safety_net < 3) {
		repeat = false;
		safety_net++;
		while (collapse_green_edges(mesh, pps->collapse_preserve_thresh, repeat));
		while (collapse_black_edges(mesh, pps->collapse_non_preserve_thresh, repeat));
		while (collapse_black_edges_off_points(mesh, pps->collapse_non_preserve_thresh, repeat));
		//while (collapse_green_edges(mesh, pps->collapse_preserve_thresh, repeat));
		while (remove_skiiny_triangles(mesh, repeat, pps->EOLon));
		//cout << "repeat" << endl;
	}
	while (remove_awkward_faces(mesh));

	// Post collapse preserve checking
	for (int i = 0; i < mesh.nodes.size(); i++) {
		if (mesh.nodes[i]->EoL) {
			for (int j = 0; j < mesh.nodes[i]->adje.size(); j++) {
				Node* n = other_node(mesh.nodes[i]->adje[j], mesh.nodes[i]);
				if (mesh.nodes[i]->on_corner && n->on_corner) {
					continue;
				}
				else if (mesh.nodes[i]->on_corner) {
					if (mesh.nodes[i]->corner_egdes[0] == n->which_edge ||
						mesh.nodes[i]->corner_egdes[1] == n->which_edge ||
						mesh.nodes[i]->corner_egdes[2] == n->which_edge) {
						mesh.nodes[i]->adje[j]->preserve = true;
					}
				}
				else if (n->on_corner) {
					if (n->corner_egdes[0] == mesh.nodes[i]->which_edge ||
						n->corner_egdes[0] == mesh.nodes[i]->which_edge ||
						n->corner_egdes[0] == mesh.nodes[i]->which_edge) {
						mesh.nodes[i]->adje[j]->preserve = true;
					}
				}
				else if (n->EoL && mesh.nodes[i]->which_edge == n->which_edge) {
					mesh.nodes[i]->adje[j]->preserve = true;
				}
			}
		}
	}

	//for (int i = 0; i < 3; i++) {
	//	while (collapse_black_edges(mesh, pps->collapse_non_preserve_thresh));
	//	while (collapse_green_edges(mesh, pps->collapse_preserve_thresh));
	//}

	//while (collapse_black_edges(mesh, pps->collapse_non_preserve_thresh));
	//while (collapse_green_edges(mesh, pps->collapse_preserve_thresh));
	//while (remove_awkward_faces());
	//reindex_nodes(mesh.nodes);
	//while (remove_skiiny_triangles(mesh));

	compute_ms_data(mesh);
	//cout << "finished edges" << endl;

	pps->PPTimer[1]->tic();
	if (pps->export_timings) {
		pps->PPTimer[0]->export_csv();
		pps->PPTimer[1]->export_csv();
	}

	return true;

}