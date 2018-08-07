#include "Collisions.h"
#include "Box.h"
#include "Points.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen\Dense>

using namespace std;
using namespace Eigen;

void CD(const Mesh& mesh, const shared_ptr<Obstacles> obs, std::vector<std::shared_ptr<btc::Collision> > &cls)
{
	MatrixXd verts2(3, mesh.nodes.size());
	MatrixXi faces2(3, mesh.faces.size());
	//VectorXi EoLs(1, mesh.nodes.size());
	VectorXi EoLs;
	EoLs.resize(mesh.nodes.size());



	for (int i = 0; i < mesh.nodes.size(); i++) {
		verts2.col(i) = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
		if (mesh.nodes[i]->EoL) EoLs(i) = 1;
	}
	for (int i = 0; i < mesh.faces.size(); i++) {
		faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}
	
	// Compute these first so they form the base of our collision list
	btc::pointTriCollision(cls, obs->cdthreshold, obs->points->pxyz, obs->points->norms, verts2, faces2, true);

	int c = cls.size();
	for (int b = 0; b < obs->num_boxes; b++) {
		vector<shared_ptr<btc::Collision> > clst;
		btc::boxTriCollision(clst, obs->cdthreshold, obs->boxes[b]->dim, obs->boxes[b]->E1, verts2, faces2, EoLs, false);
		cls.insert(cls.end(), clst.begin(), clst.end());
		// We need to augment the indices of the box geometry by the object number
		// TODO:: Internally?
		for (c; c < cls.size(); c++) {
			if (cls[c]->count1 == 1 && cls[c]->count2 == 3) {
				//cls[c]->verts1(0) = obs->points->num_points + (obs->boxes[b]->num_points * b) + cls[c]->verts1(0);
				cls[c]->verts1(0) = obs->points->num_points + (b* obs->boxes[b]->num_points) + (b* obs->boxes[b]->num_edges) + cls[c]->verts1(0);
			}
			for (int e = 0; e < cls[c]->edge1.size(); e++) {
				//cls[c]->edge1[e] = obs->points->num_points + (obs->boxes[b]->num_edges * b) + cls[c]->edge1[e];
				//cls[c]->edge1[e] = (obs->boxes[b]->num_edges * b) + cls[c]->edge1[e];
				cls[c]->edge1[e] = obs->points->num_points + (b* obs->boxes[b]->num_points) + (b* obs->boxes[b]->num_edges) + (obs->boxes[b]->num_points + cls[c]->edge1[e]);
			}
		}
	}

	
}

void CD2(const Mesh& mesh, const shared_ptr<Obstacles> obs, std::vector<std::shared_ptr<btc::Collision> > &cls)
{
	MatrixXd verts2(3, mesh.nodes.size());
	MatrixXi faces2(3, mesh.faces.size());
	VectorXi EoLs;
	EoLs.resize(mesh.nodes.size());

	for (int i = 0; i < mesh.nodes.size(); i++) {
		verts2.col(i) = Vector3d(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
		if (mesh.nodes[i]->EoL) EoLs(i) = 1;
	}
	for (int i = 0; i < mesh.faces.size(); i++) {
		faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}

	// Compute these first so they form the base of our collision list
	btc::pointTriCollision(cls, obs->cdthreshold, obs->points->pxyz, obs->points->norms, verts2, faces2, false);

	for (int b = 0; b < obs->num_boxes; b++) {
		vector<shared_ptr<btc::Collision> > clst;
		btc::boxTriCollision(clst, obs->cdthreshold, obs->boxes[b]->dim, obs->boxes[b]->E1, verts2, faces2, EoLs, false);
		cls.insert(cls.end(), clst.begin(), clst.end());
	}
}