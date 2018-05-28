#include "Collisions.h"
#include "Box.h"
#include "Points.h"

#include <Eigen\Dense>

using namespace std;
using namespace Eigen;

void CD(Mesh& mesh, shared_ptr<Obstacles> obs, std::vector<std::shared_ptr<btc::Collision> > &cls)
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
	btc::pointTriCollision(cls, obs->cdthreshold, obs->points->pxyz, obs->points->norms, verts2, faces2, false);

	int c = cls.size();
	for (int b = 0; b < obs->num_boxes; b++) {
		btc::boxTriCollision(cls, obs->cdthreshold, obs->boxes[b]->dim, obs->boxes[b]->E1, verts2, faces2, EoLs, false);
		// We need to augment the indices of the box geometry by the object number
		// TODO:: Internally?
		for (c; c < cls.size(); c++) {
			if (cls[c]->count1 == 1 && cls[c]->count2 == 3) {
				cls[c]->verts1(0) = obs->points->num_points + (8 * b) + cls[c]->verts1(0);
			}
			for (int e = 0; e < cls[c]->edge1.size(); e++) {
				cls[c]->edge1[e] = obs->points->num_points + (12 * b) + cls[c]->edge1[e];
			}
		}
	}

	
}