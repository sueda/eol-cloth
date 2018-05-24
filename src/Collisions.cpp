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

	for (int b = 0; b < obs->num_boxes; b++) {
		btc::boxTriCollision(cls, obs->cdthreshold, obs->boxes[b]->dim, obs->boxes[b]->E1, verts2, faces2, EoLs, false);
	}

	btc::pointTriCollision(cls, obs->cdthreshold, obs->points->pxyz, obs->points->norms, verts2, faces2, false);
}