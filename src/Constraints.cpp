#include "Constraints.h"

#include "Obstacles.h"
#include "Box.h"
#include "Points.h"
#include "external\ArcSim\mesh.hpp"
#include "external\ArcSim\util.hpp"

#include <iostream>

using namespace std;
using namespace Eigen;

Constraints::Constraints()
{

}

void Constraints::init(const shared_ptr<Obstacles> obs)
{
	int total_size = obs->points->num_points;
	for (int b = 0; b < obs->boxes.size(); b++) {
		total_size += (obs->boxes[b]->num_points + obs->boxes[b]->num_edges);
	}

	constraintTable.resize(total_size);
	obsTable.resize(total_size);

	for (int p = 0; p < obs->points->num_points; p++) {
		constraintTable[p].resize(3);
		constraintTable[p] = obs->points->norms.col(p);
	}

	for (int b = 0; b < obs->boxes.size(); b++) {
		for (int p = 0; p < obs->boxes[b]->num_points; p++) {
			int index = obs->points->num_points + (b* obs->boxes[b]->num_points) + (b* obs->boxes[b]->num_points) + p;
			// TODO:: Corner constraints
		}
		for (int e = 0; e < obs->boxes[b]->num_edges; e++) {
			int index = obs->points->num_points + (b* obs->boxes[b]->num_points) + (b* obs->boxes[b]->num_points) + (obs->boxes[b]->num_points + e);
			constraintTable[index].resize(6);
			constraintTable[index].segment(0, 3) = obs->boxes[b]->faceNorms.col(obs->boxes[b]->edgeFaces(0, e));
			constraintTable[index].segment(3, 3) = obs->boxes[b]->faceNorms.col(obs->boxes[b]->edgeFaces(1, e));
		}
	}
}

void Constraints::updateTable(const shared_ptr<Obstacles> obs)
{
	int total_size = obs->points->num_points;
	for (int b = 0; b < obs->boxes.size(); b++) {
		total_size += (obs->boxes[b]->num_points + obs->boxes[b]->num_edges);
	}

	constraintTable.resize(total_size);
	obsTable.resize(total_size);

	for (int p = 0; p < obs->points->num_points; p++) {
		constraintTable[p].resize(3);
		constraintTable[p] = obs->points->norms.col(p);
	}

	for (int b = 0; b < obs->boxes.size(); b++) {
		for (int p = 0; p < obs->boxes[b]->num_points; p++) {
			int index = obs->points->num_points + (b* obs->boxes[b]->num_points) + (b* obs->boxes[b]->num_edges) + p;
			// TODO:: Corner constraints
			Vector3d corner_nor = Vector3d::Zero();
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 2; j++) {
					corner_nor += obs->boxes[b]->faceNorms.col(obs->boxes[b]->edgeFaces(j, obs->boxes[b]->vertEdges1(i, p)));
				}
			}
			corner_nor /= 6.0;
			constraintTable[index].resize(3);
			constraintTable[index] = corner_nor.normalized();

		}
		for (int e = 0; e < obs->boxes[b]->num_edges; e++) {
			int index = obs->points->num_points + (b* obs->boxes[b]->num_points) + (b* obs->boxes[b]->num_edges) + (obs->boxes[b]->num_points + e);
			constraintTable[index].resize(9);
			constraintTable[index].segment(0, 3) = obs->boxes[b]->faceNorms.col(obs->boxes[b]->edgeFaces(0, e));
			constraintTable[index].segment(3, 3) = obs->boxes[b]->faceNorms.col(obs->boxes[b]->edgeFaces(1, e));
			constraintTable[index].segment(6, 3) = obs->boxes[b]->faceNorms.col(obs->boxes[b]->edgeTan(e));
		}
	}
}

typedef Eigen::Triplet<double> T;
vector<T> N_;

void Constraints::fill(const Mesh& mesh, const shared_ptr<Obstacles> obs)
{
	updateTable(obs);

	vector<T> Aeq_;
	vector<T> Aineq_;
	vector<T> beq_;
	vector<T> bineq_;

	int eqsize = 0;
	int ineqsize = 0;

	for (int n = 0; n < mesh.nodes.size(); n++) {
		if (mesh.nodes[n]->EoL) {
			if (mesh.nodes[n]->cornerID >= 0) {
				Vector3d nor = constraintTable[mesh.nodes[n]->cornerID];
				Vector3d ortho1 = Vector3d(0.0, -nor(2), nor(1));
				Vector3d ortho2 = (ortho1.cross(nor)).normalized();
				ortho2.normalize();

				Aineq_.push_back(T(ineqsize, n * 3, -nor(0)));
				Aineq_.push_back(T(ineqsize, n * 3 + 1, -nor(1)));
				Aineq_.push_back(T(ineqsize, n * 3 + 2, -nor(2)));
				ineqsize++;

				Aeq_.push_back(T(ineqsize, n * 3, ortho1(0)));
				Aeq_.push_back(T(ineqsize, n * 3 + 1, ortho1(1)));
				Aeq_.push_back(T(ineqsize, n * 3 + 2, ortho1(2)));
				eqsize++;

				Aeq_.push_back(T(ineqsize, n * 3, ortho2(0)));
				Aeq_.push_back(T(ineqsize, n * 3 + 1, ortho2(1)));
				Aeq_.push_back(T(ineqsize, n * 3 + 2, ortho2(2)));
				eqsize++;
			}
			else {
				Node* node = mesh.nodes[n];

				Aineq_.push_back(T(ineqsize, n * 3, -constraintTable[node->cdEdges[0]](0)));
				Aineq_.push_back(T(ineqsize, n * 3 + 1, -constraintTable[node->cdEdges[0]](1)));
				Aineq_.push_back(T(ineqsize, n * 3 + 2, -constraintTable[node->cdEdges[0]](2)));
				ineqsize++;

				Aineq_.push_back(T(ineqsize, n * 3, -constraintTable[node->cdEdges[0]](3)));
				Aineq_.push_back(T(ineqsize, n * 3 + 1, -constraintTable[node->cdEdges[0]](4)));
				Aineq_.push_back(T(ineqsize, n * 3 + 2, -constraintTable[node->cdEdges[0]](5)));
				ineqsize++;

				//Aeq_.push_back(T(ineqsize, n * 3, constraintTable[mesh.nodes[n]->cdEdges[0]](6)));
				//Aeq_.push_back(T(ineqsize, n * 3 + 1, constraintTable[mesh.nodes[n]->cdEdges[0]](7)));
				//Aeq_.push_back(T(ineqsize, n * 3 + 2, constraintTable[mesh.nodes[n]->cdEdges[0]](8)));
				//eqsize++;

				// If a boundary, the Eulerian constraint stops it from moving outside
				if (is_seam_or_boundary(node)) {
					// Is this sufficient enough?
					Edge* edge;
					for (int e = 0; e < node->adje.size(); e++) {
						if (is_seam_or_boundary(node->adje[e])) {
							edge = node->adje[e];
							break;
						}
					}
					Node* opp_node = other_node(edge, node);
					Vector2d orth_border = Vector2d(node->verts[0]->u[1] - opp_node->verts[0]->u[1], -node->verts[0]->u[0] - opp_node->verts[0]->u[0]).normalized(); // This should be orthogonal to the edge connecting the two nodes
					Aeq_.push_back(T(eqsize, mesh.nodes.size() * 3 + node->EoL_index * 2, orth_border(0)));
					Aeq_.push_back(T(eqsize, mesh.nodes.size() * 3 + node->EoL_index * 2 + 1, orth_border(1)));
					eqsize++;
				}
				// If internal, the Eulerian constraint forces tangential motion to realize in the Lagrangian space
				else {
					Vector2d tan_ave = Vector2d::Zero();;
					int tot_conf = 0;
					for (int e = 0; e < node->adje.size(); e++) {
						if (node->adje[e]->preserve) {
							Edge* edge = node->adje[e];
							tan_ave += Vector2d(edge->n[1]->verts[0]->u[0] - edge->n[0]->verts[0]->u[0], edge->n[1]->verts[0]->u[1] - edge->n[0]->verts[0]->u[1]).normalized();
							tot_conf++;
						}
					}
					tan_ave.normalize();
					Aeq_.push_back(T(eqsize, mesh.nodes.size() * 3 + node->EoL_index * 2, tan_ave(0)));
					Aeq_.push_back(T(eqsize, mesh.nodes.size() * 3 + node->EoL_index * 2 + 1, tan_ave(1)));
					eqsize++;
				}
			}
		}
	}
}