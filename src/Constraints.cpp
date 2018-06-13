#include "Constraints.h"

#include "Obstacles.h"
#include "Box.h"
#include "Points.h"
#include "external\ArcSim\mesh.hpp"

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
				Aineq_.push_back(T(ineqsize, n * 3, -constraintTable[mesh.nodes[n]->cdEdges[0]](0)));
				Aineq_.push_back(T(ineqsize, n * 3 + 1, -constraintTable[mesh.nodes[n]->cdEdges[0]](1)));
				Aineq_.push_back(T(ineqsize, n * 3 + 2, -constraintTable[mesh.nodes[n]->cdEdges[0]](2)));
				ineqsize++;

				Aineq_.push_back(T(ineqsize, n * 3, -constraintTable[mesh.nodes[n]->cdEdges[0]](3)));
				Aineq_.push_back(T(ineqsize, n * 3 + 1, -constraintTable[mesh.nodes[n]->cdEdges[0]](4)));
				Aineq_.push_back(T(ineqsize, n * 3 + 2, -constraintTable[mesh.nodes[n]->cdEdges[0]](5)));
				ineqsize++;

				//Aeq_.push_back(T(ineqsize, n * 3, constraintTable[mesh.nodes[n]->cdEdges[0]](6)));
				//Aeq_.push_back(T(ineqsize, n * 3 + 1, constraintTable[mesh.nodes[n]->cdEdges[0]](7)));
				//Aeq_.push_back(T(ineqsize, n * 3 + 2, constraintTable[mesh.nodes[n]->cdEdges[0]](8)));
				//eqsize++;


			}
		}
	}
}