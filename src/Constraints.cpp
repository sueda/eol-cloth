#include "Constraints.h"

#include "Obstacles.h"
#include "Box.h"
#include "Points.h"

#include <iostream>

using namespace std;
using namespace Eigen;

Constraints::Constraints()
{

}

void Constraints::init(shared_ptr<Obstacles> obs)
{
	int total_size = obs->points->num_points;
	for (int b = 0; b < obs->boxes.size(); b++) {
		total_size += (obs->boxes[b]->num_points + obs->boxes[b]->num_edges);
	}

	constraintTable.resize(total_size);

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