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

#include "nearobs.hpp"

#include "collisionutil.hpp"
#include "geometry.hpp"
#include "magic.hpp"
//#include "simulation.hpp"
#include <map>
using namespace std;

template <typename T> struct Min {
	double key;
	T val;
	Min() : key(infinity), val() {}
	void add(double key, T val) {
		if (key < this->key) {
			this->key = key;
			this->val = val;
		}
	}
};

Vec3 nearest_point(const Vec3 &x, const vector<AccelStruct*> &accs,
	double dmin);

map<Node*, Plane> nearest_obstacle_planes(const vector<Node*>& nodes,
	const vector<AccelStruct*>& obs_accs) {
	const double dmin = 10 * ::magic.repulsion_thickness;
	vector<Plane> planes(nodes.size(), Plane(Vec3(0), Vec3(0)));
#pragma omp parallel for
	for (size_t n = 0; n < nodes.size(); n++) {
		Vec3 x = nodes[n]->x;
		Vec3 p = nearest_point(x, obs_accs, dmin);
		if (p != x)
			planes[n].x0 = p;
		planes[n].n = normalize(x - p);
	}
	map<Node*, Plane> plane_map;
	for (size_t i = 0; i<planes.size(); i++)
		if (norm2(planes[i].n) != 0)
			plane_map[(Node*)nodes[i]] = planes[i];
	return plane_map;
}

struct NearPoint {
	double d;
	Vec3 x;
	NearPoint(double d, const Vec3 &x) : d(d), x(x) {}
};

void update_nearest_point(const Vec3 &x, BVHNode *node, NearPoint &p);

Vec3 nearest_point(const Vec3 &x, const vector<AccelStruct*> &accs,
	double dmin) {
	NearPoint p(dmin, x);
	for (int a = 0; a < (int)accs.size(); a++)
		if (accs[a]->root)
			update_nearest_point(x, accs[a]->root, p);
	return p.x;
}

void update_nearest_point(const Vec3 &x, const Face *face, NearPoint &p);

double point_box_distance(const Vec3 &x, const BOX &box);

void update_nearest_point(const Vec3 &x, BVHNode *node, NearPoint &p) {
	if (node->isLeaf())
		update_nearest_point(x, node->getFace(), p);
	else {
		double d = point_box_distance(x, node->_box);
		if (d >= p.d)
			return;
		update_nearest_point(x, node->getLeftChild(), p);
		update_nearest_point(x, node->getRightChild(), p);
	}
}

double point_box_distance(const Vec3 &x, const BOX &box) {
	Vec3 xp = Vec3(clamp(x[0], (double)box._dist[0], (double)box._dist[9]),
		clamp(x[1], (double)box._dist[1], (double)box._dist[10]),
		clamp(x[2], (double)box._dist[2], (double)box._dist[11]));
	return norm(x - xp);
}

void update_nearest_point(const Vec3 &x, const Face *face, NearPoint &p) {
	Vec3 n;
	double w[4];
	double d = unsigned_vf_distance(x, face->v[0]->node->x, face->v[1]->node->x,
		face->v[2]->node->x, &n, w);
	if (d < p.d) {
		p.d = d;
		p.x = -(w[1] * face->v[0]->node->x + w[2] * face->v[1]->node->x
			+ w[3] * face->v[2]->node->x);
	}
}
