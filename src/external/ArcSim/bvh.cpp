// Originally from the SELF-CCD library by Min Tang and Dinesh Manocha
// (http://gamma.cs.unc.edu/SELFCD/).
// Modified by Rahul Narain.

/*************************************************************************\

Copyright 2010 The University of North Carolina at Chapel Hill.
All Rights Reserved.

Permission to use, copy, modify and distribute this software and its
documentation for educational, research and non-profit purposes, without
fee, and without a written agreement is hereby granted, provided that the
above copyright notice and the following three paragraphs appear in all
copies.

IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.

THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

The authors may be contacted via:

US Mail:             GAMMA Research Group at UNC
Department of Computer Science
Sitterson Hall, CB #3175
University of N. Carolina
Chapel Hill, NC 27599-3175

Phone:               (919)962-1749

EMail:              geom@cs.unc.edu; tang_m@zju.edu.cn


\**************************************************************************/

#include <stdlib.h>
#include <assert.h>

#include "bvh.hpp"
//#include "collision.hpp"
#include "mesh.hpp"
#include <climits>
#include <utility>
#include <map>
using namespace std;

BOX node_box(const Node *node, bool ccd) {
	BOX box;
	box += node->x;
	if (ccd)
		box += node->x0;
	return box;
}

BOX vert_box(const Vert *vert, bool ccd) {
	return node_box(vert->node, ccd);
}

BOX edge_box(const Edge *edge, bool ccd) {
	BOX box;
	box += node_box(edge->n[0], ccd);
	box += node_box(edge->n[1], ccd);
	return box;
}

BOX face_box(const Face *face, bool ccd) {
	BOX box;
	for (int v = 0; v < 3; v++)
		box += vert_box(face->v[v], ccd);
	return box;
}

BOX dilate(const BOX &box, double d) {
	static double sqrt2 = sqrt(2);
	BOX dbox = box;
	for (int i = 0; i < 3; i++) {
		dbox._dist[i] -= d;
		dbox._dist[i + 9] += d;
	}
	for (int i = 0; i < 6; i++) {
		dbox._dist[3 + i] -= sqrt2*d;
		dbox._dist[3 + i + 9] += sqrt2*d;
	}
	return dbox;
}

bool overlap(const BOX &box0, const BOX &box1, float thickness) {
	return box0.overlaps(dilate(box1, thickness));
}

float
DeformBVHTree::refit()
{
	getRoot()->refit(_ccd);

	return 0.f;
}

BOX
DeformBVHTree::box()
{
	return getRoot()->_box;
}

inline vec3f norm(vec3f &p1, vec3f &p2, vec3f &p3)
{
	vec3f s = p2 - p1;
	vec3f t = p3 - p1;
	vec3f n = cross(s, t);
	return n;
}

void
DeformBVHNode::refit(bool ccd)
{
	if (isLeaf()) {
		_box = face_box(getFace(), ccd);
	}
	else {
		getLeftChild()->refit(ccd);
		getRightChild()->refit(ccd);

		_box = getLeftChild()->_box + getRightChild()->_box;
	}
}

bool
DeformBVHNode::find(Face *face)
{
	if (isLeaf())
		return getFace() == face;

	if (getLeftChild()->find(face))
		return true;

	if (getRightChild()->find(face))
		return true;

	return false;
}

inline float middle_xyz(char xyz, const vec3f &p1, const vec3f &p2, const vec3f &p3)
{
	float t0, t1;

	t0 = MIN(p1[xyz], p2[xyz]);
	t0 = MIN(t0, p3[xyz]);
	t1 = MAX(p1[xyz], p2[xyz]);
	t1 = MAX(t1, p3[xyz]);
	return (t0 + t1)*0.5f;
}

class aap {
public:
	char _xyz;
	float _p;

	FORCEINLINE aap(const BOX &total) {
		vec3f center = total.center();
		char xyz = 2;

		if (total.width() >= total.height() && total.width() >= total.depth()) {
			xyz = 0;
		}
		else
			if (total.height() >= total.width() && total.height() >= total.depth()) {
				xyz = 1;
			}

		_xyz = xyz;
		_p = center[xyz];
	}

	FORCEINLINE bool inside(const vec3f &mid) const {
		return mid[_xyz]>_p;
	}
};

extern float middle_xyz(char xyz, const vec3f &p1, const vec3f &p2, const vec3f &p3);

DeformBVHTree::DeformBVHTree(DeformModel &mdl, bool ccd)
{
	_mdl = &mdl;
	_ccd = ccd;

	if (!mdl.verts.empty())
		Construct();
	else
		_root = NULL;
}

void
DeformBVHTree::Construct()
{
	BOX total;
	int count;

	int num_vtx = _mdl->verts.size(),
		num_tri = _mdl->faces.size();

	for (int i = 0; i<num_vtx; i++) {
		total += _mdl->verts[i]->node->x;
		if (_ccd)
			total += _mdl->verts[i]->node->x0;
	}

	count = num_tri;

	std::map<Face*, BOX> tri_boxes;
	std::map<Face*, vec3f> tri_centers;

	aap  pln(total);

	face_buffer = new Face*[count];
	int left_idx = 0, right_idx = count;

	for (int i = 0; i<num_tri; i++) {
		Face* face = _mdl->faces[i];

		vec3f &p1 = _mdl->faces[i]->v[0]->node->x;
		vec3f &p2 = _mdl->faces[i]->v[1]->node->x;
		vec3f &p3 = _mdl->faces[i]->v[2]->node->x;
		vec3f &pp1 = _mdl->faces[i]->v[0]->node->x0;
		vec3f &pp2 = _mdl->faces[i]->v[1]->node->x0;
		vec3f &pp3 = _mdl->faces[i]->v[2]->node->x0;

		if (_ccd) {
			tri_centers[face] = vec3f(
				(middle_xyz(0, p1, p2, p3) + middle_xyz(0, pp1, pp2, pp3))*0.5f,
				(middle_xyz(1, p1, p2, p3) + middle_xyz(1, pp1, pp2, pp3))*0.5f,
				(middle_xyz(2, p1, p2, p3) + middle_xyz(2, pp1, pp2, pp3))*0.5f);
		}
		else {
			tri_centers[face] = vec3f(
				middle_xyz(0, p1, p2, p3),
				middle_xyz(1, p1, p2, p3),
				middle_xyz(2, p1, p2, p3));
		}

		if (pln.inside(tri_centers[face]))
			face_buffer[left_idx++] = _mdl->faces[i];
		else
			face_buffer[--right_idx] = _mdl->faces[i];

		tri_boxes[face] += p1;
		tri_boxes[face] += p2;
		tri_boxes[face] += p3;

		if (_ccd) {
			tri_boxes[face] += pp1;
			tri_boxes[face] += pp2;
			tri_boxes[face] += pp3;
		}
	}

	_root = new DeformBVHNode();
	_root->_box = total;
	//_root->_count = count;

	if (count == 1) {
		_root->_face = _mdl->faces[0];
		_root->_left = _root->_right = NULL;
	}
	else {
		if (left_idx == 0 || left_idx == count)
			left_idx = count / 2;

		_root->_left = new DeformBVHNode(_root, face_buffer, left_idx, tri_boxes, tri_centers);
		_root->_right = new DeformBVHNode(_root, face_buffer + left_idx, count - left_idx, tri_boxes, tri_centers);
	}
}

DeformBVHTree::~DeformBVHTree()
{
	if (!_root)
		return;
	delete _root;
	delete[] face_buffer;
}

//#################################################################
// called by root
DeformBVHNode::DeformBVHNode()
{
	_face = NULL;
	_left = _right = NULL;
	_parent = NULL;
	//_count = 0;
	_active = true;
}

DeformBVHNode::~DeformBVHNode()
{
	if (_left) delete _left;
	if (_right) delete _right;
}

// called by leaf
DeformBVHNode::DeformBVHNode(DeformBVHNode *parent, Face *face, std::map<Face*, BOX>& tri_boxes)
{
	_left = _right = NULL;
	_parent = parent;
	_face = face;
	_box = tri_boxes[face];
	//_count = 1;
	_active = true;
}

// called by nodes
DeformBVHNode::DeformBVHNode(DeformBVHNode *parent, Face **lst, unsigned int lst_num, std::map<Face*, BOX>& tri_boxes, std::map<Face*, vec3f>& tri_centers)
{
	assert(lst_num > 0);
	_left = _right = NULL;
	_parent = parent;
	_face = NULL;
	//_count = lst_num;
	_active = true;

	if (lst_num == 1) {
		_face = lst[0];
		_box = tri_boxes[lst[0]];
	}
	else { // try to split them
		for (unsigned int t = 0; t<lst_num; t++) {
			_box += tri_boxes[lst[t]];
		}

		if (lst_num == 2) { // must split it!
			_left = new DeformBVHNode(this, lst[0], tri_boxes);
			_right = new DeformBVHNode(this, lst[1], tri_boxes);
		}
		else {
			aap pln(_box);
			unsigned int left_idx = 0, right_idx = lst_num - 1;

			for (unsigned int t = 0; t<lst_num; t++) {
				if (pln.inside(tri_centers[lst[left_idx]]))
					left_idx++;
				else {// swap it
					Face *tmp = lst[left_idx];
					lst[left_idx] = lst[right_idx];
					lst[right_idx--] = tmp;
				}
			}

			int hal = lst_num / 2;
			if (left_idx == 0 || left_idx == lst_num)
			{
				_left = new DeformBVHNode(this, lst, hal, tri_boxes, tri_centers);
				_right = new DeformBVHNode(this, lst + hal, lst_num - hal, tri_boxes, tri_centers);

			}
			else {
				_left = new DeformBVHNode(this, lst, left_idx, tri_boxes, tri_centers);
				_right = new DeformBVHNode(this, lst + left_idx, lst_num - left_idx, tri_boxes, tri_centers);
			}

		}
	}
}
