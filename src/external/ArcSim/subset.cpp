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

#include "subset.hpp"
#include "mesh.hpp"
//#include "display.hpp"

using namespace std;

vector<Node*> MeshSubset::get_all_nodes() {
    vector<Node*> nodes;
    nodes.reserve(active_nodes.size() + support_nodes.size());
    nodes.insert(nodes.end(), active_nodes.begin(), active_nodes.end());
    nodes.insert(nodes.end(), support_nodes.begin(), support_nodes.end());
    return nodes;
}

vector<Vert*> MeshSubset::get_verts() {
	vector<Vert*> verts;
	for (size_t i = 0; i<active_nodes.size(); i++)
		for (size_t j = 0; j<active_nodes[i]->verts.size(); j++)
			include(active_nodes[i]->verts[j], verts);
	return verts;
}

vector<Face*> MeshSubset::get_faces() {
	vector<Face*> faces;
	for (size_t i = 0; i<active_nodes.size(); i++) {
        vector<Vert*>& verts = active_nodes[i]->verts;
		for (size_t j = 0; j<verts.size(); j++) {
			vector<Face*>& f = verts[j]->adjf;
			for (size_t k = 0; k<f.size(); k++)
				include(f[k], faces);
		}
    }
	return faces;
}

vector<Edge*> MeshSubset::get_edges() {
	vector<Edge*> edges;
	for (size_t i = 0; i<active_nodes.size(); i++)
		for (size_t j = 0; j<active_nodes[i]->adje.size(); j++)
			include(active_nodes[i]->adje[j], edges);
			
	return edges;
}

void MeshSubset::recompute_support(map<Node*,int>& acc) {
	support_nodes.clear();				
	for (vector<Node*>::iterator i = active_nodes.begin(); i != active_nodes.end(); ++i)
		for (vector<Edge*>::iterator j = (*i)->adje.begin(); j != (*i)->adje.end(); ++j) {
			Node* node = other_node(*j, *i);
			if (!acc[node]) {
				support_nodes.push_back(node);
				acc[node] = 1;
			}
		}
}

void MeshSubset::update_support() {
    map<Node*,int> acc;
    for (size_t i=0; i<active_nodes.size(); i++)
        acc[active_nodes[i]] = 1;

    recompute_support(acc);    
}

void MeshSubset::grow(int rings) {
	map<Node*,int> acc;
	for (size_t i=0; i<active_nodes.size(); i++)
		acc[active_nodes[i]] = 1;

	recompute_support(acc);
	for (int i=0; i<rings; i++) {
		for (size_t j = 0; j<support_nodes.size(); j++)
			active_nodes.push_back(support_nodes[j]);
		recompute_support(acc);
	}
}

void MeshSubset::set_flag(int flag) {
    for(size_t i=0; i<active_nodes.size(); i++)
        active_nodes[i]->flag |= flag;
}

void MeshSubset::clear_flag(int flag) {
    for(size_t i=0; i<active_nodes.size(); i++)
        active_nodes[i]->flag &= ~flag;
}

//void MeshSubset::debug() {
//	for (size_t i=0; i<support_nodes.size(); i++)
//		Annotation::add(support_nodes[i],Vec3(0,0,1));
//	for (size_t i=0; i<active_nodes.size(); i++)
//		Annotation::add(active_nodes[i],(active_nodes[i]->flag & Node::FlagMayBreak) ? Vec3(0,0.4,0) : Vec3(0.8,0,0.8));
//	wait_key();
//}