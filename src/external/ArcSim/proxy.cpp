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

#include "proxy.hpp"
//#include "simulation.hpp"
#include "magic.hpp"
#include "constraint.hpp"
#include "geometry.hpp"
#include "collisionutil.hpp"

using namespace std;

FloorProxy::FloorProxy(Mesh& mesh){
    update(mesh);
}

CollisionProxy* FloorProxy::clone(Mesh& mesh) {
    return new FloorProxy(mesh);
}

Constraint* FloorProxy::constraint(const Node* node) {
    if (!is_free(node) || node->x[1] - center.x[1] > ::magic.repulsion_thickness)
        return 0;

    IneqCon *con = new IneqCon;
    con->nodes[0] = (Node*)node;
    con->nodes[1] = &center;
    con->nodes[2] = 0;
    con->nodes[3] = 0;
    for (int n = 0; n < 4; n++)
        con->free[n] = n == 0;
    con->w[0] = 1;
    con->w[1] = -1;
    con->w[2] = 0;
    con->w[3] = 0;
    
    con->stiff = ::magic.collision_stiffness * node->a;
    con->n = Vec3(0,1,0);
    
    //con->mu = sim.obs_friction;
    return con;
}

void FloorProxy::update(Mesh& mesh) {
    double y = infinity;
    for (size_t n=0; n<mesh.nodes.size(); n++)
        if (mesh.nodes[n]->x[1] < y)
            y = mesh.nodes[n]->x[1];
    center.x = Vec3(0, y, 0);
}