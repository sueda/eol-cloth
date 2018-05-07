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

#include "constraint.hpp"

#include "magic.hpp"
using namespace std;

double EqCon::value (int *sign) {
    if (sign) *sign = 0;
    return dot(n, node->x - x);
}
MeshGrad EqCon::gradient () {MeshGrad grad; grad.push_back(MeshGradV(node,n)); return grad;}
MeshGrad EqCon::project () {return MeshGrad();}
double EqCon::energy (double value) {return stiff*sq(value)/2.;}
double EqCon::energy_grad (double value) {return stiff*value;}
double EqCon::energy_hess (double value) {return stiff;}
MeshGrad EqCon::friction (double dt, MeshHess &jac) {return MeshGrad();}
bool EqCon::contains(Node* _node) { return node == _node; }

double GlueCon::value (int *sign) {
    if (sign) *sign = 0;
    return dot(n, nodes[1]->x - nodes[0]->x);
}
MeshGrad GlueCon::gradient () {
    MeshGrad grad;
    grad.push_back(MeshGradV(nodes[0],-n));
    grad.push_back(MeshGradV(nodes[1],n));
    return grad;
}
MeshGrad GlueCon::project () {return MeshGrad();}
double GlueCon::energy (double value) {return stiff*sq(value)/2.;}
double GlueCon::energy_grad (double value) {return stiff*value;}
double GlueCon::energy_hess (double value) {return stiff;}
MeshGrad GlueCon::friction (double dt, MeshHess &jac) {return MeshGrad();}
bool GlueCon::contains(Node* _node) { return nodes[0] == _node || nodes[1] == _node; }

double IneqCon::value (int *sign) {
    if (sign)
        *sign = 1;
    double d = 0;
    for (int i = 0; i < 4; i++)
        if (nodes[i])
            d += w[i]*dot(n, nodes[i]->x);
    d -= ::magic.repulsion_thickness;
    return d;
}

MeshGrad IneqCon::gradient () {
    MeshGrad grad;
    for (int i = 0; i < 4; i++)
        if (nodes[i])
            grad.push_back(MeshGradV(nodes[i], w[i]*n));
    return grad;
}

MeshGrad IneqCon::project () {
    double d = value() + ::magic.repulsion_thickness - ::magic.projection_thickness;
    if (d >= 0)
        return MeshGrad();
    double inv_mass = 0;
    for (int i = 0; i < 4; i++)
        if (nodes[i] && free[i])
            inv_mass += sq(w[i])/nodes[i]->m;
    MeshGrad dx;
    for (int i = 0; i < 4; i++)
        if (nodes[i] && free[i])
            dx.push_back(MeshGradV(nodes[i], -(w[i]/nodes[i]->m)/inv_mass*n*d));
    return dx;
}

double violation (double value) {return std::max(-value, 0.);}

double IneqCon::energy (double value) {
    double v = violation(value);
    return stiff*v*v*v/::magic.repulsion_thickness/6;
}
double IneqCon::energy_grad (double value) {
    return -stiff*sq(violation(value))/::magic.repulsion_thickness/2;
}
double IneqCon::energy_hess (double value) {
    return stiff*violation(value)/::magic.repulsion_thickness;
}

MeshGrad IneqCon::friction (double dt, MeshHess &jac) {
    if (mu == 0)
        return MeshGrad();
    double fn = abs(energy_grad(value()));
    if (fn == 0)
        return MeshGrad();
    Vec3 v = Vec3(0);
    double inv_mass = 0;
    for (int i = 0; i < 4; i++) {
        if (nodes[i]) {
            v += w[i]*nodes[i]->v;
            if (free[i])
                inv_mass += sq(w[i])/nodes[i]->m;
        }
    }
    Mat3x3 T = Mat3x3(1) - outer(n,n);
    double vt = norm(T*v);
    double f_by_v = min(mu*fn/vt, 1/(dt*inv_mass));
    // double f_by_v = mu*fn/max(vt, 1e-1);
    MeshGrad force;
    for (int i = 0; i < 4; i++) {
        if (nodes[i] && free[i]) {
            force.push_back(MeshGradV(nodes[i], -w[i]*f_by_v*T*v));
            for (int j = 0; j < 4; j++) {
                if (free[j]) {
                    jac.push_back(MeshHessV(nodes[i], nodes[j], -w[i]*w[j]*f_by_v*T));
                }
            }
        }
    }
    return force;
}

bool IneqCon::contains(Node *_node) {
	for (int i=0; i<4; i++) {
		if (w[i] && nodes[i] == _node)
			return true;
	}
	return false;
}