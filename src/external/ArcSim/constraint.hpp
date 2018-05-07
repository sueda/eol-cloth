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

#ifndef CONSTRAINT_HPP
#define CONSTRAINT_HPP

#include "mesh.hpp"
#include "sparse.hpp"
#include "spline.hpp"
#include "util.hpp"
#include "vectors.hpp"
#include <map>
#include <vector>

struct MeshGradV { Node* node; Vec3 f; MeshGradV(Node* i, const Vec3& f) : node(i), f(f) {} };
struct MeshHessV { Node* i; Node* j; Mat3x3 J; MeshHessV(Node *i, Node* j, const Mat3x3& J) : i(i), j(j), J(J) {}};

typedef std::vector<MeshGradV> MeshGrad;
typedef std::vector<MeshHessV> MeshHess;

struct Constraint {
    virtual ~Constraint () {};
    virtual double value (int *sign=NULL) = 0;
    virtual MeshGrad gradient () = 0;
    virtual MeshGrad project () = 0;
    virtual bool contains(Node* node) = 0;
    // energy function
    virtual double energy (double value) = 0;
    virtual double energy_grad (double value) = 0;
    virtual double energy_hess (double value) = 0;
    // frictional force
    virtual MeshGrad friction (double dt, MeshHess &jac) = 0;

    // useful for debugging
    virtual void serializer(Serialize& s, const std::string& name) {};
};

struct EqCon: public Constraint {
    // n . (node->x - x) = 0
    Node *node;
    Vec3 x, n;
    double stiff;
    double value (int *sign=NULL);
    bool contains(Node *node);
    MeshGrad gradient ();
    MeshGrad project ();
    double energy (double value);
    double energy_grad (double value);
    double energy_hess (double value);
    MeshGrad friction (double dt, MeshHess &jac);
};

struct GlueCon: public Constraint {
    Node *nodes[2];
    Vec3 n;
    double stiff;
    double value (int *sign=NULL);
    bool contains(Node *node);    
    MeshGrad gradient ();
    MeshGrad project ();
    double energy (double value);
    double energy_grad (double value);
    double energy_hess (double value);
    MeshGrad friction (double dt, MeshHess &jac);
};

struct IneqCon: public Constraint {
    // n . sum(w[i] verts[i]->x) >= 0
    Node *nodes[4];
    double w[4];
    bool free[4];
    Vec3 n;
    double a; // area
    double mu; // friction
    double stiff;
    double value (int *sign=NULL);
    bool contains(Node *node);    
    MeshGrad gradient ();
    MeshGrad project ();
    double energy (double value);
    double energy_grad (double value);
    double energy_hess (double value);
    MeshGrad friction (double dt, MeshHess &jac);

};

#endif
