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

#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include "sparse.hpp"
#include "vectors.hpp"

// Problems

struct NLOpt { // nonlinear optimization problem
    // minimize objective
    int nvar;
    virtual void initialize (double *x) const = 0;
    //virtual double objective (const double *x) const = 0;
    //virtual void precompute (const double *x) const {}
    //virtual void gradient (const double *x, double *g) const = 0;
    //virtual bool hessian (const double *x, SpMat<double> &H) const {
    //    return false; // should return true if implemented
    //};
    virtual void finalize (const double *x) const = 0;
};

struct NLConOpt { // nonlinear constrained optimization problem
    // minimize objective s.t. constraints = or <= 0
    int nvar, ncon;
    virtual void initialize (double *x) const = 0;
    virtual void precompute (const double *x) const {}
    virtual double objective (const double *x) const = 0;
    virtual void obj_grad (const double *x, double *grad) const = 0; // set
    virtual double constraint (const double *x, int j, int &sign) const = 0;
    virtual void con_grad (const double *x, int j, double factor,
                           double *grad) const = 0; // add factor*gradient
    virtual void finalize (const double *x) const = 0;
};

// Algorithms

struct OptOptions {
    int _max_iter;
    double _eps_x, _eps_f, _eps_g;
    OptOptions (): _max_iter(100), _eps_x(1e-6), _eps_f(1e-12), _eps_g(1e-6) {}
    // Named parameter idiom
    // http://www.parashift.com/c++-faq-lite/named-parameter-idiom.html
    OptOptions &max_iter (int n) {_max_iter = n; return *this;}
    OptOptions &eps_x (double e) {_eps_x = e; return *this;}
    OptOptions &eps_f (double e) {_eps_f = e; return *this;}
    OptOptions &eps_g (double e) {_eps_g = e; return *this;}
    int max_iter () {return _max_iter;}
    double eps_x () {return _eps_x;}
    double eps_f () {return _eps_f;}
    double eps_g () {return _eps_g;}
};

void l_bfgs_method (const NLOpt &problem,
                    OptOptions opts=OptOptions(),
                    bool verbose=false);

void line_search_newtons_method (const NLOpt &problem,
                                 OptOptions opts=OptOptions(),
                                 bool verbose=false);

void nonlinear_conjugate_gradient_method (const NLOpt &problem,
                                          OptOptions opts=OptOptions(),
                                          bool verbose=false);

void trust_region_method (const NLOpt &problem,
                          OptOptions opts=OptOptions(),
                          bool verbose=false);

void augmented_lagrangian_method (const NLConOpt &problem,
                                  OptOptions opts=OptOptions(),
                                  bool verbose=false);

// convenience functions for when optimization variables are Vec3-valued

inline Vec3 get_subvec (const double *x, int i) {
    return Vec3(x[i*3+0], x[i*3+1], x[i*3+2]);}
inline void set_subvec (double *x, int i, const Vec3 &xi) {
    for (int j = 0; j < 3; j++) x[i*3+j] = xi[j];}
inline void add_subvec (double *x, int i, const Vec3 &xi) {
    for (int j = 0; j < 3; j++) x[i*3+j] += xi[j];}

template <int n> Vec<n> get_subvec (const double *x, int i) {
    Vec<n> v; for (int j = 0; j < n; j++) v[j] = x[i*n+j]; return v;}
template <int n> void set_subvec (double *x, int i, const Vec<n> &xi) {
    for (int j = 0; j < n; j++) x[i*n+j] = xi[j];}
template <int n> void add_subvec (double *x, int i, const Vec<n> &xi) {
    for (int j = 0; j < n; j++) x[i*n+j] += xi[j];}

//inline Mat3x3 get_submat (SpMat<double> &A, int i, int j) {
//    Mat3x3 Aij;
//    for (int ii = 0; ii < 3; ii++) for (int jj = 0; jj < 3; jj++)
//        Aij(ii,jj) = A(i*3+ii, j*3+jj);
//    return Aij;
//}
//inline void set_submat (SpMat<double> &A, int i, int j, const Mat3x3 &Aij) {
//    for (int ii = 0; ii < 3; ii++) for (int jj = 0; jj < 3; jj++)
//        A(i*3+ii, j*3+jj) = Aij(ii,jj);
//}
//inline void add_submat (SpMat<double> &A, int i, int j, const Mat3x3 &Aij) {
//    for (int ii = 0; ii < 3; ii++) for (int jj = 0; jj < 3; jj++)
//        A(i*3+ii, j*3+jj) += Aij(ii,jj);
//}


#endif
