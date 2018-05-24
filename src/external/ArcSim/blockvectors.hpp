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

#ifndef BLOCKVECTORS_HPP
#define BLOCKVECTORS_HPP

#include "vectors.hpp"

template <int m, int n, typename T>
Vec<m*n,T> mat_to_vec (const Mat<m,n,T> &A) {
    Vec<m*n,T> a;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            a[i+j*m] = A(i,j);
    return a;
}

template <int m, int n, typename T>
Mat<m,n,T> vec_to_mat (const Vec<m*n,T> &a) {
    Mat<m,n,T> A;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            A(i,j) = a[i+j*m];
    return A;
}

template <int bn, int m, int n, typename T> Mat<m*bn,n*bn,T> blockdiag (const Mat<m,n,T> &A) {
    Mat<m*bn,n*bn,T> B = 0;
    for (int b = 0; b < bn; b++)
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                B(b*m+i, b*n+j) = A(i,j);
    return B;
}

template <int m, int n> Mat<m*n,m*n,double> transpose () {
    Mat<m*n,m*n,double> T = 0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            T(n*i+j, i+j*m) = 1;
    return T;
}

template <int n> Mat<n*(n+1)/2,n*n> symmetrize ();

template <> inline Mat<3,4> symmetrize<2> () {
    Mat<3,4> S = Mat<3,4>(0);
    S(0,0) = 1.f;
    S(1,3) = 1.f;
    S(2,1) = S(2,2) = 1/2.f;
    return S;
}

#endif
