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

#include "vectors.hpp"
#include "blockvectors.hpp"
#include "util.hpp"
#include <cmath>
#include <cfloat>
using namespace std;

// LAPACK stuff
extern "C" {

#ifdef _WIN32
#define LAPACKE_dgesvd dgesvd_
#define LAPACKE_dsyev  dsyev_
#endif

#define lapack_int int
#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102
	lapack_int LAPACKE_dsyev(int matrix_order, char jobz, char uplo, lapack_int n,
		double* a, lapack_int lda, double* w);
	lapack_int LAPACKE_dgesvd(int matrix_order, char jobu, char jobvt,
		lapack_int m, lapack_int n, double* a,
		lapack_int lda, double* s, double* u, lapack_int ldu,
		double* vt, lapack_int ldvt, double* superb);

}

template <int n> Vec<n> eigen_values(const Mat<n, n> &A) {
	Vec<n*n> a = mat_to_vec(A);
	Vec<n> w;
	int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'U', n, &a[0], n, &w[0]);
	if (info != 0)
		cout << "LAPACKE_dsyev failed with return value " << info << " on matrix " << A << endl;
	for (int i = 0; i < n / 2; i++) {
		swap(w[i], w[n - i - 1]);
	}
	return w;
}

template <int n> Eig<n> eigen_decomposition(const Mat<n, n> &A) {
	Eig<n> eig;
	Vec<n*n> a = mat_to_vec(A);
	Vec<n> &w = eig.l;
	int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, &a[0], n, &w[0]);
	if (info != 0)
		cout << "LAPACKE_dsyev failed with return value " << info << " on matrix " << A << endl;
	// SSYEV overwrites a with the eigenvectors
	eig.Q = vec_to_mat<n, n>(a);
	for (int i = 0; i < n / 2; i++) {
		swap(eig.l[i], eig.l[n - i - 1]);
		swap(eig.Q.col(i), eig.Q.col(n - i - 1));
	}
	return eig;
}

int dsyevc3(const Mat3x3& A, Vec3& w) {
	double m, c1, c0;

	// Determine coefficients of characteristic poynomial. We write
	//       | a   d   f  |
	//  A =  | d*  b   e  |
	//       | f*  e*  c  |
	double de = A(0, 1) * A(1, 2);                                    // d * e
	double dd = sqr(A(0, 1));                                         // d^2
	double ee = sqr(A(1, 2));                                         // e^2
	double ff = sqr(A(0, 2));                                         // f^2
	m = A(0, 0) + A(1, 1) + A(2, 2);
	c1 = (A(0, 0)*A(1, 1) + A(0, 0)*A(2, 2) + A(1, 1)*A(2, 2))        // a*b + a*c + b*c - d^2 - e^2 - f^2
		- (dd + ee + ff);
	c0 = A(2, 2)*dd + A(0, 0)*ee + A(1, 1)*ff - A(0, 0)*A(1, 1)*A(2, 2)
		- 2.0 * A(0, 2)*de;                                     // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

	double p, sqrt_p, q, c, s, phi;
	p = sqr(m) - 3.0*c1;
	q = m*(p - (3.0 / 2.0)*c1) - (27.0 / 2.0)*c0;
	sqrt_p = sqrt(fabs(p));

	phi = 27.0 * (0.25*sqr(c1)*(p - c1) + c0*(q + 27.0 / 4.0*c0));
	phi = (1.0 / 3.0) * atan2(sqrt(fabs(phi)), q);

	c = sqrt_p*cos(phi);
	s = (1.0 / 1.73205080756887729352744634151)*sqrt_p*sin(phi);

	w[1] = (1.0 / 3.0)*(m - c);
	w[2] = w[1] + s;
	w[0] = w[1] + c;
	w[1] -= s;

	return 0;
}

template<> Vec3 eigen_values<3>(const Mat3x3& A) {
	Vec3 w;
	dsyevc3(A, w);

	// sort eigenvalues
	if (w[1] > w[0]) swap(w[0], w[1]);
	if (w[2] > w[0]) swap(w[0], w[2]);
	if (w[2] > w[1]) swap(w[1], w[2]);
	return w;
}

// http://www.mpi-hd.mpg.de/personalhomes/globes/3x3
template<> Eig<3> eigen_decomposition<3>(const Mat3x3 &B) {
	Eig<3> e;
	Mat3x3& Q = e.Q;
	Mat3x3 A = B;

	double norm;          // Squared norm or inverse norm of current eigenvector
	double n0, n1;        // Norm of first and second columns of A
	double n0tmp, n1tmp;  // "Templates" for the calculation of n0/n1 - saves a few FLOPS
	double thresh;        // Small number used as threshold for floating point comparisons
	double error;         // Estimated maximum roundoff error in some steps
	double wmax;          // The eigenvalue of maximum modulus
	double f, t;          // Intermediate storage
	int i, j;             // Loop counters

						  // Calculate eigenvalues
	dsyevc3(A, e.l);

	wmax = fabs(e.l[0]);
	if ((t = fabs(e.l[1])) > wmax)
		wmax = t;
	if ((t = fabs(e.l[2])) > wmax)
		wmax = t;
	thresh = sqr(8.0 * DBL_EPSILON * wmax);

	// Prepare calculation of eigenvectors
	n0tmp = sqr(A(0, 1)) + sqr(A(0, 2));
	n1tmp = sqr(A(0, 1)) + sqr(A(1, 2));
	Q(0, 1) = A(0, 1)*A(1, 2) - A(0, 2)*A(1, 1);
	Q(1, 1) = A(0, 2)*A(0, 1) - A(1, 2)*A(0, 0);
	Q(2, 1) = sqr(A(0, 1));

	// Calculate first eigenvector by the formula
	//   v[0] = (A - e.l[0]).e1 x (A - e.l[0]).e2
	A(0, 0) -= e.l[0];
	A(1, 1) -= e.l[0];
	Q(0, 0) = Q(0, 1) + A(0, 2)*e.l[0];
	Q(1, 0) = Q(1, 1) + A(1, 2)*e.l[0];
	Q(2, 0) = A(0, 0)*A(1, 1) - Q(2, 1);
	norm = sqr(Q(0, 0)) + sqr(Q(1, 0)) + sqr(Q(2, 0));
	n0 = n0tmp + sqr(A(0, 0));
	n1 = n1tmp + sqr(A(1, 1));
	error = n0 * n1;

	if (n0 <= thresh)         // If the first column is zero, then (1,0,0) is an eigenvector
	{
		Q(0, 0) = 1.0;
		Q(1, 0) = 0.0;
		Q(2, 0) = 0.0;
	}
	else if (n1 <= thresh)    // If the second column is zero, then (0,1,0) is an eigenvector
	{
		Q(0, 0) = 0.0;
		Q(1, 0) = 1.0;
		Q(2, 0) = 0.0;
	}
	else if (norm < sqr(64.0 * DBL_EPSILON) * error)
	{                         // If angle between A[0] and A[1] is too small, don't use
		t = sqr(A(0, 1));       // cross product, but calculate v ~ (1, -A0/A1, 0)
		f = -A(0, 0) / A(0, 1);
		if (sqr(A(1, 1)) > t)
		{
			t = sqr(A(1, 1));
			f = -A(0, 1) / A(1, 1);
		}
		if (sqr(A(1, 2)) > t)
			f = -A(0, 2) / A(1, 2);
		norm = 1.0 / sqrt(1 + sqr(f));
		Q(0, 0) = norm;
		Q(1, 0) = f * norm;
		Q(2, 0) = 0.0;
	}
	else                      // This is the standard branch
	{
		norm = sqrt(1.0 / norm);
		for (j = 0; j < 3; j++)
			Q(j, 0) = Q(j, 0) * norm;
	}


	// Prepare calculation of second eigenvector
	t = e.l[0] - e.l[1];
	if (fabs(t) > 8.0 * DBL_EPSILON * wmax)
	{
		// For non-degenerate eigenvalue, calculate second eigenvector by the formula
		//   v[1] = (A - e.l[1]).e1 x (A - e.l[1]).e2
		A(0, 0) += t;
		A(1, 1) += t;
		Q(0, 1) = Q(0, 1) + A(0, 2)*e.l[1];
		Q(1, 1) = Q(1, 1) + A(1, 2)*e.l[1];
		Q(2, 1) = A(0, 0)*A(1, 1) - Q(2, 1);
		norm = sqr(Q(0, 1)) + sqr(Q(1, 1)) + sqr(Q(2, 1));
		n0 = n0tmp + sqr(A(0, 0));
		n1 = n1tmp + sqr(A(1, 1));
		error = n0 * n1;

		if (n0 <= thresh)       // If the first column is zero, then (1,0,0) is an eigenvector
		{
			Q(0, 1) = 1.0;
			Q(1, 1) = 0.0;
			Q(2, 1) = 0.0;
		}
		else if (n1 <= thresh)  // If the second column is zero, then (0,1,0) is an eigenvector
		{
			Q(0, 1) = 0.0;
			Q(1, 1) = 1.0;
			Q(2, 1) = 0.0;
		}
		else if (norm < sqr(64.0 * DBL_EPSILON) * error)
		{                       // If angle between A[0] and A[1] is too small, don't use
			t = sqr(A(0, 1));     // cross product, but calculate v ~ (1, -A0/A1, 0)
			f = -A(0, 0) / A(0, 1);
			if (sqr(A(1, 1)) > t)
			{
				t = sqr(A(1, 1));
				f = -A(0, 1) / A(1, 1);
			}
			if (sqr(A(1, 2)) > t)
				f = -A(0, 2) / A(1, 2);
			norm = 1.0 / sqrt(1 + sqr(f));
			Q(0, 1) = norm;
			Q(1, 1) = f * norm;
			Q(2, 1) = 0.0;
		}
		else
		{
			norm = sqrt(1.0 / norm);
			for (j = 0; j < 3; j++)
				Q(j, 1) = Q(j, 1) * norm;
		}
	}
	else
	{
		// For degenerate eigenvalue, calculate second eigenvector according to
		//   v[1] = v[0] x (A - e.l[1]).e[i]
		//   
		// This would really get to complicated if we could not assume all of A to
		// contain meaningful values.
		A(1, 0) = A(0, 1);
		A(2, 0) = A(0, 2);
		A(2, 1) = A(1, 2);
		A(0, 0) += e.l[0];
		A(1, 1) += e.l[0];
		for (i = 0; i < 3; i++)
		{
			A(i, i) -= e.l[1];
			n0 = sqr(A(0, i)) + sqr(A(1, i)) + sqr(A(2, i));
			if (n0 > thresh)
			{
				Q(0, 1) = Q(1, 0)*A(2, i) - Q(2, 0)*A(1, i);
				Q(1, 1) = Q(2, 0)*A(0, i) - Q(0, 0)*A(2, i);
				Q(2, 1) = Q(0, 0)*A(1, i) - Q(1, 0)*A(0, i);
				norm = sqr(Q(0, 1)) + sqr(Q(1, 1)) + sqr(Q(2, 1));
				if (norm > sqr(256.0 * DBL_EPSILON) * n0) // Accept cross product only if the angle between
				{                                         // the two vectors was not too small
					norm = sqrt(1.0 / norm);
					for (j = 0; j < 3; j++)
						Q(j, 1) = Q(j, 1) * norm;
					break;
				}
			}
		}

		if (i == 3)    // This means that any vector orthogonal to v[0] is an EV.
		{
			for (j = 0; j < 3; j++)
				if (Q(j, 0) != 0.0)                                   // Find nonzero element of v[0] ...
				{                                                     // ... and swap it with the next one
					norm = 1.0 / sqrt(sqr(Q(j, 0)) + sqr(Q((j + 1) % 3, 0)));
					Q(j, 1) = Q((j + 1) % 3, 0) * norm;
					Q((j + 1) % 3, 1) = -Q(j, 0) * norm;
					Q((j + 2) % 3, 1) = 0.0;
					break;
				}
		}
	}

	// Calculate third eigenvector according to
	//   v[2] = v[0] x v[1]
	Q(0, 2) = Q(1, 0)*Q(2, 1) - Q(2, 0)*Q(1, 1);
	Q(1, 2) = Q(2, 0)*Q(0, 1) - Q(0, 0)*Q(2, 1);
	Q(2, 2) = Q(0, 0)*Q(1, 1) - Q(1, 0)*Q(0, 1);

	// sort eigenvectors
	if (e.l[1] > e.l[0]) { swap(e.Q.col(0), e.Q.col(1)); swap(e.l[0], e.l[1]); }
	if (e.l[2] > e.l[0]) { swap(e.Q.col(0), e.Q.col(2)); swap(e.l[0], e.l[2]); }
	if (e.l[2] > e.l[1]) { swap(e.Q.col(1), e.Q.col(2)); swap(e.l[1], e.l[2]); }

	return e;
}

template<> Vec2 eigen_values<2>(const Mat2x2& A) {
	double a = A(0, 0), b = A(1, 0), d = A(1, 1); // A(1,0) == A(0,1)
	double amd = a - d;
	double apd = a + d;
	double b2 = b * b;
	double det = sqrt(4 * b2 + amd*amd);
	double l1 = 0.5 * (apd + det);
	double l2 = 0.5 * (apd - det);
	return Vec2(l1, l2);
}

template<> Eig<2> eigen_decomposition<2>(const Mat2x2 &A) {
	// http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
	// http://en.wikipedia.org/wiki/Eigenvalue_algorithm
	Eig<2> eig;
	double a = A(0, 0), b = A(1, 0), d = A(1, 1); // A(1,0) == A(0,1)
	double amd = a - d;
	double apd = a + d;
	double b2 = b * b;
	double det = sqrt(4 * b2 + amd*amd);
	double l1 = 0.5 * (apd + det);
	double l2 = 0.5 * (apd - det);

	eig.l[0] = l1;
	eig.l[1] = l2;

	double v0, v1, vn;
	if (b) {
		v0 = l1 - d;
		v1 = b;
		vn = sqrt(v0*v0 + b2);
		eig.Q(0, 0) = v0 / vn;
		eig.Q(1, 0) = v1 / vn;

		v0 = l2 - d;
		vn = sqrt(v0*v0 + b2);
		eig.Q(0, 1) = v0 / vn;
		eig.Q(1, 1) = v1 / vn;
	}
	else if (a >= d) {
		eig.Q(0, 0) = 1;
		eig.Q(1, 0) = 0;
		eig.Q(0, 1) = 0;
		eig.Q(1, 1) = 1;
	}
	else {
		eig.Q(0, 0) = 0;
		eig.Q(1, 0) = 1;
		eig.Q(0, 1) = 1;
		eig.Q(1, 1) = 0;
	}

	return eig;
}

//template <int m, int n> SVD<m, n> singular_value_decomposition(const Mat<m, n> &A) {
//	SVD<m, n> svd;
//	Vec<m*n> a = mat_to_vec(A);
//	Vec<m*m> u;
//	Vec<n> &s = svd.s;
//	Vec<n*n> vt;
//	Vec<n> superb;
//	int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'A', 'A', m, n, &a[0], m,
//		&s[0], &u[0], m, &vt[0], n, &superb[0]);
//	if (info != 0)
//		cout << "LAPACKE_dgesvd failed with return value " << info << " on matrix " << A << endl;
//	svd.U = vec_to_mat<m, m>(u);
//	svd.Vt = vec_to_mat<n, n>(vt);
//	return svd;
//}

template SVD<3, 3> singular_value_decomposition<3, 3>(const Mat<3, 3> &);

template<> SVD<3, 2> singular_value_decomposition<3, 2>(const Mat<3, 2> &A) {
	//SVD<3,2> svd0 = singular_value_decomposition0(A);
	SVD<3, 2> svd;
	const Vec<3>& c0 = A.col(0);
	const Vec<3>& c1 = A.col(1);
	double a0 = dot(c0, c0); //        |a0  b|
	double b = dot(c0, c1); // A*A' = |     |
	double a1 = dot(c1, c1); //        |b  a1|
	double am = a0 - a1;
	double ap = a0 + a1;
	double det = sqrt(am*am + 4 * b*b);
	// eigen values
	double ev0 = sqrt(0.5 * (ap + det));
	double ev1 = sqrt(0.5 * (ap - det));
	svd.s[0] = ev0;
	svd.s[1] = ev1;
	// http://en.wikipedia.org/wiki/Trigonometric_identities
	double sina, cosa;
	if (b == 0) {
		sina = 0;
		cosa = 1;
	}
	else {
		double tana = (am - det) / (2 * b);
		cosa = 1.0 / sqrt(1 + tana*tana);
		sina = tana * cosa;
	}
	// 2x2
	svd.Vt(0, 0) = -cosa;
	svd.Vt(1, 0) = sina;
	svd.Vt(0, 1) = sina;
	svd.Vt(1, 1) = cosa;
	// 3x3
	double t00 = -cosa / ev0;
	double t10 = sina / ev0;
	svd.U(0, 0) = t00*c0[0] + t10*c1[0];
	svd.U(1, 0) = t00*c0[1] + t10*c1[1];
	svd.U(2, 0) = t00*c0[2] + t10*c1[2];
	double t01 = sina / ev1;
	double t11 = cosa / ev1;
	svd.U(0, 1) = t01*c0[0] + t11*c1[0];
	svd.U(1, 1) = t01*c0[1] + t11*c1[1];
	svd.U(2, 1) = t01*c0[2] + t11*c1[2];
	svd.U.col(2) = cross(svd.U.col(0), svd.U.col(1));
	return svd;
}

template <int n> Mat<n, n> get_positive(const Mat<n, n> &A) {
	Eig<n> eig = eigen_decomposition(A);
	for (int i = 0; i < n; i++)
		eig.l[i] = max(eig.l[i], 0.);
	return eig.Q*diag(eig.l)*eig.Q.t();
}
template Mat2x2 get_positive<2>(const Mat2x2&);
template Mat3x3 get_positive<3>(const Mat3x3&);

template<> Vec2 solve_symmetric(const Mat2x2& A, const Vec2& b) {
	double div = sq(A(0, 1)) - A(1, 1)*A(0, 0);
	if (fabs(div)< 1e-14) {
		cout << A << endl;
		cout << div << endl;
		cout << "singular matrix" << endl; exit(1);
	}
	return Vec2(b[1] * A(0, 1) - b[0] * A(1, 1), b[0] * A(0, 1) - b[1] * A(0, 0)) / div;
}

template<> Vec3 solve_symmetric(const Mat3x3& A, const Vec3& b) {
	double t13 = A(1, 2)*A(1, 2);
	double t14 = A(0, 2)*A(0, 2);
	double t15 = A(0, 0)*t13;
	double t16 = A(1, 1)*t14;
	double t17 = A(0, 1)*A(0, 1);
	double t18 = A(2, 2)*t17;
	double t21 = A(0, 1)*A(0, 2)*A(1, 2)*2.0;
	double t22 = A(0, 0)*A(1, 1)*A(2, 2);
	double t19 = t15 + t16 + t18 - t21 - t22;
	if (fabs(t19) == 0) {
		cout << A << endl << "singular matrix" << endl; exit(1);
	}
	double t20 = 1.0 / t19;
	return Vec3(t20*(t13*b[0] + A(0, 2)*(A(1, 1)*b[2] - A(1, 2)*b[1]) - A(0, 1)*(A(1, 2)*b[2] - A(2, 2)*b[1]) - A(1, 1)*A(2, 2)*b[0]),
		t20*(t14*b[1] + A(1, 2)*(A(0, 0)*b[2] - A(0, 2)*b[0]) - A(0, 1)*(A(0, 2)*b[2] - A(2, 2)*b[0]) - A(0, 0)*A(2, 2)*b[1]),
		t20*(t17*b[2] + A(1, 2)*(A(0, 0)*b[1] - A(0, 1)*b[0]) - A(0, 2)*(A(0, 1)*b[1] - A(1, 1)*b[0]) - A(0, 0)*A(1, 1)*b[2]));
}

template <int m, int n> Vec<n> solve_llsq(const Mat<m, n> &A, const Vec<m>& b) {
	Mat<n, n> M;
	Vec<n> y;
	for (int i = 0; i<n; i++) {
		y[i] = dot(b, A.col(i));
		for (int j = 0; j<n; j++)
			M(i, j) = dot(A.col(i), A.col(j));
	}
	if (norm_F(M) == 0) {
		cout << "llsq: normF = 0 " << endl;
		exit(1);
	}
	return solve_symmetric(M, y);
}

template Vec3 solve_llsq(const Mat<6, 3>&, const Vec<6>&);