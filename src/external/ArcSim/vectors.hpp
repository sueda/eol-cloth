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

#ifndef VECTORS_HPP
#define VECTORS_HPP

#include "winport.hpp" // aa: windows bindings, etc

#include <cmath>
#include <iostream>
//#include <zlib.h>


// aa: if defined, AVX SIMD doubles will be used for vector math.
//#define _AVX

// aa: force alignment for global new/delete operators (in all files that include this one)
// aa: force alignment for global new/delete operators (in all files that include this one)
#if defined(_WIN32)
#define __align(sz) __declspec(align(sz))
inline void* malloc_align(size_t size, size_t alignment = 32) { return _aligned_malloc(size, alignment); }
inline void aligned_free(void *ptr) { _aligned_free(ptr); }
#else
// This is for Linux, Apple OS will require a separate treatment...
#define __align(sz) __attribute__((aligned(sz)))
#include <stdlib.h>
#include <assert.h>
inline void* malloc_align(size_t size, size_t alignment = 32) {
	void* ptr;
	if (0 != posix_memalign(&ptr, alignment, size)) {
		ptr = NULL;
	}
	assert(ptr != NULL);
	return ptr;
}
inline void aligned_free(void *ptr) { free(ptr); }
#endif

#if defined(_AVX)
inline void* operator new(size_t sz) { return malloc_align(sz); }
inline void* operator new[](size_t sz) { return malloc_align(sz); }
inline void  operator delete(void *ptr) { aligned_free(ptr); }
inline void  operator delete[](void *ptr) { aligned_free(ptr); }
#endif

inline double sq(double x) { return x*x; }

template <bool cond> struct static_assertion_failure;
template <> struct static_assertion_failure<true> { void operator() () {} };
#define static_assert(cond) static_assertion_failure<cond>();

#define tpl template <int n, typename T>
#define VecnT Vec<n,T>

template <int n, typename T = double> class Vec {
private:
#if defined(_AVX)
	__align(32) T c[n % 4 ? 4 * (1 + n / 4) : n];
#else
	T c[n];
#endif
public:
	Vec() { for (int i = 0; i < n; i++) c[i] = 0; }
	explicit Vec(T x) { for (int i = 0; i < n; i++) c[i] = x; }
	explicit Vec(T x, T y) { static_assert(n == 2); c[0] = x; c[1] = y; }
	explicit Vec(T x, T y, T z) { static_assert(n == 3); c[0] = x; c[1] = y; c[2] = z; }
	explicit Vec(T x, T y, T z, T w) { static_assert(n == 4); c[0] = x; c[1] = y; c[2] = z; c[3] = w; }
	explicit Vec(T v0, T v1, T v2, T v3, T v4, T v5) { static_assert(n == 6); c[0] = v0; c[1] = v1; c[2] = v2; c[3] = v3; c[4] = v4; c[5] = v5; }
	T &operator[] (int i) { return c[i]; }
	const T &operator[] (int i) const { return c[i]; }
};
tpl VecnT operator+ (const VecnT &u) { return u; }
tpl VecnT operator+ (const VecnT &u, const VecnT &v) { VecnT w; for (int i = 0; i < n; i++) w[i] = u[i] + v[i]; return w; }
tpl VecnT &operator+= (VecnT &u, const VecnT &v) { return u = u + v; }
tpl VecnT operator- (const VecnT &u) { VecnT v; for (int i = 0; i < n; i++) v[i] = -u[i]; return v; }
tpl VecnT operator- (const VecnT &u, const VecnT &v) { return u + (-v); }
tpl VecnT &operator-= (VecnT &u, const VecnT &v) { return u = u - v; }
tpl VecnT operator* (const T &a, const VecnT &u) { VecnT v; for (int i = 0; i < n; i++) v[i] = a*u[i]; return v; }
tpl VecnT operator* (const VecnT &u, const T &a) { return a*u; }
tpl VecnT &operator*= (VecnT &u, const T &a) { return u = u*a; }
tpl VecnT operator/ (const VecnT &u, const T &a) { return u*(1 / a); }
tpl VecnT &operator/= (VecnT &u, const T &a) { return u = u / a; }
tpl bool operator== (const VecnT &u, const VecnT &v) { for (int i = 0; i<n; ++i) if (u[i] != v[i]) return false; return true; }
tpl bool operator!= (const VecnT &u, const VecnT &v) { return !(u == v); }
tpl T dot(const VecnT &u, const VecnT &v) { T d = 0; for (int i = 0; i < n; i++) d += u[i] * v[i]; return d; }
tpl T norm2(const VecnT &u) { return dot(u, u); }
tpl T norm(const VecnT &u) { return sqrt(norm2(u)); }
tpl VecnT normalize(const VecnT &u) { T m = norm(u); return m == 0 ? VecnT(0) : u / m; }
tpl std::ostream &operator<< (std::ostream &out, const VecnT &u) { out << "("; for (int i = 0; i < n; i++) out << (i == 0 ? "" : ", ") << u[i]; out << ")"; return out; }
template <typename T> Vec<3, T> cross(const Vec<3, T> &u, const Vec<3, T> &v) { Vec<3, T> w; w[0] = u[1] * v[2] - u[2] * v[1]; w[1] = u[2] * v[0] - u[0] * v[2]; w[2] = u[0] * v[1] - u[1] * v[0]; return w; }
template <typename T> T stp(const Vec<3, T> &u, const Vec<3, T> &v, const Vec<3, T> &w) { return dot(u, cross(v, w)); }
template <typename T> bool right_handed(const Vec<3, T> &u, const Vec<3, T> &v, const Vec<3, T> &w) { return stp(u, v, w) >= 0; }
template <int m, int n, typename T> Vec<m, T> project(const VecnT &u) { Vec<m, T> v; for (int i = 0; i < m; i++) v[i] = (i<n) ? u[i] : 0; return v; }
template <typename T> Vec<2, T> perp(const Vec<2, T> &u) { return Vec<2, T>(-u[1], u[0]); }
inline Vec<2> reduce_xy(const Vec<3>& v) { return Vec<2>(v[0], v[1]); }
inline Vec<3> expand_xy(const Vec<2>& v) { return Vec<3>(v[0], v[1], 0); }
//tpl void serializer_vec(gzFile fp, VecnT& v, bool save) { for (int i = 0; i<3; i++) serializer(fp, v[i], save); }
tpl inline bool is_bullshit(const VecnT& v) { for (int i = 0; i<n; i++) { if (v[i] > 1e100 || v[i] < -1e100 || v[i] != v[i]) return true; } return false; }

#if defined(_AVX)
#if !defined(_WIN32)
#include <x86intrin.h>
#endif
template<> inline Vec<3, double> operator+<3, double>(const Vec<3, double> &u, const Vec<3, double> &v) { __m256d r = _mm256_add_pd((__m256d&)u, (__m256d&)v); return (Vec<3, double>&)r; }
template<> inline Vec<3, double> operator-<3, double>(const Vec<3, double> &u, const Vec<3, double> &v) { __m256d r = _mm256_sub_pd((__m256d&)u, (__m256d&)v); return (Vec<3, double>&)r; }
template<> inline Vec<3, double> operator*<3, double>(const double &a, const Vec<3, double> &v) { __m256d r = _mm256_mul_pd(_mm256_set1_pd(a), (__m256d&)v); return (Vec<3, double>&)r; }
template<> inline Vec<3, double> operator/<3, double>(const Vec<3, double> &u, const double &a) { __m256d r = _mm256_div_pd((__m256d&)u, _mm256_set1_pd(a)); return (Vec<3, double>&)r; }

template<> inline Vec<3, double>& operator+=<3, double>(Vec<3, double> &r, const Vec<3, double> &v) { (__m256d&)r = _mm256_add_pd((__m256d&)r, (__m256d&)v); return r; }
template<> inline Vec<3, double>& operator-=<3, double>(Vec<3, double> &r, const Vec<3, double> &v) { (__m256d&)r = _mm256_sub_pd((__m256d&)r, (__m256d&)v); return r; }

// There are no base templates for these functions, we will define inline versions.
inline Vec<3, double>& operator*=(Vec<3, double> &r, const Vec<3, double> &v) { (__m256d&)r = _mm256_mul_pd((__m256d&)r, (__m256d&)v); return r; }
inline Vec<3, double>& operator/=(Vec<3, double> &r, const Vec<3, double> &v) { (__m256d&)r = _mm256_div_pd((__m256d&)r, (__m256d&)v); return r; }
inline Vec<3, double>  operator*(const Vec<3, double> &u, const Vec<3, double> &v) { __m256d r = _mm256_mul_pd((__m256d&)u, (__m256d&)v); return (Vec<3, double>&)r; }
inline Vec<3, double>  operator/(const Vec<3, double> &u, const Vec<3, double> &v) { __m256d r = _mm256_div_pd((__m256d&)u, (__m256d&)v); return (Vec<3, double>&)r; }
#endif

#undef tpl
#undef VecnT

typedef Vec<2> Vec2;
typedef Vec<3> Vec3;

#define tpl template <int m, int n, typename T>
#define MatmnT Mat<m,n,T>
#define MatnmT Mat<n,m,T>
#define MatnnT Mat<n,n,T>
#define VecmT Vec<m,T>
#define VecnT Vec<n,T>

//aa: transposed matrix functionality
template <int m, int n, typename T = double> class MatTransposed;

template <int m, int n, typename T = double> class Mat {
private:
	VecmT c[n];
public:
	Mat() { for (int j = 0; j < n; j++) c[j] = VecmT(0); }
	explicit Mat(T x) { for (int j = 0; j < n; j++) { c[j] = VecmT(0); if (j < m) c[j][j] = x; } }
	explicit Mat(VecmT x, VecmT y) { static_assert(n == 2); c[0] = x; c[1] = y; }
	explicit Mat(VecmT x, VecmT y, VecmT z) { static_assert(n == 3); c[0] = x; c[1] = y; c[2] = z; }
	explicit Mat(VecmT x, VecmT y, VecmT z, VecmT w) { static_assert(n == 4); c[0] = x; c[1] = y; c[2] = z; c[3] = w; }
	//aa:    static Mat rows (VecnT x, VecnT y) {return Mat<n,2,T>(x,y).t();}
	//aa:    static Mat rows (VecnT x, VecnT y, VecnT z) {return Mat<n,3,T>(x,y,z).t();}
	static Mat rows(VecnT x, VecnT y) { Mat<2, n, T> M; for (int i = 0; i < n; i++) { M.col(i)[0] = x[i]; M.col(i)[1] = y[i]; } return M; }
	static Mat rows(VecnT x, VecnT y, VecnT z) { Mat<3, n, T> M; for (int i = 0; i < n; i++) { M.col(i)[0] = x[i]; M.col(i)[1] = y[i]; M.col(i)[2] = z[i]; } return M; }
	static Mat rows(VecnT x, VecnT y, VecnT z, VecnT w) { Mat<4, n, T> M; for (int i = 0; i < n; i++) { M.col(i)[0] = x[i]; M.col(i)[1] = y[i]; M.col(i)[2] = z[i]; M.col(i)[3] = w[i]; } return M; }
	VecnT row(int i) const { VecnT R; for (int col = 0; col < n; ++col) { R[col] = c[col][i]; } return R; }
	void set_row(int i, const VecnT& v) { for (int col = 0; col < n; ++col) c[col][i] = v[col]; }

	inline T &operator() (int i, int j) { return c[j][i]; }
	inline const T &operator() (int i, int j) const { return c[j][i]; }
	inline VecmT &col(int j) { return c[j]; }
	inline const VecmT &col(int j) const { return c[j]; }
	MatnmT t() const { return transpose(*this); }
	// const MatTransposed<m,n,T>& t () const {return reinterpret_cast<const MatTransposed<m,n,T>&>(*this);}
	MatmnT inv() const { return inverse(*this); }
};
tpl MatmnT operator+ (const MatmnT &A) { return A; }
tpl MatmnT operator+ (const MatmnT &A, const MatmnT &B) { MatmnT C; for (int j = 0; j < n; j++) C.col(j) = A.col(j) + B.col(j); return C; }
tpl MatmnT &operator+= (MatmnT &A, const MatmnT &B) { return A = A + B; }
tpl MatmnT operator- (const MatmnT &A) { MatmnT B; for (int j = 0; j < n; j++) B.col(j) = -A.col(j); return B; }
tpl MatmnT operator- (const MatmnT &A, const MatmnT &B) { return A + (-B); }
tpl MatmnT &operator-= (MatmnT &A, const MatmnT &B) { return A = A - B; }
tpl MatmnT operator* (const T &a, const MatmnT &A) { MatmnT B; for (int j = 0; j < n; j++) B.col(j) = a*A.col(j); return B; }
tpl MatmnT operator* (const MatmnT &A, const T &a) { return a*A; }
tpl MatmnT &operator*= (MatmnT &A, const T &a) { return A = A*a; }
tpl MatmnT operator/ (const MatmnT &A, const T &a) { return A*(1 / a); }
tpl MatmnT &operator/= (MatmnT &A, const T &a) { return A = A / a; }
tpl VecmT operator* (const MatmnT &A, const VecnT &u) { VecmT v = VecmT(0); for (int j = 0; j < n; j++) v += A.col(j)*u[j]; return v; }
template <int m, int n, int o, typename T> Mat<m, o, T> operator* (const Mat<m, n, T> &A, const Mat<n, o, T> &B) { Mat<m, o, T> C; for (int k = 0; k < o; k++) C.col(k) = A*B.col(k); return C; }
tpl MatmnT *operator*= (const MatmnT &A, const MatnnT &B) { return A = A*B; }
tpl MatnmT transpose(const MatmnT &A) { MatnmT B; for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) B(j, i) = A(i, j); return B; }
template <int n, typename T> VecnT diag(const MatnnT &A) { VecnT u; for (int j = 0; j < n; j++) u[j] = A(j, j); return u; }
template <int n, typename T> T trace(const MatnnT &A) { T t = 0; for (int j = 0; j < n; j++) t += A(j, j); return t; }
template <typename T> T det(const Mat<2, 2, T> &A) { return A(0, 0)*A(1, 1) - A(0, 1)*A(1, 0); }
template <typename T> T det(const Mat<3, 3, T> &A) { return stp(A.col(0), A.col(1), A.col(2)); }
template <typename T> Mat<2, 2, T> inverse(const Mat<2, 2, T> &A) { return Mat<2, 2, T>(Vec<2, T>(A(1, 1), -A(1, 0)), Vec<2, T>(-A(0, 1), A(0, 0))) / det(A); }
template <typename T> T wedge(const Vec<2, T> &u, const Vec<2, T> &v) { return u[0] * v[1] - u[1] * v[0]; }
template <typename T> Mat<3, 3, T> inverse(const Mat<3, 3, T> &A) { return Mat<3, 3, T>(cross(A.col(1), A.col(2)), cross(A.col(2), A.col(0)), cross(A.col(0), A.col(1))).t() / det(A); }
template <int n, typename T> MatnnT diag(const VecnT &u) { MatnnT A = MatnnT(0); for (int j = 0; j < n; j++) A(j, j) = u[j]; return A; }
tpl MatmnT outer(const VecmT &u, const VecnT &v) { MatmnT A; for (int j = 0; j < n; j++) A.col(j) = u*v[j]; return A; }
tpl T inner(const MatmnT &a, const MatmnT& b) { T r = 0; for (int j = 0; j<n; j++) for (int i = 0; i<m; i++) r += a.col(j)[i] * b.col(j)[i]; return r; }
tpl std::ostream &operator<< (std::ostream &out, const MatmnT &A) { MatnmT At = transpose(A); out << "(" << std::endl; for (int i = 0; i < m; i++) out << "    " << At.col(i) << (i + 1 == m ? "" : ",") << std::endl; out << ")"; return out; }
inline Mat<2, 2> reduce_xy(const Mat<3, 3>& M) { return Mat<2, 2>(Vec2(M(0, 0), M(0, 1)), Vec2(M(1, 0), M(1, 1))); }
inline Mat<3, 3> expand_xy(const Mat<2, 2>& M) { return Mat<3, 3>(Vec3(M(0, 0), M(0, 1), 0), Vec3(M(1, 0), M(1, 1), 0), Vec3(0, 0, 0)); }
tpl MatmnT max(const MatmnT& a, const MatmnT& b) { MatmnT c; for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) c(i, j) = std::max(a(i, j), b(i, j)); return c; }

// Frobenius norm
tpl T norm2_F(const MatmnT &A) { T a = 0; for (int j = 0; j < n; j++) a += norm2(A.col(j)); return a; }
tpl T norm_F(const MatmnT &A) { return sqrt(norm2_F(A)); }

template <int m1, int n1, int m2, int n2, typename T> Mat<m1, n1, T> project(const Mat<m2, n2, T> &A) { Mat<m1, n1, T> B; for (int j = 0; j < n1; j++) B.col(j) = (j<n2) ? project<m1>(A.col(j)) : Vec<m1, T>(0); return B; }

#undef tpl
#undef MatmnT
#undef MatnnT
#undef VecmT
#undef VecnT

typedef Mat<2, 2> Mat2x2;
typedef Mat<3, 3> Mat3x3;
typedef Mat<3, 2> Mat3x2;
typedef Mat<2, 3> Mat2x3;

template <int n> struct Eig {
	Mat<n, n> Q;
	Vec<n> l;
};

template <int n> Vec<n> eigen_values(const Mat<n, n>& A);
template <int n> Eig<n> eigen_decomposition(const Mat<n, n> &A);

template <int m, int n> struct SVD {
	Mat<m, m> U;
	Vec<n> s;
	Mat<n, n> Vt;
};

template <int m, int n> Vec<n> solve_llsq(const Mat<m, n> &A, const Vec<m>& b);
template <int n> Vec<n> solve_symmetric(const Mat<n, n>& A, const Vec<n>& b);
template <int n> Mat<n, n> get_positive(const Mat<n, n> &A);

template <int m, int n> SVD<m, n> singular_value_decomposition(const Mat<m, n> &A);
template<> SVD<3, 2> singular_value_decomposition<3, 2>(const Mat<3, 2> &A);


template <int m, int n, int o, typename T> Mat<m, o, T> operator* (const Mat<m, n, T> &A, const MatTransposed<o, n, T> &B);
template <int m, int n, typename T> Vec<n, T> operator* (const MatTransposed<m, n, T> &A, const Vec<m, T> &u);
template <int m, int n, int o, typename T> Mat<m, o, T> operator* (const MatTransposed<n, m, T> &A, const Mat<n, o, T> &B);

//aa: transposed matrix functionality
template <int m, int n, typename T> class MatTransposed : protected Mat<m, n, T> {
	friend Vec<n, T> operator*<> (const MatTransposed<m, n, T> &A, const Vec<m, T> &u);
	template <int m1, int n1, int o, typename T1>  friend Mat<m1, o, T1> operator* (const Mat<m1, n1, T1> &A, const MatTransposed<o, n1, T1> &B);
	template <int m1, int n1, int o, typename T1> friend Mat<m1, o, T1> operator* (const MatTransposed<n1, m1, T1> &A, const Mat<n1, o, T1> &B);
public:

	const Mat<m, n, T>& t() const { return static_cast<const Mat<m, n, T>&>(*this); }
};

template <int m, int n, typename T> Vec<n, T> operator* (const MatTransposed<m, n, T> &A, const Vec<m, T> &u)
{
	Vec<n, T> v;
	for (int j = 0; j < n; j++) v[j] = dot(A.col(j), u);
	return v;
}

template <int m, int n, int o, typename T> Mat<m, o, T> operator* (const Mat<m, n, T> &A, const MatTransposed<o, n, T> &B)
{
	Mat<m, o, T> C;
	for (int k = 0; k < o; k++) {
		//C.col(k) = 0; //!!!
		C.col(k) = A.col(0)*B.col(0)[k];
		for (int i = 1; i < n; i++) {
			C.col(k) += A.col(i)*B.col(i)[k];
		}
	}
	return C;
}

template <int m, int n, int o, typename T> Mat<m, o, T> operator* (const MatTransposed<n, m, T> &A, const Mat<n, o, T> &B)
{
	Mat<m, o, T> C;
	for (int k = 0; k < o; k++) {
		for (int j = 0; j < m; j++)
			C.col(k)[j] = dot(A.col(j), B.col(k)); //!!!
	}
	return C;
}

Eig<3> eigen3(const Mat3x3 &B);


/*
template <int m, int n, int o, typename T> Mat<m,o,T> operator* (const MatTransposed<n,m,T> &A, const MatTransposed<o,n,T> &B)
{
Mat<m,o,T> C;
for (int k = 0; k < o; k++) {
for(int j = 0; j < m; j++) {
C.col(k)[j] = 0; //!!!
for(int i = 0; i < n; i++)
C.col(k)[j] += A.col(j)[i]*B.col(i)[k];
}
}
return C;
}
*/
#undef static_assert

#endif
