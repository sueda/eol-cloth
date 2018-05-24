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

#ifndef TRANSFORMATION_HPP
#define TRANSFORMATION_HPP

//#include "spline.hpp"
#include "vectors.hpp"
#include <iostream>

// Transform the mesh
struct Quaternion {
    double s;
    Vec3 v;
    Vec3 rotate (const Vec3 &point) const;
    static Quaternion from_axisangle(const Vec3 &axis, double angle);
    std::pair<Vec3, double> to_axisangle() const;
    Quaternion operator+(const Quaternion& q) const;
    Quaternion operator-(const Quaternion& q) const;
    Quaternion operator-() const;
    Quaternion operator*(const Quaternion& q) const;
    Quaternion operator*(double scalar) const;
    Quaternion operator/(double scalar) const;
};

Quaternion normalize (const Quaternion &q);
Quaternion inverse(const Quaternion &q);
double norm2(const Quaternion &q);
inline std::ostream &operator<< (std::ostream &out, const Quaternion &q) {out << "(" << q.s << ", " << q.v << ")"; return out;}

struct Transformation {
    Vec3 translation;
    double scale;
    Quaternion rotation;
    Transformation (double factor=1);
    Vec3 apply (const Vec3 &point) const;
    Vec3 apply_vec (const Vec3 &vec) const;
    Transformation operator+(const Transformation& t) const;
    Transformation operator-(const Transformation& t) const;
    Transformation operator*(const Transformation& t) const;
    Transformation operator*(double scalar) const;
    Transformation operator/(double scalar) const;
};

namespace tform {
	Transformation identity();
}
Transformation inverse(const Transformation &tr);
inline std::ostream &operator<< (std::ostream &out, const Transformation &t) {out << "(translation: " << t.translation << ", rotation: " << t.rotation << ", scale: " << t.scale << ")"; return out;}

//typedef Spline<Transformation> Motion;
typedef std::pair<Transformation,Transformation> DTransformation;

//void clean_up_quaternions (Motion &motion); // remove sign flips

//Transformation get_trans (const Motion &motion, double t);
//DTransformation get_dtrans (const Motion &motion, double t);
//Vec3 apply_dtrans (const DTransformation &dT, const Vec3 &x0, Vec3 *vel=NULL);
//Vec3 apply_dtrans_vec (const DTransformation &dT, const Vec3 &v0);

#endif
