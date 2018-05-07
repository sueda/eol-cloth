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

#include "transformation.hpp"

using namespace std;

Transformation tform::identity () {
    return Transformation();
}

Transformation inverse(const Transformation &tr) {
    Transformation in;
    in.scale = 1.f / tr.scale;
    in.rotation = inverse(tr.rotation);
    in.translation = Vec3(0.) -
        in.rotation.rotate(
            in.scale * (
                tr.translation
            )
        );
    return in;
}

Quaternion inverse(const Quaternion &q) {
    Quaternion in;
    double divisor = norm2(q);
    in.s = q.s / divisor;
    in.v = -q.v / divisor;
    return in;
}

Quaternion Quaternion::from_axisangle(const Vec3 &axis, double angle) {
    Quaternion q;
    if (angle == 0) {
        q.s = 1;
        q.v = Vec3(0);
    } else {
        q.s = cos(angle/2);
        q.v = sin(angle/2)*normalize(axis);
    }
    return q;
}

pair<Vec3, double> Quaternion::to_axisangle() const {
    double angle = 2 * acos(s);
    Vec3 axis;
    if(angle == 0) {
        axis = Vec3(1);
    } else {
        axis = v / sqrt(1.0-s*s);
    }
    return pair<Vec3, double>(axis, angle);
}

Transformation::Transformation(double factor) {
    translation = Vec3(0);
    scale = factor;
    rotation = Quaternion::from_axisangle(Vec3(1), 0)*factor;
}

Transformation Transformation::operator-(const Transformation& other) const {
    Transformation t;
    t.scale = this->scale - other.scale;
    t.translation = this->translation - other.translation;
    t.rotation = this->rotation - other.rotation;
    return t;
}

Transformation Transformation::operator+(const Transformation& other) const {
    Transformation t;
    t.scale = this->scale + other.scale;
    t.translation = this->translation + other.translation;
    t.rotation = this->rotation + other.rotation;
    return t;
}

Transformation Transformation::operator*(const Transformation& other) const {
    Transformation t;
    t.scale = this->scale * other.scale;
    t.translation = this->translation + 
                    this->rotation.rotate(other.translation * this->scale);
    t.rotation = this->rotation * other.rotation;
    return t;
}

Transformation Transformation::operator*(double s) const {
    Transformation t;
    t.scale = this->scale * s;
    t.translation = this->translation * s;
    t.rotation = this->rotation * s;
    return t;
}

Transformation Transformation::operator/(double s) const {
    return (*this)*(1./s);
}

Quaternion Quaternion::operator+(const Quaternion& other) const {
    Quaternion q;
    q.v = this->v + other.v;
    q.s = this->s + other.s;
    return q;
}

Quaternion Quaternion::operator-(const Quaternion& other) const {
    Quaternion q;
    q.v = this->v - other.v;
    q.s = this->s - other.s;
    return q;
}

Quaternion Quaternion::operator-() const {
    Quaternion q;
    q.v = -this->v;
    q.s = -this->s;
    return q;
}

Quaternion Quaternion::operator*(const Quaternion& other) const {
    Quaternion q;
    q.v = (this->s * other.v) + (other.s * this->v) +
               cross(this->v, other.v);
    q.s = (this->s * other.s) - dot(this->v, other.v);
    return q;
}

Quaternion Quaternion::operator*(double s) const {
    Quaternion q;
    q.v = this->v * s;
    q.s = this->s * s;
    return q;
}

Quaternion Quaternion::operator/(double s) const {
    return (*this)*(1./s);
}

Vec3 Quaternion::rotate (const Vec3 &x) const {
    return x*(sq(s) - dot(v,v)) +
           2.*v*dot(v,x) + 2.*cross(v,x)*s;
}

Vec3 Transformation::apply (const Vec3 &x) const {
    return translation + scale*rotation.rotate(x);
}

Vec3 Transformation::apply_vec (const Vec3 &v) const {
    return rotation.rotate(v);
}

double norm2(const Quaternion &q) {
    return sq(q.s) + norm2(q.v);
}

Quaternion normalize (const Quaternion &q) {
    double norm = sqrt(norm2(q));
    Quaternion p;
    p.s = q.s/norm;
    p.v = q.v/norm;
    return p;
}

void clean_up_quaternions (Motion &motion) {
    for (int p = 1; p < motion.points.size(); p++) {
        const Quaternion &q0 = motion.points[p-1].x.rotation;
        Quaternion &q1 = motion.points[p].x.rotation;
        double d = dot(q0.v, q1.v) + q0.s*q1.s;
        if (d < 0)
            q1 = -q1;
    }
}

Transformation get_trans (const Motion &motion, double t) {
    Transformation T = motion.pos(t);
    T.rotation = normalize(T.rotation);
    return T;
}

DTransformation get_dtrans (const Motion &motion, double t) {
    Transformation T = motion.pos(t), dT = motion.vel(t);
    Quaternion q = T.rotation, dq = dT.rotation;
    double qq = sq(q.s) + norm2(q.v),
           qdq = q.s*dq.s + dot(q.v, dq.v);
    double normq = sqrt(qq);
    T.rotation = q/normq;
    dT.rotation = dq/normq - q/normq*qdq/qq;
    return make_pair(T, dT);
}

Vec3 apply_dtrans (const DTransformation &dtrans, const Vec3 &x0, Vec3 *vel) {
    const Transformation &T = dtrans.first, &dT = dtrans.second;
    Vec3 x = T.apply(x0);
    if (vel) {
        Vec3 w = 2.*(dT.rotation*inverse(T.rotation)).v;
        *vel = dT.translation + dT.scale*T.rotation.rotate(x0)
             + T.scale*cross(w, T.rotation.rotate(x0));
    }
    return x;
}

Vec3 apply_dtrans_vec (const DTransformation &dtrans, const Vec3 &v0) {
    return dtrans.first.apply_vec(v0);
}
