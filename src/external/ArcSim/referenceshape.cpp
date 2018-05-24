#include "referenceshape.hpp"

using namespace std;

ReferenceSphere::ReferenceSphere(const Vec3& _center, double _radius) : center(_center), radius(_radius) {
}

// assume mesh is a discretized sphere; find center and radius using LSQ
ReferenceSphere::ReferenceSphere(const Mesh& mesh) {
    Mat3x3 M(0);
    Vec3 b(0);
    for (size_t i=0; i<mesh.nodes.size(); i++) {
        Vec3 pi = mesh.nodes[i]->x;
        for (size_t j=0; j<mesh.nodes.size(); j++) {
            Vec3 pj = mesh.nodes[j]->x;
            for (int k=0; k<3; k++) M.col(k) += pi[k]*(pj-pi);
            b += 0.5*pi*(norm2(pj)-norm2(pi));
        }
    }
    center = M.inv() * b;
    double R2 = 0;
    for (size_t i=0; i<mesh.nodes.size(); i++) {
        R2 += norm2(center - mesh.nodes[i]->x);
    }
    radius = sqrt(R2/mesh.nodes.size());    
}

Vec3 ReferenceSphere::closest_point(const Vec3& p) {
    return center + radius * normalize(p - center);
}

bool ReferenceSphere::raycast(Vec3& p, const Vec3& dir) {
    p = normalize(p-center) * radius + center;
    return true;
}

ReferenceLinear::ReferenceLinear(const Mesh& mesh) {    
}


Vec3 ReferenceLinear::closest_point(const Vec3& p) {
    return p;
}

bool ReferenceLinear::raycast(Vec3& p, const Vec3& dir) {
    return true;
}

ReferenceMesh::ReferenceMesh(const Mesh& mesh, const string& filename) {  
    cout << "Reference mesh lookup not implemented here." << endl;
    exit(1);
}

Vec3 ReferenceMesh::closest_point(const Vec3& p) {
    return p;
}

bool ReferenceMesh::raycast(Vec3& p, const Vec3& dir) {
    return false;
}
