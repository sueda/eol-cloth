#ifndef REFERENCESHAPE_HPP
#define REFERENCESHAPE_HPP

#include "mesh.hpp"

class ReferenceShape {
public:
    virtual Vec3 closest_point (const Vec3& p) = 0;
    virtual bool raycast(Vec3& p, const Vec3& dir) = 0;    
};

class ReferenceSphere : public ReferenceShape {
public:
    ReferenceSphere(const Vec3& center, double radius);
    ReferenceSphere(const Mesh& mesh);
    
    Vec3 closest_point (const Vec3& p);
    bool raycast(Vec3& p, const Vec3& dir);
private:
    Vec3 center;
    double radius;
};

class ReferenceLinear : public ReferenceShape {
public:
    ReferenceLinear(const Mesh& mesh);
    
    Vec3 closest_point (const Vec3& p);
    bool raycast(Vec3& p, const Vec3& dir);
};

class ReferenceMesh : public ReferenceShape {
public:
    ReferenceMesh(const Mesh& mesh, const std::string& filename);
    
    Vec3 closest_point (const Vec3& p);
    bool raycast(Vec3& p, const Vec3& dir);
};


#endif