#include "conversions.h"

#include "external\ArcSim\vectors.hpp"

using namespace Eigen;

Vec2 e2v(const Vector2d v)
{
	return Vec2(v(0), v(1));
}

Vec3 e2v(const Vector3d v)
{
	return Vec3(v(0), v(1), v(2));
}

Vector2d v2e(const Vec2 v)
{
	return Vector2d(v[0], v[1]);
}

Vector2d v322e(const Vec3 v)
{
	return Vector2d(v[0], v[1]);
}

Vector3d v2e(const Vec3 v)
{
	return Vector3d(v[0], v[1], v[2]);
}