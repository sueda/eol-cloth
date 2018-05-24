#pragma once
#ifndef __cONVERSIONS__
#define __CONVERSIONS__

#include "external\ArcSim\vectors.hpp"
#include <Eigen\Core>

Vec2 e2v(const Eigen::Vector2d);
Vec3 e2v(const Eigen::Vector3d);

Eigen::Vector2d v2e(const Vec2);
Eigen::Vector3d v2e(const Vec3);

#endif