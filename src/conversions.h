#pragma once
#ifndef __CONVERSIONS__
#define __CONVERSIONS__

#include "external\ArcSim\vectors.hpp"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen\Core>

// These conversions translate ArcSim vectors to Eigen vectors

// Eigen 2 vec
Vec2 e2v(const Eigen::Vector2d);
Vec3 e2v(const Eigen::Vector3d);

// Vec 2 Eigen
Eigen::Vector2d v2e(const Vec2);
Eigen::Vector3d v2e(const Vec3);

// Vec 3D to Eigen 2D
Eigen::Vector2d v322e(const Vec3 v);

#endif