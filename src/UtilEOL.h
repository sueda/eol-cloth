#pragma once
#ifndef __UTILEOL__
#define __UTILEOL__

#include "external\ArcSim\mesh.hpp"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

Eigen::MatrixXd deform_grad(const Face *f);

Eigen::MatrixXd deform_grad_v(const Vert* v);

#endif