#pragma once
#ifndef __Preprocessor__
#define __Preprocessor__

#include "external\ArcSim\mesh.hpp"
#include "boxTriCollision.h"

void preprocess(Mesh& mesh, const Eigen::MatrixXd &boundaries, std::vector<std::shared_ptr<btc::Collision> > cls);

void preprocessPart(Mesh& mesh, const Eigen::MatrixXd &boundaries, const std::vector<std::shared_ptr<btc::Collision> > cls, int &part);

void preprocessClean(Mesh& mesh);

#endif