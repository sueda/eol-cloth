#pragma once
#ifndef __Preprocessor__
#define __Preprocessor__

#include "external\ArcSim\mesh.hpp"
#include "boxTriCollision.h"

void preprocess(Mesh& mesh, std::vector<std::shared_ptr<btc::Collision> > cls);

void preprocessPart(Mesh& mesh, std::vector<std::shared_ptr<btc::Collision> > cls, int &part);

void preprocessClean(Mesh& mesh);

#endif