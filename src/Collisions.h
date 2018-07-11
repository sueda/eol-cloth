#pragma once
#ifndef __Collisions__
#define __Collisions__

#include "external\ArcSim\mesh.hpp"
#include "boxTriCollision.h"
#include "Obstacles.h"

void CD(const Mesh& mesh, const std::shared_ptr<Obstacles> obs, std::vector<std::shared_ptr<btc::Collision> > &cls);

void CD2(const Mesh& mesh, const std::shared_ptr<Obstacles> obs, std::vector<std::shared_ptr<btc::Collision> > &cls);

#endif