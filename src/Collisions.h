#pragma once
#ifndef __Collisions__
#define __Collisions__

#include "external\ArcSim\mesh.hpp"
#include "boxTriCollision.h"
#include "Obstacles.h"

void CD(Mesh& mesh, std::shared_ptr<Obstacles> obs, std::vector<std::shared_ptr<btc::Collision> > &cls);

#endif