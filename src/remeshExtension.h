#pragma once
#ifndef __remeshExtension__
#define __remeshExtension__

#include "external\ArcSim\remesh.hpp"

RemeshOp split_edgeForced(Edge* edge, double d);

RemeshOp collapse_edgeForced(Edge* edge, int i);

RemeshOp split_face(Face* face, Vec3 b);

#endif