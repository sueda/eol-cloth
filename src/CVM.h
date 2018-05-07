#pragma once
#ifndef __CVM__
#define __CVM__

#include "mesh.hpp"
#include "boxTriCollision.h"

struct CVM
{
	//std::shared_ptr<btc::Collision> coll;
	//int* vert_index;
	Vert* vert;
	int which_weight;
	//CVM(std::shared_ptr<btc::Collision> c, Vert* v) : coll(c), vert(v), which_weight(0) {

	//}
	CVM(Vert* v) : vert(v), which_weight(0) {}

	//CVM(std::shared_ptr<btc::Collision> c, Vert* v, int ww) : coll(c), vert(v), which_weight(ww) {

	//}
	CVM(Vert* v, int ww) : vert(v) , which_weight(ww) {}
};

#endif