#pragma once
#ifndef __Scene__
#define __Scene__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include "boxTriCollision.h"

class Cloth;
class Obstacles;
class Shape;

#ifdef EOLC_ONLINE
class MatrixStack;
class Program;
#endif // EOLC_ONLINE

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		Scene();
	virtual ~Scene() {};
	
	void load(const std::string &RESOURCE_DIR);
	void init(const bool online, const bool exportObjs);
	void reset();
	void step();
	void partialStep();

	

#ifdef EOLC_ONLINE
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
#endif // EOLC_ONLINE

	double h;

	bool EOLon;

	int part;

	std::shared_ptr<Cloth> cloth;
	std::shared_ptr<Obstacles> obs;
	std::vector<std::shared_ptr<btc::Collision> > cls;
	
private:

	double t;

};

#endif
