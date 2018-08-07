#pragma once
#ifndef __Scene__
#define __Scene__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include "boxTriCollision.h"
#include "BrenderManager.h"

class Cloth;
class Obstacles;
class Constraints;
class Shape;
class GeneralizedSolver;

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
	void init(const bool& online, const bool& exportObjs, const std::string& OUTPUT_DIR);
	void reset();
	void step(const bool& online, const bool& exportObjs);
	void partialStep();

#ifdef EOLC_ONLINE
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
#endif // EOLC_ONLINE

	double h;

	Eigen::Vector3d grav;

	bool EOLon;
	bool REMESHon;

	int part;

	std::shared_ptr<GeneralizedSolver> GS;

	std::shared_ptr<Cloth> cloth;
	std::shared_ptr<Obstacles> obs;
	std::vector<std::shared_ptr<btc::Collision> > cls;
	
private:

	double t;

	// Export
	BrenderManager *brender;

};

#endif
