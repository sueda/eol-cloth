#pragma once
#ifndef __Spring__
#define __Spring__

#include <memory>

//#include <Eigen/Dense>

class Particle;

class Spring
{
public:
	Spring(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);
	virtual ~Spring();
	
	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	std::shared_ptr<Particle> p2;
	//Eigen::Vector2d xy0;
	//Eigen::Vector2d xy1;
	//Eigen::Vector2d xy2;
	double E;
};

#endif
