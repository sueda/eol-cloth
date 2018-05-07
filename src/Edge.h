#pragma once
#ifndef __Edge__
#define __Edge__

#include <memory>

//#include <Eigen/Dense>

class Particle;

class Edge_old
{
public:
	Edge_old(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, std::shared_ptr<Particle> p3);
	virtual ~Edge_old();

	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	std::shared_ptr<Particle> p2;
	std::shared_ptr<Particle> p3;
};

#endif