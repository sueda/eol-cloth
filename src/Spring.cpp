#include "Spring.h"
#include "Particle.h"

using namespace std;
using namespace Eigen;

Spring::Spring(shared_ptr<Particle> p0, shared_ptr<Particle> p1, shared_ptr<Particle> p2) :
	E(1.0)
{
	assert(p0);
	assert(p1);
	assert(p2);
	assert(p0 != p1);
	assert(p0 != p2);
	assert(p1 != p2);
	this->p0 = p0;
	this->p1 = p1;
	this->p2 = p2;
	//this->xy0 = Eigen::Vector2d(p0->x(0),p0->x(0));
	//this->xy1 = Eigen::Vector2d(p1->x(0), p1->x(0));
	//this->xy2 = Eigen::Vector2d(p2->x(0), p2->x(0));
}

Spring::~Spring()
{
	
}
