#include "Edge.h"
#include "Particle.h"

using namespace std;

Edge_old::Edge_old(shared_ptr<Particle> p0, shared_ptr<Particle> p1, shared_ptr<Particle> p2, shared_ptr<Particle> p3)
{
	assert(p0);
	assert(p1);
	assert(p2);
	assert(p3);
	assert(p0 != p1);
	assert(p0 != p2);
	assert(p0 != p3);
	assert(p1 != p2);
	assert(p1 != p3);
	assert(p2 != p3);
	this->p0 = p0;
	this->p1 = p1;
	this->p2 = p2;
	this->p3 = p3;
}

Edge_old::~Edge_old()
{

}
