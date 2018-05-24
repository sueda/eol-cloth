#pragma once
#ifndef __Obstacles__
#define __Obstacles__

#include <vector>
#include <memory>

class Points;
class Box;
class Shape;

#ifdef EOLC_ONLINE
class MatrixStack;
class Program;
#endif // EOLC_ONLINE

class Obstacles
{
public:

	Obstacles();
	virtual ~Obstacles() {};

	int num_boxes;
	double cdthreshold;
	std::shared_ptr<Points> points;
	std::vector<std::shared_ptr<Box> > boxes;
	std::vector<std::shared_ptr<Shape> > shapes;

	void load(const std::string &RESOURCE_DIR);
	void step(double h);

#ifdef EOLC_ONLINE
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void init();
#endif // EOLC_ONLINE



private:

};

#endif