#include "Obstacles.h"
#include "Points.h"
#include "Box.h"
#include "Shape.h"

#ifdef EOLC_ONLINE
#include "online\MatrixStack.h"
#include "online\Program.h"
#endif // EOLC_ONLINE

using namespace std;
using namespace Eigen;

Obstacles::Obstacles() : 
	num_boxes(0)
{
	points = std::make_shared<Points>();
}

void Obstacles::load(const std::string &RESOURCE_DIR)
{
	auto box_shape = make_shared<Shape>();
	box_shape->loadMesh(RESOURCE_DIR + "box.obj");
	shapes.push_back(box_shape);
}

void Obstacles::step(double h)
{
	for (int i = 0; i < boxes.size(); i++) {
		boxes[i]->step(h);
	}
}

void Obstacles::addExport(BrenderManager *brender)
{
	for (int b = 0; b < num_boxes; b++) {
		brender->add(boxes[b]);
	}
}

#ifdef EOLC_ONLINE
void Obstacles::init()
{
	for (int i = 0; i < boxes.size(); i++) {
		boxes[i]->init();
	}
}

void Obstacles::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const
{
	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(0.0, 1.0, 1.0).data());
	for (int i = 0; i < boxes.size(); i++) {
		boxes[i]->draw(MV,p);
	}
}

void Obstacles::drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const
{
	//for (int i = 0; i < boxes.size(); i++) {
	//	//boxes[i]->drawSimple(MV, p);
	//}
	points->drawSimple(MV, p);
}
#endif // EOLC_ONLINE

