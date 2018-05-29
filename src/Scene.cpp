#include "Scene.h"
#include "Cloth.h"
#include "Obstacles.h"
#include "Shape.h"
#include "Preprocessor.h"
#include "Collisions.h"
#include "Constraints.h"
#include "matlabOutputs.h"

#include "external\ArcSim\dynamicremesh.hpp"

#ifdef EOLC_ONLINE
#include "online\MatrixStack.h"
#include "online\Program.h"
#endif // EOLC_ONLINE


using namespace std;
using namespace Eigen;

Scene::Scene() : 
	h(0.005),
	part(0)
{
	cloth = make_shared<Cloth>();
	obs = make_shared<Obstacles>();
	consts = make_shared<Constraints>();
}

void Scene::load(const string &RESOURCE_DIR)
{
	obs->load(RESOURCE_DIR);
}

void Scene::init(const bool online, const bool exportObjs)
{
#ifdef EOLC_ONLINE
	if (online) {
		cloth->init();
		obs->init();
	}
#endif

	consts->init(obs);

	if (exportObjs) {

	}
}

void Scene::step()
{
	if (part != 0) {
		cout << "Please finish the partial step before making a full step" << endl;
		return;
	}
	dynamic_remesh(cloth->mesh);
	set_indices(cloth->mesh);
	CD(cloth->mesh, obs, cls);
	//dynamic_remesh(cloth->mesh);
	preprocess(cloth->mesh, cls);
	//reindex_nodes(cloth->mesh.nodes);
	set_indices(cloth->mesh);
	cout << "pre" << endl;
	dynamic_remesh(cloth->mesh);
	set_indices(cloth->mesh);
	//preprocessClean(cloth->mesh);
	//set_indices(cloth->mesh);
	cloth->updateBuffers();
	obs->step(h);
	cls.clear();
	mesh2m(cloth->mesh, "mesh.m", true);
	cout << "step" << endl;
}

void Scene::partialStep()
{
	if (part == 0) {
		dynamic_remesh(cloth->mesh);
		set_indices(cloth->mesh);
		CD(cloth->mesh, obs, cls);
		cout << "CD" << endl;
	}
	else if (part >= 1 && part < 8) {
		preprocessPart(cloth->mesh, cls, part);
		set_indices(cloth->mesh);
		cloth->updateBuffers();
		mesh2m(cloth->mesh, "mesh.m", true);
	}
	else if (part == 8) {
		dynamic_remesh(cloth->mesh);
		set_indices(cloth->mesh);
		cloth->updateBuffers();
	}
	else if (part == 9) {
		cloth->updateBuffers();
		cout << "Finished step" << endl;
		cls.clear();
		part = -1;
	}
	part++;
}

#ifdef EOLC_ONLINE
void Scene::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const
{
	obs->draw(MV,p);

	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.5, 0.5).data());
	glUniform3fv(p->getUniform("kdBack"), 1, Vector3f(1.0, 0.5, 0.5).data());
	cloth->draw(MV, p);
}

void Scene::drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const
{
	obs->drawSimple(MV, p);
	cloth->drawSimple(MV, p);
}
#endif // EOLC_ONLINE
