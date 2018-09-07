#include "Scene.h"
#include "Cloth.h"
#include "Obstacles.h"
#include "Shape.h"
#include "Preprocessor.h"
#include "Collisions.h"
#include "Constraints.h"
#include "GeneralizedSolver.h"
#include "matlabOutputs.h"

#include "external\ArcSim\dynamicremesh.hpp"

#ifdef EOLC_ONLINE
#include "online\MatrixStack.h"
#include "online\Program.h"
#endif // EOLC_ONLINE


using namespace std;
using namespace Eigen;

Scene::Scene() : 
	t(0.0),
	h(0.005),
	grav(Vector3d(0.0,0.0,-9.8)),
	part(0)
{
	cloth = make_shared<Cloth>();
	obs = make_shared<Obstacles>();
	GS = make_shared<GeneralizedSolver>();
}

void Scene::load(const string &RESOURCE_DIR)
{
	obs->load(RESOURCE_DIR);
}

void Scene::init(const bool& online, const bool& exportObjs, const string& OUTPUT_DIR)
{
#ifdef EOLC_ONLINE
	if (online) {
		cloth->init();
		obs->init();
	}
#endif

	if (REMESHon) {
		dynamic_remesh(cloth->mesh);
		set_indices(cloth->mesh);
	}

	cloth->consts->init(obs);

	if (exportObjs) {
		brender = BrenderManager::getInstance();
		brender->setExportDir(OUTPUT_DIR);
		brender->add(cloth);
		obs->addExport(brender);
		brender->exportBrender(t);
	}
}

void printstate(Mesh& mesh)
{
	for (int n = 0; n < mesh.nodes.size(); n++) {
		if (mesh.nodes[n]->EoL) {
			cout << mesh.nodes[n]->EoL_state << endl;
		}
	}
}

void Scene::step(const bool& online, const bool& exportObjs)
{
	cout << "Sim time: " << t << endl;

	if (part != 0) {
		cout << "Please finish the partial step before making a full step" << endl;
		return;
	}
	cloth->updateFix(t);
	if (EOLon) {
		cloth->updatePreviousMesh();
		CD(cloth->mesh, obs, cls);
		preprocess(cloth->mesh, cloth->boundaries, cls);
		//cout << "pre" << endl;
	}
	if (REMESHon) {
		dynamic_remesh(cloth->mesh);
		set_indices(cloth->mesh);
	}
	cloth->step(GS, obs, grav, h, REMESHon, online);
	obs->step(h);
	cls.clear();
	//mesh2m(cloth->mesh, "mesh.m", true);
	if (exportObjs) brender->exportBrender(t);
	//cout << "step" << endl;
	t += h;
}

void Scene::partialStep()
{
	if (part == 0) {
		cloth->updatePreviousMesh();
		cloth->updateFix(t);
		CD(cloth->mesh, obs, cls);
		cout << "CD" << endl;
	}
	else if (part >= 1 && part < 8) {
		preprocessPart(cloth->mesh, cloth->boundaries, cls, part);
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
		cloth->step(GS, obs, grav, h, REMESHon, true);
		obs->step(h);
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
