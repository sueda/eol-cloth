#include "parseParams.h"

#include "external\Json\json.h"
#include "external\ArcSim\util.hpp"
#include "external\ArcSim\vectors.hpp"
#include "external\ArcSim\referenceshape.hpp"

#include "Cloth.h"
#include "Constraints.h"
#include "FixedList.h"
#include "Obstacles.h"
#include "Points.h"
#include "Box.h"
#include "Shape.h"
#include "GeneralizedSolver.h"

#include <Eigen\Core>

#include <iostream>
#include <vector>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
using namespace Eigen;

// Borrowed from ArcSim

void complain(const Json::Value &json, const string &expected) {
	cout << "Expected " << expected << ", found " << json << " instead" << endl;
	abort();
}

void complainM(int msize, int gotsize) {
	cout << "Expected " << msize << "elements, found " << gotsize << " instead" << endl;
	abort();
}

template <typename T> void parse(vector<T> &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	v.resize(json.size());
	for (int i = 0; i < json.size(); i++)
		parse(v[i], json[i]);
}

template <typename T> void parse(T &x, const Json::Value &json, const T &x0) {
	if (json.isNull())
		x = x0;
	else
		parse(x, json);
}

void parse(bool &b, const Json::Value &json) {
	if (!json.isBool()) complain(json, "boolean");
	b = json.asBool();
}
void parse(int &n, const Json::Value &json) {
	if (!json.isIntegral()) complain(json, "integer");
	n = json.asInt();
}
void parse(double &x, const Json::Value &json) {
	if (!json.isNumeric()) complain(json, "real");
	x = json.asDouble();
}
void parse(string &s, const Json::Value &json) {
	if (!json.isString()) complain(json, "string");
	s = json.asString();
}

struct Range {
	double &min, &max;
	Range(double &min, double &max) : min(min), max(max) {}
};

void parse(Range range, const Json::Value &json, Vec2 range0, const string key) {
	if (json.isNull()) {
		range.min = range0[0];
		range.max = range0[1];
		return;
	}
	if (json.size() != 2) {
		cout << "The array value " << key << "was defined unproperly" << endl;
		cout << "Do so with the JSON key:" << endl;
		cout << "	\"" << key << "\": [#,#]" << endl;
		abort();
	}
	parse(range.min, json[0u]);
	parse(range.max, json[1]);
}

// End borrowed from ArcSim

void parse(Vector2i &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	if (json.size() != 2) complain(json, "array of size" + to_string(2));
	for (int i = 0; i < v.size(); i++) {
		v(i) = json[i].asInt();
	}
}

void parse(Vector2d &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	if (json.size() != 2) complain(json, "array of size" + to_string(2));
	for (int i = 0; i < v.size(); i++) {
		v(i) = json[i].asDouble();
	}
}

void parse(Vector3d &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	if (json.size() != 3) complain(json, "array of size" + to_string(3));
	for (int i = 0; i < v.size(); i++) {
		v(i) = json[i].asDouble();
	}
}

void parse(VectorXd &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	if (v.size() != json.size()) complain(json, "array of size " + to_string(v.size()));
	for (int i = 0; i < v.size(); i++) {
		v(i) = json[i].asDouble();
	}
}

void load_genset(const shared_ptr<genSet> genset, const string &JSON_FILE)
{
	Json::Value json;
	Json::Reader reader;
	ifstream file(JSON_FILE.c_str());
	bool parsingSuccessful = reader.parse(file, json);
	if (!parsingSuccessful) {
		fprintf(stderr, "Error reading file: %s\n", JSON_FILE.c_str());
		fprintf(stderr, "%s", reader.getFormatedErrorMessages().c_str());
		abort();
	}
	file.close();

	parse(genset->online, json["online"], false);
	parse(genset->exportObjs, json["exportObjs"], false);
	parse(genset->exportTimings, json["exportTimings"], false);
	parse(genset->RESOURCE_DIR, json["RESOURCE_DIR"], string(""));
	if (genset->RESOURCE_DIR == "") {
		cout << "Resource directory was not specified." << endl;
		cout << "Do so with the JSON key:" << endl;
		cout << "	\"RESOURCE_DIR\": <path-to-resource-dir>" << endl;
		abort();
	}
	if (genset->exportObjs) {
		parse(genset->OUTPUT_DIR, json["OUTPUT_DIR"], string(""));
		if (genset->OUTPUT_DIR == "") {
			cout << "Export set to on, but no output directory specified." << endl;
			cout << "Do so with the JSON key:" << endl;
			cout << "	\"OUTPUT_DIR\": <path-to-output-dir>" << endl;
			abort();
		}
	}

	genset->RESOURCE_DIR += string("/"); // Just in case

	genset->printGenSet();
}

void load_solver(shared_ptr<GeneralizedSolver> gs, const Json::Value& json)
{
	if (!json.isString()) complain(json, "string");
	string which = json.asString();
	if (which == "mosek") {
#ifndef EOLC_MOSEK
		cout << "ERROR:" << endl;
		cout << "You've specified using the mosek solver, but you haven't built with mosek support" << endl;
		abort();
#endif // !EOLC_MOSEK
		gs->whichSolver = GeneralizedSolver::Mosek;
	}
	else if (which == "gurobi") {
#ifndef EOLC_GUROBI
		cout << "ERROR:" << endl;
		cout << "You've specified using the gurobi solver, but you haven't built with gurobi support" << endl;
		abort();
#endif // !EOLC_GUROBI
		gs->whichSolver = GeneralizedSolver::Gurobi;
	}
	else {
		cout << "Unrecognized solver:" << endl;
		cout << "	\"" << which << "\"" << endl;
		cout << "Check spelling" << endl;
		abort();
	}
}

void load_matset(Material& material, const Json::Value& json)
{
	parse(material.density, json["density"], 0.05);
	parse(material.e, json["youngs"], 50.0);
	parse(material.nu, json["poissons"], 0.01);
	parse(material.beta, json["stiffness"], 1.0e-5);
	parse(Range(material.dampingA, material.dampingB),
		json["damping"], Vec2(0.0, 1.0), "damping");
}

// Borrowed from ArcSim
void load_remeshset(Remeshing& remeshing, const Json::Value& json)
{
	parse(remeshing.refine_angle, json["refine_angle"], infinity);
	parse(remeshing.refine_compression, json["refine_compression"], infinity);
	parse(remeshing.refine_velocity, json["refine_velocity"], infinity);
	parse(Range(remeshing.size_min, remeshing.size_max),
		json["size"], Vec2(-infinity, infinity), "size");
	parse(remeshing.aspect_min, json["aspect_min"], -infinity);
}

void load_fixedset(vector<shared_ptr<FixedList> > &fsv, const Json::Value& json)
{
	VectorXd nothing(6);
	nothing << -1.0, -1.0, -1.0, -1.0, -1.0, -1.0;
	if (json.isMember("0")) {
		auto fs = make_shared<FixedList>();
		parse(fs->when, json["0"]["when"], 0.0);
		parse(fs->c1, json["0"]["corner1"], nothing);
		parse(fs->c2, json["0"]["corner2"], nothing);
		parse(fs->c3, json["0"]["corner3"], nothing);
		parse(fs->c4, json["0"]["corner4"], nothing);
		fsv.push_back(fs);
	}
	else {
		auto fs = make_shared<FixedList>();
		parse(fs->when, json["when"], 0.0);
		parse(fs->c1, json["corner1"], nothing);
		parse(fs->c2, json["corner2"], nothing);
		parse(fs->c3, json["corner3"], nothing);
		parse(fs->c4, json["corner4"], nothing);
		fsv.push_back(fs);
		return;
	}
	int i = 1;
	while(true) {
		if (!json.isMember(to_string(i))) break;
		auto fs = make_shared<FixedList>();
		parse(fs->when, json[to_string(i)]["when"], 0.0);
		parse(fs->c1, json[to_string(i)]["corner1"], nothing);
		parse(fs->c2, json[to_string(i)]["corner2"], nothing);
		parse(fs->c3, json[to_string(i)]["corner3"], nothing);
		parse(fs->c4, json[to_string(i)]["corner4"], nothing);
		fsv.push_back(fs);
		i++;
	}
}

// Individual object settings loaders

void load_defclothset(shared_ptr<Cloth> cloth, const Json::Value& json)
{
	Vector2i res;
	parse(res, json["initial_cloth_res"], Vector2i(2, 2));
	VectorXd p00(5), p01(5), p10(5), p11(5);
	VectorXd dp00(5), dp01(5), dp10(5), dp11(5);
	dp00 << 0.0, 0.0, 0.0, 0.0, 0.0;
	dp10 << 1.0, 0.0, 0.0, 1.0, 0.0;
	dp01 << 0.0, 1.0, 0.0, 0.0, 1.0;
	dp11 << 1.0, 1.0, 0.0, 1.0, 1.0;
	parse(p00, json["corner1"], dp00);
	parse(p01, json["corner2"], dp01);
	parse(p10, json["corner3"], dp10);
	parse(p11, json["corner4"], dp11);
	cloth->build(res, p00, p01, p10, p11);
	cloth->mesh.parent = cloth;
	cloth->mesh.ref = new ReferenceLinear(cloth->mesh);
}

void load_clothset(shared_ptr<Cloth> cloth, const Json::Value& json)
{
	if (json.isMember("cloth_obj") && json["cloth_obj"] != "") {

	}
	else if (json.isMember("init")) {
		load_defclothset(cloth, json["init"]);
	}

	if (json.isMember("Material")) {
		load_matset(cloth->material, json["Material"]);
	}
	else {
		load_matset(cloth->material, json);
	}

	if (json.isMember("Remeshing")) {
		load_remeshset(cloth->remeshing, json["Remeshing"]);
	}
	else {
		load_remeshset(cloth->remeshing, json);
	}

	if (json.isMember("Fixed")) {
		load_fixedset(cloth->fs, json["Fixed"]);
		// I don't like how I did this
		Vector2i res;
		parse(res, json["init"]["initial_cloth_res"], Vector2i(2, 2));
		for (int f = 0; f < cloth->fs.size(); f++) {
			cloth->fs[f]->c1i = 0;
			cloth->fs[f]->c2i = res(0) * (res(1) - 1);
			cloth->fs[f]->c3i = res(0) * res(1) - 1;
			cloth->fs[f]->c4i = res(0) - 1;
		}
	}
}

void load_pointsset(shared_ptr<Points> p, const Json::Value& json)
{
	if (!json.isArray()) complain(json, "array");
	if (p->num_points != 0) {
		MatrixXd pxyzfromfile = p->pxyz;
		MatrixXd normsfromfile = p->norms;
		p->pxyz.resize(3, p->num_points + json.size());
		p->norms.resize(3, p->num_points + json.size());
		p->pxyz.block(0, 0, 3, p->num_points) = pxyzfromfile;
		p->norms.block(0, 0, 3, p->num_points) = normsfromfile;
		for (int i = 0; i < json.size(); i++) {
			if (!json[i].isArray()) complain(json, "array");
			if (json[i].size() != 6) {
				cout << "Invalid point define." << endl;
				cout << "Do so with the JSON key:" << endl;
				cout << "	\"points\": [x, y, z, nx, ny, nz]" << endl;
				abort();
			}
			p->pxyz.block(0, i + p->num_points, 3, 1) = Vector3d(json[i][0].asDouble(), json[i][1].asDouble(), json[i][2].asDouble());
			p->norms.block(0, i + p->num_points, 3, 1) = Vector3d(json[i][3].asDouble(), json[i][4].asDouble(), json[i][5].asDouble());
		}
		p->num_points += json.size();
	}
	else {
		p->num_points = json.size();
		p->pxyz.resize(3, json.size());
		p->norms.resize(3, json.size());
		for (int i = 0; i < json.size(); i++) {
			if (!json[i].isArray()) complain(json, "array");
			if (json[i].size() != 6) {
				cout << "Invalid point define." << endl;
				cout << "Points are formated as:" << endl;
				cout << "	[x, y, z, nx, ny, nz]" << endl;
				abort();
			}
			p->pxyz.block(0, i, 3, 1) = Vector3d(json[i][0].asDouble(), json[i][1].asDouble(), json[i][2].asDouble());
			p->norms.block(0, i, 3, 1) = Vector3d(json[i][3].asDouble(), json[i][4].asDouble(), json[i][5].asDouble());
		}
	}
}

void load_boxset(vector<shared_ptr<Box> > &boxes, shared_ptr<Shape> shape, const Json::Value& json)
{
	if (!json.isArray()) complain(json, "array");
	for (int i = 0; i < json.size(); i++) {
		if (json[i].size() != 12) {
			cout << "Invalid box define." << endl;
			cout << "Boxes are formated as:" << endl;
			cout << "	[scalex, scaley, scalez, x, y, z, agvx, angvy, angvz, vx, vy, vz]" << endl;
			abort();
		}
		auto b = make_shared<Box>(shape,"Box" + to_string(i));
		b->dim = Vector3d(json[i][0].asDouble(), json[i][1].asDouble(), json[i][2].asDouble());
		b->E1.block<3,1>(0,3) = Vector3d(json[i][3].asDouble(), json[i][4].asDouble(), json[i][5].asDouble());
		b->E1inv = b->E1.inverse();
		b->v << json[i][6].asDouble(), json[i][7].asDouble(), json[i][8].asDouble(), json[i][9].asDouble(), json[i][10].asDouble(), json[i][11].asDouble();
		boxes.push_back(b);
	}
}

void load_obsset(shared_ptr<Obstacles> obs, const Json::Value& json)
{
	if(json.isMember("threshold")) parse(obs->cdthreshold, json["threshold"], 5e-3);

	if (json.isMember("points_file")); // TODO

	if (json.isMember("points")) load_pointsset(obs->points, json["points"]);

	if (json.isMember("box_file")); // TODO

	if (json.isMember("boxes")) load_boxset(obs->boxes, obs->shapes[0], json["boxes"]);

	obs->num_boxes = obs->boxes.size();
}

void printSimSet(shared_ptr<Scene> scene);

void load_simset(shared_ptr<Scene> scene, const string &JSON_FILE)
{
	Json::Value json;
	Json::Reader reader;
	ifstream file(JSON_FILE.c_str());
	bool parsingSuccessful = reader.parse(file, json);
	if (!parsingSuccessful) {
		fprintf(stderr, "Error reading file: %s\n", JSON_FILE.c_str());
		fprintf(stderr, "%s", reader.getFormatedErrorMessages().c_str());
		abort();
	}
	file.close();

	if (json.isMember("solver")) load_solver(scene->GS, json["solver"]);

	parse(scene->h, json["timestep"], 0.005);
	parse(scene->grav, json["gravity"], Vector3d(0.0, 0.0, -9.8));
	parse(scene->REMESHon, json["REMESH"], false);
	parse(scene->EOLon, json["EOL"], false);
	if (scene->EOLon) scene->REMESHon = true;

	if (json.isMember("Cloth")) load_clothset(scene->cloth, json["Cloth"]);

	if (json.isMember("Obstacles")) load_obsset(scene->obs, json["Obstacles"]);

	printSimSet(scene);
}

string printSimBool(bool b) {
	if (b) return "True";
	return "False";
}

void printSimSet(shared_ptr<Scene> scene)
{
	cout << "Simulation settings" << endl;
	cout << "	Timestep: " << scene->h << endl;
	cout << "	REMESH: " << printSimBool(scene->REMESHon) << endl;
	cout << "	EOL: " << printSimBool(scene->EOLon) << endl;
	cout << "	Cloth:" << endl;
	cout << "		cloth_obj: " << "" << endl;
	cout << "		Material: " << endl;
	cout << "			density: " << scene->cloth->material.density << endl;
	cout << "			youngs: " << scene->cloth->material.e << endl;
	cout << "			poissons: " << scene->cloth->material.nu << endl;
	cout << "			stiffness: " << scene->cloth->material.beta << endl;
	cout << "			damping: [" << scene->cloth->material.dampingA << ", " << scene->cloth->material.dampingA << "]" << endl;
	cout << "		Remeshing: " << endl;
	cout << "			refine_angle: " << scene->cloth->remeshing.refine_angle << endl;
	cout << "			refine_compression: " << scene->cloth->remeshing.refine_compression << endl;
	cout << "			refine_velocity: " << scene->cloth->remeshing.refine_velocity << endl;
	cout << "			size: [" << scene->cloth->remeshing.size_min << ", " << scene->cloth->remeshing.size_max << "]" << endl;
	cout << "			aspect_min: " << scene->cloth->remeshing.aspect_min << endl;
	cout << "	Obstacles:" << endl;
	cout << "		threshold: " << scene->obs->cdthreshold << endl;
	//cout << "		points_file: " << "" << endl;
	cout << "		total_points: " << scene->obs->points->num_points << endl;
	//cout << "		box_file: " << "" << endl;
	cout << "		total_boxes" << scene->obs->num_boxes << endl;

}