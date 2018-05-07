#include <iostream>
#include <vector>
#include <omp.h>
#include <fstream>
#include <sstream>

#include "Scene.h"
#include "Cloth.h"
#include "Box.h"
#include "PreProcess.h"
#include "Environment_Obj.h"
#include "Shape.h"
#include "Program.h"
#include "boxTriCollision.h"
#include "ChronoTimer.h"

// JSON parsing
#include "json.h"

// ArcSim
#include "mesh.hpp"
#include "dynamicremesh.hpp"
#include "referenceshape.hpp"
#include "util.hpp"
#include "io.hpp"
#include "nearobs.hpp"
#include "geometry.hpp"

using namespace std;
using namespace Eigen;

void reorient_MS(Mesh& mesh) {
	Plane plane = plane_fit<MS>(mesh);
	Mat3x3 M = local_base(plane.n);
	for (size_t i = 0; i<mesh.verts.size(); i++) {
		mesh.verts[i]->u = M * (mesh.verts[i]->u - plane.x0) + plane.x0;
		if (fabs(mesh.verts[i]->u[2]) > 1e-4) {
			cout << "error: DDE materials only support planar rest shapes." << endl;
			exit(1);
		}
	}
}

void reproject_all(Mesh& mesh) {
	for (size_t i = 0; i<mesh.verts.size(); i++) {
		Vert* vert = mesh.verts[i];
		vert->u = mesh.ref->closest_point(vert->u);
	}
}

Scene::Scene() :
	t(0.0),
	h(1e-3),
	grav(0.0, 0.0, 0.0)
{
	MasterTimer.push_back(make_shared<ChronoTimer>("PrePr"));
	MasterTimer.push_back(make_shared<ChronoTimer>("ArcSim"));
	MasterTimer.push_back(make_shared<ChronoTimer>("Physics"));
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR, const string &JSON_FILE)
{
	// Units: meters, kilograms, seconds 

	step_count = 0;
	how_many_boxes = 2;

	cloth = make_shared<Cloth>();
	if (collon) {
		box_shape = make_shared<Shape>();
		box_shape->loadMesh(RESOURCE_DIR + "box.obj");
		//for (int i = 0; i < how_many_boxes; i++) {
		//	auto boxShape = make_shared<Shape>();
		//	boxShape->loadMesh(RESOURCE_DIR + "box.obj");
		//	boxes.push_back(make_shared<Box>(boxShape));
		//	load_obj(boxes.back()->mesh, RESOURCE_DIR + "box1_65.obj");
		//	box_mesh_vec.push_back(&boxes.back()->mesh);
		//}
	}

	pps = make_shared<pp_settings>();

	load_json(JSON_FILE);

	if (cloth->points && use_points_file) {
		parse_points_file(cloth->verts1, RESOURCE_DIR + "points.p");
		parse_points_file(pps->verts1, RESOURCE_DIR + "points.p");
	}

	//Vector3d x00(pos(0)-dim(0)/2, pos(1) + dim(1) / 2, pos(2));
	//Vector3d x01(pos(0) + dim(0)/2, pos(1) + dim(1) / 2, pos(2));
	//Vector3d x10(pos(0) - dim(0)/2, pos(1) - dim(1) / 2, pos(2));
	//Vector3d x11(pos(0) + dim(0)/2, pos(1) - dim(1) / 2, pos(2));
	Vector3d x00(pos(0), pos(1), pos(2));
	Vector3d x01, x10, x11;
	if (vert_cloth) {
		x01 = Vector3d(pos(0) + dim(0), pos(1), pos(2));
		x10 = Vector3d(pos(0), pos(1), pos(2) + dim(1));
		x11 = Vector3d(pos(0) + dim(0), pos(1), pos(2) + dim(1));
	}
	else {
		x01 = Vector3d(pos(0) + dim(0), pos(1), pos(2));
		x10 = Vector3d(pos(0), pos(1) + dim(1), pos(2));
		x11 = Vector3d(pos(0) + dim(0), pos(1) + dim(1), pos(2));
	}

	cloth->rows = res(0);
	cloth->cols = res(1);
	cloth->setup(x00, x01, x10, x11);
	//cloth->setup(RESOURCE_DIR + "cloth2k.obj");
	cloth->mesh.parent = cloth;

	cloth->mesh.ref = new ReferenceLinear(cloth->mesh);
	reorient_MS(cloth->mesh);
	reproject_all(cloth->mesh);

	//for (int i = 0; i < cloth->mesh.edges.size(); i++) {
	//	cloth->mesh.edges[i]->preserve = true;
	//}

	//cloth->compute_masses();
	if (remeshon && staticon) {
		static_remesh(cloth->mesh);
		cloth->updateBuffers();
	}

	cloth->pull_free = false;
	cloth->move_wire = true;
	pps->move_points = true;

	pps->PPTimer.push_back(make_shared<ChronoTimer>("CD1"));
	pps->PPTimer.push_back(make_shared<ChronoTimer>("PrePr"));
	if (cloth->export_timings) {
		MasterTimer[0]->export_csv_header();
	}

	// Help with debugging
	north = make_shared<Shape>();
	north->loadMesh(RESOURCE_DIR + "north.obj");
	auto nor = make_shared<Env_obj>(north);
	env_objs.push_back(nor);
	nor->x = Vector3d(-0.2, 0.0, 1.0);

	east = make_shared<Shape>();
	east->loadMesh(RESOURCE_DIR + "east.obj");
	auto eas = make_shared<Env_obj>(east);
	env_objs.push_back(eas);
	eas->x = Vector3d(-0.7, 0.0, 0.1);

	south = make_shared<Shape>();
	south->loadMesh(RESOURCE_DIR + "south.obj");
	auto sou = make_shared<Env_obj>(south);
	env_objs.push_back(sou);
	sou->x = Vector3d(-0.18, 0.0, -0.7);

	west = make_shared<Shape>();
	west->loadMesh(RESOURCE_DIR + "west.obj");
	auto wes = make_shared<Env_obj>(west);
	env_objs.push_back(wes);
	wes->x = Vector3d(0.7, 0.0, -0.3);
}

void Scene::init()
{
	if (collon) {
		for (int i = 0; i < boxes.size(); i++) {
			boxes[i]->init();
		}
	}
	cloth->init();
	north->init();
	east->init();
	south->init();
	west->init();

	// Export
	if (export_objs) {
		brender = BrenderManager::getInstance();
		brender->setExportDir("objs");
		brender->add(cloth);
		//brender->add(boxes[0]);

		//brender_preprocessed = BrenderManager::getInstance();
		//brender_preprocessed->setExportDir("objs_Preprocessed");
		//brender_preprocessed->add(cloth);
		//brender_preprocessed->add(box);

		//brender_ArcSimed = BrenderManager::getInstance();
		//brender_ArcSimed->setExportDir("objs_ArcSimed");
		//brender_ArcSimed->add(cloth);
		//brender_ArcSimed->add(box);
	}
}

void Scene::tare()
{
	cloth->tare();
}

void Scene::reset()
{
	t = 0.0;
	cloth->reset();
}

MatrixXd deform_grad2(const Face *f) {
	MatrixXd Dx(3, 2);
	Dx(0, 0) = f->v[1]->node->x[0] - f->v[0]->node->x[0];
	Dx(0, 1) = f->v[2]->node->x[0] - f->v[0]->node->x[0];
	Dx(1, 0) = f->v[1]->node->x[1] - f->v[0]->node->x[1];
	Dx(1, 1) = f->v[2]->node->x[1] - f->v[0]->node->x[1];
	Dx(2, 0) = f->v[1]->node->x[2] - f->v[0]->node->x[2];
	Dx(2, 1) = f->v[2]->node->x[2] - f->v[0]->node->x[2];
	Matrix2d DX;
	DX(0, 0) = f->v[1]->u[0] - f->v[0]->u[0];
	DX(0, 1) = f->v[2]->u[0] - f->v[0]->u[0];
	DX(1, 0) = f->v[1]->u[1] - f->v[0]->u[1];
	DX(1, 1) = f->v[2]->u[1] - f->v[0]->u[1];
	return Dx * DX.inverse();
}

void Scene::step()
{
	t += h;
	cloth->step = step_count;
	
	if (step_count == 0) cloth->pull_free = true;
	//if (step_count == cloth->pull_step) cloth->energyadd = 1.05;
	//if (step_count == cloth->pull_step) pps->nudge = true;
	//if (step_count == 1077) {
	//	for (int i = 0; i < cloth->mesh.nodes.size(); i++) {
	//		if (cloth->mesh.nodes[i]->EoL) {
	//			cloth->mesh.nodes[i]->preserve = false;
	//			MatrixXd F = MatrixXd::Zero(3, 2);
	//			Vector3d v = Vector3d(cloth->mesh.nodes[i]->v[0], cloth->mesh.nodes[i]->v[1], cloth->mesh.nodes[i]->v[2]);;
	//			Vector2d V = Vector2d(cloth->mesh.nodes[i]->verts[0]->v[0], cloth->mesh.nodes[i]->verts[0]->v[1]);
	//			for (int j = 0; j < cloth->mesh.nodes[i]->verts[0]->adjf.size(); j++) {
	//				F += incedent_angle(cloth->mesh.nodes[i]->verts[0], cloth->mesh.nodes[i]->verts[0]->adjf[j]) * deform_grad2(cloth->mesh.nodes[i]->verts[0]->adjf[j]);
	//			}
	//			is_seam_or_boundary(cloth->mesh.nodes[i]) ? F *= (1 / M_PI) : F *= (1 / (2 * M_PI));
	//			//F /= mesh.nodes[i]->verts[0]->adjf.size();
	//			Vector3d newV = v - F*V;
	//			cloth->mesh.nodes[i]->v = Vec3(newV(0), newV(1), newV(2));
	//			cloth->mesh.nodes[i]->verts[0]->v = Vec3(0.0, 0.0, 0.0);
	//		}
	//	}
	//	for (int i = 0; i < cloth->mesh.edges.size(); i++) {
	//		cloth->mesh.edges[i]->preserve = false;
	//	}
	//	edge_preserveon = false;
	//	EoLon = false;
	//}
	//if (step_count == 10) remeshon = false;
	//if (step_count == 990) {
	//	EoLon = false;
	//	edge_preserveon = false;
	//}
	//if (step_count == 20) {
	//	boxes[0]->v(2) = 5.0;
	//}
	//if (step_count == 500) {
	//	//boxes[0]->v(2) = 0.0;
	//	boxes[0]->v(5) = 0.0;
	//}
	//if (step_count == 670) {
	//	boxes[0]->v(2) = 0.0;
	//	boxes[0]->v(5) = 0.0;
	//}
	//if (step_count == 700) {
	//	boxes[0]->v(1) = 5.0;
	//}
	if (step_count == 1000) { 
		cloth->move_wire = false; 
		pps->move_points = false;
		boxes[0]->v(5) = 0.0;
	}
	cloth->side_pull_change_step = 0;
	if (step_count >= 500) {
		//if (step_count == 250) {
		//	cloth->mesh.nodes[0]->v = Vec3(0);
		//	cloth->mesh.nodes[1]->v = Vec3(0);
		//	cloth->mesh.nodes[2]->v = Vec3(0);
		//	cloth->mesh.nodes[3]->v = Vec3(0);
		//	for (int i = 0; i < cloth->mesh.nodes.size(); i++) {
		//		if (is_seam_or_boundary(cloth->mesh.nodes[i])) cloth->mesh.nodes[i]->v = Vec3(0);
		//	}
		//}
		cloth->fixed.block<1, 3>(4, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(4, 3) = Vector3d(0, 1, 0);
		cloth->fixed.block<1, 3>(5, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(5, 3) = Vector3d(0, -1, 0);
		cloth->fixed.block<1, 3>(6, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(6, 3) = Vector3d(0, -1, 0);
		cloth->fixed.block<1, 3>(7, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(7, 3) = Vector3d(0, 1, 0);
	}
	if (step_count >= 1300) {
		if (step_count == 1300) {
			cloth->mesh.nodes[0]->v = Vec3(0);
			cloth->mesh.nodes[1]->v = Vec3(0);
			cloth->mesh.nodes[2]->v = Vec3(0);
			cloth->mesh.nodes[3]->v = Vec3(0);
			for (int i = 0; i < cloth->mesh.nodes.size(); i++) {
				if (is_seam_or_boundary(cloth->mesh.nodes[i])) cloth->mesh.nodes[i]->v = Vec3(0);
			}
		}
		cloth->fixed.block<1, 3>(4, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(4, 3) = Vector3d(0, 0, 0);
		cloth->fixed.block<1, 3>(5, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(5, 3) = Vector3d(0, 0, 0);
		cloth->fixed.block<1, 3>(6, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(6, 3) = Vector3d(0, 0, 0);
		cloth->fixed.block<1, 3>(7, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(7, 3) = Vector3d(0, 0, 0);
	}
	//if (step_count >= 900) {
	//	if (step_count == 900) {
	//		cloth->mesh.nodes[0]->v = Vec3(0);
	//		cloth->mesh.nodes[1]->v = Vec3(0);
	//		cloth->mesh.nodes[2]->v = Vec3(0);
	//		cloth->mesh.nodes[3]->v = Vec3(0);
	//		for (int i = 0; i < cloth->mesh.nodes.size(); i++) {
	//			if (is_seam_or_boundary(cloth->mesh.nodes[i])) cloth->mesh.nodes[i]->v = Vec3(0);
	//		}
	//	}
	//	cloth->fixed.block<1, 3>(4, 0) = Vector3d(1, 1, 1);
	//	cloth->fixed.block<1, 3>(4, 3) = Vector3d(0, 0, 0);
	//	cloth->fixed.block<1, 3>(5, 0) = Vector3d(1, 1, 1);
	//	cloth->fixed.block<1, 3>(5, 3) = Vector3d(0, 0, 0);
	//	cloth->fixed.block<1, 3>(6, 0) = Vector3d(1, 1, 1);
	//	cloth->fixed.block<1, 3>(6, 3) = Vector3d(0, 0, 0);
	//	cloth->fixed.block<1, 3>(7, 0) = Vector3d(1, 1, 1);
	//	cloth->fixed.block<1, 3>(7, 3) = Vector3d(0, 0, 0);
	//}
	if (step_count >= 1400) {
		if (step_count == 1400) {
			cloth->mesh.nodes[0]->v = Vec3(0);
			cloth->mesh.nodes[1]->v = Vec3(0);
			cloth->mesh.nodes[2]->v = Vec3(0);
			cloth->mesh.nodes[3]->v = Vec3(0);
			for (int i = 0; i < cloth->mesh.nodes.size(); i++) {
				if (is_seam_or_boundary(cloth->mesh.nodes[i])) cloth->mesh.nodes[i]->v = Vec3(0);
			}
		}
		cloth->fixed.block<1, 3>(4, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(4, 3) = Vector3d(0, 0, -1);
		cloth->fixed.block<1, 3>(5, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(5, 3) = Vector3d(0, 0, 1);
		cloth->fixed.block<1, 3>(6, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(6, 3) = Vector3d(0, 0, 1);
		cloth->fixed.block<1, 3>(7, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(7, 3) = Vector3d(0, 0, -1);
	}
	if (step_count >= 2000) {
		if (step_count == 2000) {
			cloth->mesh.nodes[0]->v = Vec3(0);
			cloth->mesh.nodes[1]->v = Vec3(0);
			cloth->mesh.nodes[2]->v = Vec3(0);
			cloth->mesh.nodes[3]->v = Vec3(0);
			for (int i = 0; i < cloth->mesh.nodes.size(); i++) {
				if (is_seam_or_boundary(cloth->mesh.nodes[i])) cloth->mesh.nodes[i]->v = Vec3(0);
			}
		}
		cloth->fixed.block<1, 3>(4, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(4, 3) = Vector3d(0, 0, 0);
		cloth->fixed.block<1, 3>(5, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(5, 3) = Vector3d(0, 0, 0);
		cloth->fixed.block<1, 3>(6, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(6, 3) = Vector3d(0, 0, 0);
		cloth->fixed.block<1, 3>(7, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(7, 3) = Vector3d(0, 0, 0);
	}
	if (step_count >= 2100) {
		if (step_count == 2100) {
			cloth->mesh.nodes[0]->v = Vec3(0);
			cloth->mesh.nodes[1]->v = Vec3(0);
			cloth->mesh.nodes[2]->v = Vec3(0);
			cloth->mesh.nodes[3]->v = Vec3(0);
			for (int i = 0; i < cloth->mesh.nodes.size(); i++) {
				if (is_seam_or_boundary(cloth->mesh.nodes[i])) cloth->mesh.nodes[i]->v = Vec3(0);
			}
		}
		cloth->fixed.block<1, 3>(4, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(4, 3) = Vector3d(0, 0, 1);
		cloth->fixed.block<1, 3>(5, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(5, 3) = Vector3d(0, 0, -1);
		cloth->fixed.block<1, 3>(6, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(6, 3) = Vector3d(0, 0, -1);
		cloth->fixed.block<1, 3>(7, 0) = Vector3d(1, 1, 1);
		cloth->fixed.block<1, 3>(7, 3) = Vector3d(0, 0, 1);
	}
	//if (step_count >= 950) {
	//	if (step_count == 950) {
	//		cloth->mesh.nodes[0]->v = Vec3(0);
	//		cloth->mesh.nodes[1]->v = Vec3(0);
	//		cloth->mesh.nodes[2]->v = Vec3(0);
	//		cloth->mesh.nodes[3]->v = Vec3(0);
	//		for (int i = 0; i < cloth->mesh.nodes.size(); i++) {
	//			if (is_seam_or_boundary(cloth->mesh.nodes[i])) cloth->mesh.nodes[i]->v = Vec3(0);
	//		}
	//	}
	//	cloth->fixed.block<1, 3>(4, 0) = Vector3d(1, 1, 1);
	//	cloth->fixed.block<1, 3>(4, 3) = Vector3d(0, 0, 0);
	//	cloth->fixed.block<1, 3>(5, 0) = Vector3d(1, 1, 1);
	//	cloth->fixed.block<1, 3>(5, 3) = Vector3d(0, 0, 0);
	//	cloth->fixed.block<1, 3>(6, 0) = Vector3d(1, 1, 1);
	//	cloth->fixed.block<1, 3>(6, 3) = Vector3d(0, 0, 0);
	//	cloth->fixed.block<1, 3>(7, 0) = Vector3d(1, 1, 1);
	//	cloth->fixed.block<1, 3>(7, 3) = Vector3d(0, 0, 0);
	//}
	//cout << cloth->fixed << endl;

	//if (step_count > 800) {
	//	//cloth->fixed()
	//	//grav = Vector3d(0.0,0.0,0.0);
	//	cloth->damping = Vector2d(0.0, 0.0);
	//}

	if ((remeshon && dynamicon) || edge_preserveon) {
		//vector<Plane> planes = nearest_obstacle_planes(cloth->mesh, boxes[0]->mesh_vec);

		//vector<AccelStruct*> obs_accs = create_accel_structs(box_mesh_vec, false);
		//map<Node*, Plane> planes = nearest_obstacle_planes(cloth->mesh.nodes, obs_accs);
		//destroy_accel_structs(obs_accs);
		//dynamic_remesh(cloth->mesh, planes);
		//reindex_nodes(cloth->mesh.nodes);

		if (EoLon) {
			delete_mesh(cloth->last_mesh);
			cloth->last_mesh = deep_copy(cloth->mesh);
		}

		// Collision Remesh
		if (edge_preserveon) {
			MasterTimer[0]->tic();
			if (collisionRemesh(cloth->mesh, boxes, cloth->collisions, cloth->EOLverts, pps)) {
				//reorient_MS(cloth->mesh);
				//reproject_all(cloth->mesh);
				reindex_nodes(cloth->mesh.nodes, EoLon);
			}
			MasterTimer[0]->toc();
			MasterTimer[0]->print();
		}

		//cloth->updateBuffers();
		//if (export_objs) brender_preprocessed->exportBrender(t);

		if (pps->export_postscript) export2D(cloth->mesh, "collapsed.ps");

		if (remeshon) {
			MasterTimer[1]->tic();
			vector<AccelStruct*> obs_accs = create_accel_structs(box_mesh_vec, false);
			map<Node*, Plane> planes = nearest_obstacle_planes(cloth->mesh.nodes, obs_accs);
			destroy_accel_structs(obs_accs);
			dynamic_remesh(cloth->mesh, planes);
			reindex_nodes(cloth->mesh.nodes, EoLon);
			MasterTimer[1]->toc();
			MasterTimer[1]->print();
			if (cloth->export_timings) MasterTimer[1]->export_csv();
		}

		//cloth->updateBuffers();
		//if (export_objs) brender_ArcSimed->exportBrender(t);

		//for (int i = 0; i < cloth->mesh.faces.size(); i++) {
		//	cout << aspect(cloth->mesh.faces[i]) << endl;;
		//}

		if (pps->export_postscript) export2D(cloth->mesh, "post_remesh.ps");
	}

	MasterTimer[2]->tic();
	if (EoLon) {
		bool success = cloth->stepEoL(h, grav, boxes, t, false);
		if (!success) {
			cout << "EQUALIYT" << endl;
			cloth->stepEoL(h, grav, boxes, t, true);
		}
	}
	else {
		cloth->stepL(h, grav, boxes, t);
	}
	MasterTimer[2]->toc();
	MasterTimer[2]->print();

	updateEdgeWeights(cloth->mesh, boxes); // Edge weights

	if (pps->export_postscript) export2D(cloth->mesh, "post_physics.ps");

	cloth->updateBuffers();

	for (int i = 0; i < boxes.size(); i++) {
		boxes[i]->step(h);
	}

	step_count++;

	if (cloth->export_timings) MasterTimer[1]->export_csv_endline();

	if (export_objs) brender->exportBrender(t);
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(0.0, 1.0, 1.0).data());
	for (int i = 0; i < boxes.size(); i++) {
		boxes[i]->draw(MV, prog);
	}

	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 0.5, 0.5).data());
	glUniform3fv(prog->getUniform("kdBack"), 1, Vector3f(1.0, 0.5, 0.5).data());
	for (int i = 0; i < (int)env_objs.size(); i++) {
		//env_objs[i]->draw(MV, prog);
	}

	cloth->draw(MV, prog, true);
}

void Scene::parse_points_file(MatrixXd &m, const std::string &POINTS_FILE)
{
	string line;
	int lines = 0;
	ifstream file(POINTS_FILE);
	if (file.is_open()) {
		while (getline(file, line)) {
			lines++;
		}
		file.close();
	}
	ifstream file2(POINTS_FILE);
	m.resize(3, lines);
	lines = 0;
	if (file2.is_open()) {
		while (getline(file2, line)) {
			Vector3d p = Vector3d::Zero();;
			int i = 0;
			istringstream iss(line);
			for (string s; iss >> s; ) {
				p(i) = stod(s);
				i++;
			}
			m.col(lines) = p;
			lines++;
		}
		file2.close();
	}
	cloth->norms1.resize(3, lines);
	pps->norms1.resize(3, lines);
	Vector3d n(0.0, 0.0, 1.0);
	for (int i = 0; i < lines; i++) {
		cloth->norms1.col(i) = n;
		pps->norms1.col(i) = n;
	}
	cloth->pmove1 = 0.0;
	pps->pmove1 = 0.0;
}

//-----------------------------------

void complain(const Json::Value &json, const string &expected) {
	cout << "Expected " << expected << ", found " << json << " instead" << endl;
	abort();
}

void complainM(int msize, int gotsize) {
	cout << "Expected " << msize << "elements, found " << gotsize << " instead" << endl;
	abort();
}

template <int n> void parse(Vec<n> &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	assert(json.size() == n);
	for (int i = 0; i < n; i++)
		parse(v[i], json[i]);
}

template <typename T> void parse(T &x, const Json::Value &json, const T &x0) {
	if (json.isNull())
		x = x0;
	else
		parse(x, json);
}

// TODO:: Template these
void parse(Vector2i &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	for (int i = 0; i < v.size(); i++) {
		v(i) = json[i].asInt();
	}
}

void parse(Vector2d &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	for (int i = 0; i < v.size(); i++) {
		v(i) = json[i].asDouble();
	}
}

void parse(Vector3d &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	for (int i = 0; i < v.size(); i++) {
		v(i) = json[i].asDouble();
	}
}

void parse(MatrixXd &m, const Json::Value &json, int count) {
	if (!json.isArray() && !json.isString()) complain(json, "array");
	if (json.isArray()) {
		if (json.size() != 6) complainM(6, json.size());
		for (int i = 0; i < json.size(); i++) {
			m(count, i) = json[i].asDouble();
		}
	}
	else if (json.isString()) {
		for (int i = 0; i < 6; i++) {
			m(count, i) = -1;
		}
	}
}

void parse_points(MatrixXd &m, const Json::Value &json, int count) {
	if (!json.isArray() && !json.isString()) complain(json, "array");
	if (json.isArray()) {
		if (json.size() != 3) complainM(3, json.size());
		for (int i = 0; i < json.size(); i++) {
			m(i, count) = json[i].asDouble();
		}
	}
}

void parse_transform(Matrix4d &r, const Json::Value &json, const Json::Value &json2) {
	if (!json.isArray()) complain(json, "array");
	double tx = (M_PI * json[0].asDouble()) / 180;
	double ty = (M_PI * json[1].asDouble()) / 180;
	double tz = (M_PI * json[2].asDouble()) / 180;
	r(0, 0) = cos(ty) * cos(tz);
	r(0, 1) = -cos(ty) * sin(tz);
	r(0, 2) = sin(ty);
	r(0, 3) = json2[0].asDouble();
	r(1, 0) = sin(tx) * sin(ty) * cos(tz) + cos(tx) * sin(tz);
	r(1, 1) = -sin(tx) * sin(ty) * sin(tz) + cos(tx) * cos(tz);
	r(1, 2) = -sin(tx) * cos(ty);
	r(1, 3) = json2[1].asDouble();
	r(2, 0) = -cos(tx) * sin(ty) * cos(tz) + sin(tx) * sin(tz);
	r(2, 1) = cos(tx) * sin(ty) * sin(tz) + sin(tx) * cos(tz);
	r(2, 2) = cos(tx) * cos(ty);
	r(2, 3) = json2[2].asDouble();
	r(3, 0) = 0;
	r(3, 1) = 0;
	r(3, 2) = 0;
	r(3, 3) = 1;
}

void parse_box_v(VectorXd &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	if (json.size() != 6) complainM(6, json.size());
	for (int i = 0; i < v.size(); i++) {
		v(i) = json[i].asDouble();
	}
}

void parse_box_rot(Vector3d &r, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	for (int i = 0; i < r.size(); i++) {
		r(i) = (M_PI * json[i].asDouble()) / 180;
	}
}
//

template <typename T> void parse(vector<T> &v, const Json::Value &json) {
	if (!json.isArray()) complain(json, "array");
	v.resize(json.size());
	for (int i = 0; i < json.size(); i++)
		parse(v[i], json[i]);
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

void parse(Range range, const Json::Value &json, Vec2 range0) {
	if (json.isNull()) {
		range.min = range0[0];
		range.max = range0[1];
		return;
	}
	assert(json.size() == 2);
	parse(range.min, json[0u]);
	parse(range.max, json[1]);
}

void Scene::load_json(const std::string &JSON_FILE)
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
	parse(dim, json["initial_cloth_dim"]);
	parse(pps->bounds, json["initial_cloth_dim"]);
	parse(res, json["initial_cloth_res"]);
	parse(pos, json["initial_cloth_pos"]);
	parse(cloth->material.density, json["density"]);
	parse(vert_cloth, json["vert_cloth"]);
	parse(cloth->material.edge_density, json["edge_density"]);
	parse(cloth->e, json["youngs"]);
	parse(cloth->nu, json["poissons"]);
	parse(cloth->beta, json["stiffness"]);
	parse(cloth->damping, json["damping"]);
	parse(grav, json["gravity"]);
	parse(h, json["timestep"]);
	parse(pps->h, json["timestep"]);

	parse(cloth->export_timings, json["Export_timing"]);
	parse(pps->export_timings, json["Export_timing"]);

	parse(EoLon, json["EoL_on"]);
	parse(pps->EOLon, json["EoL_on"]);

	parse(cloth->friction, json["Friction"]);
	parse(cloth->alpha, json["alpha"]);
	parse(cloth->alpha_change, json["alpha_change"]);

	parse(remeshon, json["remeshing_on"]);
	parse(staticon, json["static_on"]);
	parse(dynamicon, json["dynamic_on"]);
	if (remeshon) {
		parse(cloth->remeshing.refine_angle, json["refine_angle"], infinity);
		parse(cloth->remeshing.refine_compression, json["refine_compression"], infinity);
		parse(cloth->remeshing.refine_velocity, json["refine_velocity"], infinity);
		parse(Range(cloth->remeshing.size_min, cloth->remeshing.size_max),
			json["size"], Vec2(-infinity, infinity));
		parse(cloth->remeshing.aspect_min, json["aspect_min"], -infinity);
	}

	parse(fixon, json["fixed_on"]);
	parse(cloth->fixon, json["fixed_on"]);
	if (fixon) {
		parse(cloth->pull_step, json["pull_step"]);
		cloth->fixed = MatrixXd::Zero(8, 6);
		parse(cloth->fixed, json["N"], 0);
		parse(cloth->fixed, json["E"], 1);
		parse(cloth->fixed, json["S"], 2);
		parse(cloth->fixed, json["W"], 3);
		parse(cloth->fixed, json["NE"], 4);
		parse(cloth->fixed, json["SE"], 5);
		parse(cloth->fixed, json["SW"], 6);
		parse(cloth->fixed, json["NW"], 7);
	}

	parse(collon, json["collision_on"]);
	parse(cloth->coll, json["collision_on"]);
	if (collon) {
		parse(how_many_boxes, json["how_many_boxes"]);
		for (int i = 0; i < how_many_boxes; i++) {
			boxes.push_back(make_shared<Box>(box_shape));
		}
		// Boxes
		// 1
		parse(boxes[0]->dim, json["box_dim"]);
		parse_transform(boxes[0]->E1, json["rot_xyz"], json["box_xyz"]);
		boxes[0]->E1inv = boxes[0]->E1.inverse();
		parse_box_rot(boxes[0]->rot, json["rot_xyz"]);
		parse(boxes[0]->x, json["box_xyz"]);
		parse_box_v(boxes[0]->v, json["box_v"]);
		parse(boxes[0]->thresh, json["box_threshold"]);
		// 2
		if (how_many_boxes > 1) {
			parse(boxes[1]->dim, json["box_dim2"]);
			parse_transform(boxes[1]->E1, json["rot_xyz2"], json["box_xyz2"]);
			parse_box_rot(boxes[1]->rot, json["rot_xyz2"]);
			parse(boxes[1]->x, json["box_xyz2"]);
			parse(boxes[1]->thresh, json["box_threshold2"]);
		}
		// 3
		if (how_many_boxes > 2) {
			parse(boxes[2]->dim, json["box_dim3"]);
			parse_transform(boxes[2]->E1, json["rot_xyz3"], json["box_xyz3"]);
			parse_box_rot(boxes[2]->rot, json["rot_xyz3"]);
			parse(boxes[2]->x, json["box_xyz3"]);
			parse(boxes[2]->thresh, json["box_threshold3"]);
		}
		// 4
		if (how_many_boxes > 3) {
			parse(boxes[3]->dim, json["box_dim4"]);
			parse_transform(boxes[3]->E1, json["rot_xyz4"], json["box_xyz4"]);
			parse_box_rot(boxes[3]->rot, json["rot_xyz4"]);
			parse(boxes[3]->x, json["box_xyz4"]);
			parse(boxes[3]->thresh, json["box_threshold4"]);
		}

		// Points
		int how_many_points;
		parse(cloth->points, json["use_points"]);
		parse(pps->points, json["use_points"]);
		parse(use_points_file, json["use_points_file"]);
		if (!use_points_file) {
			parse(how_many_points, json["how_many_points"]);
			cloth->verts1.resize(3, how_many_points);
			pps->verts1.resize(3, how_many_points);
			cloth->norms1.resize(3, how_many_points);
			pps->norms1.resize(3, how_many_points);
			for (int i = 0; i < how_many_points; i++) {
				string point = "point" + to_string(i);
				string norm = "norm" + to_string(i);
				string pmove = "pmove" + to_string(i);
				parse_points(cloth->verts1, json[point], i);
				parse_points(pps->verts1, json[point], i);
				parse_points(cloth->norms1, json[norm], i);
				parse_points(pps->norms1, json[norm], i);
				parse(cloth->pmove1, json[pmove]);
				parse(pps->pmove1, json[pmove]);
			}
		}

		parse(cloth->wire, json["use_wire"]);
		parse(pps->wire, json["use_wire"]);
		if (cloth->points && cloth->wire) {
			cout << "Points and Wire on" << endl;
			exit(EXIT_FAILURE);
		}
	}

	parse(edge_preserveon, json["preserve_edge"]);
	parse(cloth->preserve_edge, json["preserve_edge"]);
	if (edge_preserveon) {
		parse(pps->collapse_non_preserve_thresh, json["collapse_non_preserve_thresh"]);
		parse(pps->collapse_preserve_thresh, json["collapse_preserve_thresh"]);
	}

	parse(cloth->matlab_debug_physics, json["matlab_debug_physics"]);
	if (edge_preserveon) parse(pps->matlab_debug_collision, json["matlab_debug_collision"]);

	if (edge_preserveon) parse(pps->export_postscript, json["export_postscript"]);
	parse(export_objs, json["export_objs"]);
}