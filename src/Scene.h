#pragma once
#ifndef __Scene__
#define __Scene__

#include <vector>
#include <memory>
#include <string>

// ArcSim
#include "mesh.hpp"
#include "BrenderManager.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Cloth;
class Particle;
class MatrixStack;
class Program;
class Shape;
class Box;
class Env_obj;
class Mesh;
class ChronoTimer;
class pp_settings;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();

	bool EoLon;
	
	void load(const std::string &RESOURCE_DIR, const std::string &JSON_FILE);
	void init();
	void tare();
	void reset();
	void step();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
	
	double getTime() const { return t; }

	std::shared_ptr<Cloth> cloth;
	std::vector<std::shared_ptr<Box> > boxes;
	std::shared_ptr<pp_settings> pps;

	std::vector<std::shared_ptr<ChronoTimer> > MasterTimer;
	
private:

	void load_json(const std::string &JSON_FILE);
	void parse_points_file(Eigen::MatrixXd &m, const std::string &POINTS_FILE);

	bool remeshon;
	bool staticon;
	bool dynamicon;
	bool collon;
	bool fixon;
	bool vert_cloth;
	bool edge_preserveon;
	bool export_objs;
	bool use_points_file;

	double t;
	double h;

	int step_count;
	int how_many_boxes;

	Eigen::Vector3d grav;
	Eigen::Vector2i res;
	Eigen::Vector2d dim;
	Eigen::Vector3d pos;

	std::shared_ptr<Shape> box_shape;
	std::shared_ptr<Shape> north;
	std::shared_ptr<Shape> south;
	std::shared_ptr<Shape> east;
	std::shared_ptr<Shape> west;

	std::vector< std::shared_ptr<Env_obj> > env_objs;
	std::vector<Mesh*> box_mesh_vec;

	// Export
	BrenderManager *brender;
	BrenderManager *brender_preprocessed;
	BrenderManager *brender_ArcSimed;
};

#endif
