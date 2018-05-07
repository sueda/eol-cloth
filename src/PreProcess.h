#pragma once
#ifndef __PreProcess__
#define __PreProcess__

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include "boxTriCollision.h"
#include "CVM.h"

#include "mesh.hpp"

class Box;
class ChronoTimer;

class pp_settings
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	double safety_margin;
	double collapse_non_preserve_thresh;
	double collapse_preserve_thresh;
	bool matlab_debug_collision;
	bool export_postscript;
	bool EOLon;
	bool wire;
	bool points;
	bool move_points;
	bool export_timings;
	bool nudge;

	Eigen::Vector2d bounds;

	double h;
	double pmove1;
	Eigen::MatrixXd verts1;
	Eigen::MatrixXd norms1;

	std::vector<std::shared_ptr<ChronoTimer> > PPTimer;

	bool once_hack;

	pp_settings() :
		EOLon(false),
		matlab_debug_collision(false),
		export_postscript(false),
		safety_margin(0.01),
		once_hack(false),
		nudge(false) {}

	pp_settings(double cn, double c, bool eolon, bool mdc, bool ep) :
		collapse_non_preserve_thresh(cn),
		collapse_preserve_thresh(c),
		EOLon(eolon),
		matlab_debug_collision(mdc),
		export_postscript(ep),
		safety_margin(0.01),
		once_hack(false),
		nudge(false)
	{}
};

void export2D(const Mesh& mesh, std::string file_name);

void updateEdgeWeights(Mesh& mesh, std::vector<std::shared_ptr<Box> > box);

bool collisionRemesh(Mesh& mesh, std::vector<std::shared_ptr<Box> > box, std::vector<std::shared_ptr<btc::Collision> > collisions, Eigen::VectorXi EOLS, std::shared_ptr<pp_settings> pps);

bool collisionRemeshAlt(Mesh& mesh, std::vector<std::shared_ptr<Box> > box, std::vector<std::shared_ptr<btc::Collision> > collisions, Eigen::VectorXi EOLS, std::shared_ptr<pp_settings> pps);

#endif