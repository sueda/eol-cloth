#pragma once
#ifndef __Cloth__
#define __Cloth__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "boxTriCollision.h"
#include "CVM.h"
#include "Brenderable.h"

// ArcSim
#include "mesh.hpp"

class MatrixStack;
class Program;
class Box;
class ChronoTimer;

extern struct Material {
	double density; // area density
	double edge_density;
	//StretchingSamples dde_stretching;
	//BendingData dde_bending;
	double damping; // stiffness-proportional damping coefficient
	double strain_min, strain_max; // strain limits
	double yield_curv, weakening; // plasticity parameters
	double yield_stretch, plastic_flow, plastic_limit;
	bool use_dde; // use DDE material files
	double thickness;
	double alt_stretching, alt_bending, alt_poisson; // alternative material model
	double toughness, fracture_bend_thickness; // fracture toughness
};

extern struct Remeshing {
	double refine_angle, refine_compression, refine_velocity;
	double size_min, size_max, size_uniform; // size limits
	double aspect_min; // aspect ratio control
	double refine_fracture;
};

class Cloth : public Brenderable
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Mesh mesh;
	Remeshing remeshing;
	Material material;

	Mesh last_mesh;

	std::vector<std::shared_ptr<btc::Collision> > collisions;
	//std::vector<std::shared_ptr<CVM> > collisions_passed;
	
	Cloth();
	virtual ~Cloth();

	void setup(const Eigen::Vector3d &x00,
		const Eigen::Vector3d &x01,
		const Eigen::Vector3d &x10,
		const Eigen::Vector3d &x11);
	void setup(const std::string &filename);

	int rows;
	int cols;
	int step;
	int pull_step;
	int side_pull_change_step;
	bool pull_free;
	bool move_wire;

	bool coll;
	bool fixon;
	bool preserve_edge;
	bool matlab_debug_physics;
	bool friction;
	bool wire;
	bool points;
	//bool matlab_debug_collision;
	//bool export_postscript;
	bool export_timings;

	double Xmin, Xmax, Ymin, Ymax;
	double density;
	double e;
	double nu;
	double beta;
	double alpha;
	double alpha_change;
	Eigen::Vector2d damping;

	//double collapse_non_preserve_thresh;
	//double collapse_preserve_thresh;

	Eigen::MatrixXd verts1;
	Eigen::MatrixXd norms1;
	double pmove1;

	Eigen::VectorXi EOLverts;

	Eigen::MatrixXd fixed;
	
	void tare();
	void reset();
	void updatePosNor();
	void updateBuffers();
	void stepL(double h, const Eigen::Vector3d &grav, std::vector<std::shared_ptr<Box> > box, double t);
	bool stepEoL(double h, const Eigen::Vector3d &grav, std::vector<std::shared_ptr<Box> > box, double t, bool equality);
	//bool collisionRemesh(std::shared_ptr<Box> box);
	//void export2D(std::string file_name);
	
	void init();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p, bool on2D) const;

	// Exporting
	int getBrenderCount() const;
	std::vector<std::string> getBrenderNames() const;
	void exportBrender(std::vector< std::shared_ptr< std::ofstream > > outfiles) const;

	bool tmpbool;

	double energyadd;
	
private:

	double safety_margin;
	
	Eigen::VectorXd v_old;
	Eigen::VectorXd v;
	Eigen::VectorXd f;
	
	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> posBuf2D;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
	unsigned eleBufID;
	unsigned posBufID;
	unsigned posBuf2DID;
	unsigned norBufID;
	unsigned texBufID;

	std::vector<std::shared_ptr<ChronoTimer> > ClothTimer;
};

#endif
