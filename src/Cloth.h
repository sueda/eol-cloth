#pragma once
#ifndef __Cloth__
#define __Cloth__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

// ArcSim
#include "external/ArcSim//mesh.hpp"

#include "Brenderable.h"

class Obstacles;
class FixedList;
class Constraints;
class Forces;
class GeneralizedSolver;

#ifdef EOLC_ONLINE
class MatrixStack;
class Program;
#endif

extern struct Material {
	double density; // area density
	double e;
	double nu;
	double beta;
	double dampingA;
	double dampingB;
};

extern struct Remeshing {
	double refine_angle, refine_compression, refine_velocity;
	double size_min, size_max; // size limits
	double aspect_min; // aspect ratio control
};

class Cloth : public Brenderable
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Mesh mesh;
	Remeshing remeshing;
	Material material;

	Mesh last_mesh;

	std::vector<std::shared_ptr<FixedList> > fs;

	std::shared_ptr<Constraints> consts;
	std::shared_ptr<Forces> myForces;
	
	Cloth();
	virtual ~Cloth() {};

	void build(const Eigen::Vector2i res,
		const Eigen::Vector3d &p00,
		const Eigen::Vector3d &p01,
		const Eigen::Vector3d &p10,
		const Eigen::Vector3d &p11);
	void build(const std::string &filename);
	void updatePosNor();
	void updateBuffers();

	void updatePreviousMesh();
	void velocityTransfer();
	void step(std::shared_ptr<GeneralizedSolver> gs, std::shared_ptr<Obstacles> obs, const Eigen::Vector3d& grav, double h, const bool& REMESHon, const bool& online);
	void solve(std::shared_ptr<GeneralizedSolver> gs, double h);
	void updateFix(double t);

#ifdef EOLC_ONLINE
	void init();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void drawSimple(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
#endif // EOLC_ONLINE

	// Exporting
	int getBrenderCount() const;
	std::vector<std::string> getBrenderNames() const;
	void exportBrender(std::vector< std::shared_ptr< std::ofstream > > outfiles) const;

private:

	int fsindex;
	
	Eigen::VectorXd v_old;
	Eigen::VectorXd v;
	Eigen::VectorXd f;
	
	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

#ifdef EOLC_ONLINE
	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
#endif // EOLC_ONLINE

};

#endif
