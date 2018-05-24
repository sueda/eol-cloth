#pragma once
#ifndef _SHAPE_H_
#define _SHAPE_H_

#include <string>
#include <vector>
#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#ifdef EOLC_ONLINE
class Program;
#endif

/**
 * A shape defined by a list of triangles
 * - posBuf should be of length 3*ntris
 * - norBuf should be of length 3*ntris (if normals are available)
 * - texBuf should be of length 2*ntris (if texture coords are available)
 * posBufID, norBufID, and texBufID are OpenGL buffer identifiers.
 */
class Shape
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Shape();
	virtual ~Shape();
	void loadMesh(const std::string &meshName);
	void exportBrender(const Eigen::Matrix4d& E, std::ofstream& outfile) const;

#ifdef EOLC_ONLINE
	void init();
	void draw(const std::shared_ptr<Program> prog) const;
#endif // EOLC_ONLINE

	
private:
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

#endif
