#pragma once
#ifndef BOXTRICOLLISION_H_
#define BOXTRICOLLISION_H_

#include <memory>
#include <vector>
#include <Eigen/Dense>

// Throughout this code, `1` refers to the box and `2` refers to the cloth.
// For example, `pos1` is the collision point on the box, and `pos2` is the
// collision point on the cloth.

namespace btc
{

	/**
	* The vertex ordering is as follows:
	*
	*      x2
	*     /  \
	*    / t0 \
	*   /      \
	* x0--edge--x1
	*   \      /
	*    \ t1 /
	*     \  /
	*      x3
	*/
	class Edge
	{
	public:
		Edge();
		virtual ~Edge();

		Eigen::Vector4i verts;
		Eigen::Vector2i faces;
		bool internal;
		double angle;
		Eigen::Vector3d normals[2];
	};

	void createEdges(
		std::vector<std::shared_ptr<Edge> > &edges, // output
		const Eigen::MatrixXi &faces, // input
		const Eigen::MatrixXd &verts  // input
	);

	class Collision
	{
	public:
		Collision();
		virtual ~Collision();

		// The collision can be Vert-Face, Face-Vert, or Edge-Edge.
		// Vert-Face:
		//   count1 = 1
		//   count2 = 3
		//   verts1 = [i -1 -1], where i is the vertex index for box
		//   verts2 = [i j k], where i, j, k are the vertex indices for cloth
		//   weights1 = [a 0 0], where a is the box vertex weight for i
		//   weights2 = [a b c], where a, b, c, are the cloth vertex weights for i, j, k
		//   tri1 = -1
		//   tri2 = i, where i is the triangle index for cloth
		//   edge1 = [a b c ...], wher a, b, c, etc. are the edge indices of the surrounding edges
		//   edge2 = -1
		// Edge-Edge:
		//   count1 = 2
		//   count2 = 2
		//   verts1 = [i j -1], where i, j are the vertex indices for box
		//   verts2 = [i j -1], where i, j are the vertex indices for cloth
		//   weights1 = [a b 0], where a, b are the box vertex weight for i, j
		//   weights2 = [a b 0], where a, b are the cloth vertex weights for i, j
		//   tri1 = -1
		//   tri2 = -1
		//   edge1 = [i], where i is the edge index for box
		//   edge2 = i, where i is the edge index for cloth
		// Face-Vert:
		//   count1 = 3
		//   count2 = 1
		//   verts1 = [i j k], where i, j, k are the vertex indices for box
		//   verts2 = [i -1 -1], where i is the vertex index for cloth
		//   weights1 = [a b c], where a, b, c, are the cloth vertex weights for i, j, k
		//   weights2 = [a 0 0], where a is the box vertex weight for i
		//   tri1 = i, where i is the triangle index for box
		//   tri2 = -1
		//   edge1 = []
		//   edge2 = -1

		double dist; // distance
		Eigen::Vector3d nor1; // normal
		Eigen::Vector3d nor2; // normal on cloth
		Eigen::Vector3d pos1; // position on box
		Eigen::Vector3d pos2; // position on cloth
		Eigen::Vector3d pos1_; // position on box offset to be slightly inside
		int count1; // Collision type for box:   1=vert, 2=edge, 3=face
		int count2; // Collision type for cloth: 1=vert, 2=edge, 3=face
		Eigen::Vector3i verts1; // The collided vertices on box
		Eigen::Vector3i verts2; // The collided vertices on cloth
		Eigen::Vector3d weights1; // The vertex weights on box
		Eigen::Vector3d weights2; // The vertex weights on cloth
		int tri1; // Triangle index for box
		int tri2; // Triangle index for cloth
		std::vector<int> edge1; // all surrounding edge indices for he box
		int edge2; // edge index for cloth
	};

	/**
	* Main function
	* OUTPUT
	*   collisions: a vector of detected collisions
	* INPUT
	*   threshold: detection threshold (e.g., 1e-5)
	*   whd1: box dimensions
	*   E1: box frame
	*   verts2: 3xn matrix of cloth vertices
	*   faces2: 3xm matrix of cloth faces
	*/
	void boxTriCollision(
		std::vector<std::shared_ptr<Collision> > &collisions,
		double threshold,
		const Eigen::Vector3d &whd1,
		const Eigen::Matrix4d &E1,
		const Eigen::MatrixXd &verts2,
		const Eigen::MatrixXi &faces2,
		const Eigen::VectorXi &isEOL2,
		bool justVerts2);

	/**
	* Same as above, but doesn't create the edge structure. If the edge structure
	* is already created, then call this function. Otherwise, call the function
	* above, which creates the edge structure based on verts2 and faces2.
	*/
	void boxTriCollision(
		std::vector<std::shared_ptr<Collision> > &collisions,
		double threshold,
		const Eigen::Vector3d &whd1,
		const Eigen::Matrix4d &E1,
		const Eigen::MatrixXd &verts2,
		const Eigen::MatrixXi &faces2,
		const Eigen::VectorXi &isEOL2,
		bool justVerts2,
		const std::vector<std::shared_ptr<Edge> > &edges2);

	///////////////////////////////////////////////////////////////////////////////

	void pointTriCollision(
		std::vector<std::shared_ptr<Collision> > &collisions,
		double threshold,
		const Eigen::MatrixXd &verts1,
		const Eigen::MatrixXd &norms1,
		const Eigen::MatrixXd &verts2,
		const Eigen::MatrixXi &faces2,
		bool justVerts2);

}

#endif
