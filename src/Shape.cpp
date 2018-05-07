#include "Shape.h"
#include <iostream>

#include "GLSL.h"
#include "Program.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

using namespace std;

Shape::Shape() :
	posBufID(0),
	norBufID(0),
	texBufID(0)
{
}

Shape::~Shape()
{
}

void Shape::loadMesh(const string &meshName)
{
	// Load geometry
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	string errStr;
	bool rc = tinyobj::LoadObj(&attrib, &shapes, &materials, &errStr, meshName.c_str());
	if(!rc) {
		cerr << errStr << endl;
	} else {
		// Some OBJ files have different indices for vertex positions, normals,
		// and texture coordinates. For example, a cube corner vertex may have
		// three different normals. Here, we are going to duplicate all such
		// vertices.
		// Loop over shapes
		for(size_t s = 0; s < shapes.size(); s++) {
			// Loop over faces (polygons)
			size_t index_offset = 0;
			for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
				size_t fv = shapes[s].mesh.num_face_vertices[f];
				// Loop over vertices in the face.
				for(size_t v = 0; v < fv; v++) {
					// access to vertex
					tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+0]);
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+1]);
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+2]);
					if(!attrib.normals.empty()) {
						norBuf.push_back(attrib.normals[3*idx.normal_index+0]);
						norBuf.push_back(attrib.normals[3*idx.normal_index+1]);
						norBuf.push_back(attrib.normals[3*idx.normal_index+2]);
					}
					if(!attrib.texcoords.empty()) {
						texBuf.push_back(attrib.texcoords[2*idx.texcoord_index+0]);
						texBuf.push_back(attrib.texcoords[2*idx.texcoord_index+1]);
					}
				}
				index_offset += fv;
				// per-face material (IGNORE)
				shapes[s].mesh.material_ids[f];
			}
		}
	}
}

void Shape::init()
{
	// Send the position array to the GPU
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_STATIC_DRAW);
	
	// Send the normal array to the GPU
	if(!norBuf.empty()) {
		glGenBuffers(1, &norBufID);
		glBindBuffer(GL_ARRAY_BUFFER, norBufID);
		glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_STATIC_DRAW);
	}
	
	// Send the texture array to the GPU
	if(!texBuf.empty()) {
		glGenBuffers(1, &texBufID);
		glBindBuffer(GL_ARRAY_BUFFER, texBufID);
		glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_STATIC_DRAW);
	}
	
	// Unbind the arrays
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	GLSL::checkError(GET_FILE_LINE);
}

void Shape::draw(const shared_ptr<Program> prog) const
{
	GLSL::checkError(GET_FILE_LINE);
	// Bind position buffer
	int h_pos = prog->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	
	// Bind normal buffer
	int h_nor = prog->getAttribute("aNor");
	if(h_nor != -1 && norBufID != 0) {
		glEnableVertexAttribArray(h_nor);
		glBindBuffer(GL_ARRAY_BUFFER, norBufID);
		glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	}
	
	// Bind texcoords buffer
	//int h_tex = prog->getAttribute("aTex");
	//if(h_tex != -1 && texBufID != 0) {
	//	glEnableVertexAttribArray(h_tex);
	//	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	//	glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	//}
	
	// Draw
	int count = posBuf.size()/3; // number of indices to be rendered
	glDrawArrays(GL_TRIANGLES, 0, count);
	
	// Disable and unbind
	//if(h_tex != -1) {
	//	glDisableVertexAttribArray(h_tex);
	//}
	if(h_nor != -1) {
		glDisableVertexAttribArray(h_nor);
	}
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	GLSL::checkError(GET_FILE_LINE);
}

void Shape::exportBrender(const Eigen::Matrix4d &E, std::ofstream& outfile) const {
	//Make transformation conversions
	// //multiply E by each p(x,y,z) where p is wrt Sphere
	//can i multiply each value in posbuf-->i dont think so.. 
	// -->create a vertex4d for each set of 3 values in posbuf, 
	//		then multiply and use result as value to export



	//vertex positions
	for (int i = 0; i < posBuf.size(); i = i + 3) {
		Eigen::Vector4d pos(posBuf[i], posBuf[i + 1], posBuf[i + 2], 1.0); //double check if it should be zero or 1
																		   // apply transformation
		Eigen::Vector4d trans_pos = E * pos;
		char vert[50];
		//sprintf(vert, "v %f %f %f\n", posBuf[i], posBuf[i + 1], posBuf[i + 2]);
		sprintf(vert, "v %f %f %f\n", trans_pos(0), trans_pos(1), trans_pos(2));
		outfile << vert;
	}
	//texture coordinates
	for (int i = 0; i < texBuf.size(); i = i + 2) {
		char vtex[50];
		sprintf(vtex, "vt %f %f\n", texBuf[i], texBuf[i + 1]);
		outfile << vtex;
	}
	//normal vectors
	for (int i = 0; i < norBuf.size(); i = i + 3) {
		Eigen::Vector4d nor(norBuf[i], norBuf[i + 1], norBuf[i + 2], 0.0);
		// apply transformation
		Eigen::Vector4d trans_nor = E * nor;
		char norm[50];
		sprintf(norm, "vn %f %f %f\n", trans_nor(0), trans_nor(1), trans_nor(2));
		outfile << norm;
	}
	//faces--Using Triangle Strips

	//
	//face
	//f vertex/texture/normal1 vertex/texture/normal2 vertex/texture/normal3 
	// posbuf holds all vertices. each face has 3 vertices.
	for (int i = 1; i < (posBuf.size() / 3); i = i + 3) {
		char face[50];
		int f1, f2, f3;
		f1 = i;
		f2 = i + 1;
		f3 = i + 2;

		sprintf(face, "f %i/%i/%i %i/%i/%i %i/%i/%i\n", f1, f1, f1, f2, f2, f2, f3, f3, f3);
		outfile << face;
	}

	// for (int j = 0; j < rows*2-2; j = j+2) {
	// 	///one row:
	// 	for (int i = 0; i < cols*2-2; i++) {
	// 		char facetri[50];
	// 		int strt = cols*j;
	// 		int v1, v2, v3;
	// 		v1 = eleBuf[strt+ i] + 1;
	// 		v2 = eleBuf[strt+ i + 1] + 1;
	// 		v3 = eleBuf[strt+ i + 2] + 1;

	// 		sprintf(facetri, "f %i/%i/%i %i/%i/%i %i/%i/%i\n", v1, v1, v1, v2, v2, v2, v3, v3, v3);
	// 		outfile << facetri;
	// 	}
	// }
}