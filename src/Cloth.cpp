#include "Cloth.h"
#include "conversions.h"

#include "external\ArcSim\mesh.hpp"
#include "external\ArcSim\io.hpp"

#ifdef EOLC_ONLINE
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "online/MatrixStack.h"
#include "online/Program.h"
#include "online/GLSL.h"
#endif // EOLC_ONLINE


using namespace std;
using namespace Eigen;

void Cloth::build(const Vector2i res,
	const Vector3d &x00,
	const Vector3d &x01,
	const Vector3d &x10,
	const Vector3d &x11)
{
	for (int i = 0; i < res(0); ++i) {
		double u = i / (res(0) - 1.0);
		Vector3d x0 = (1 - u)*x00 + u*x10;
		Vector3d x1 = (1 - u)*x01 + u*x11;
		for (int j = 0; j < res(1); ++j) {
			double v = j / (res(1) - 1.0);
			Vector3d x = (1 - v)*x0 + v*x1;
			mesh.add(new Vert(Vec3(x(0),x(1),0.0), Vec3(0)));
			mesh.add(new Node(e2v(x), e2v(x), Vec3(0), 0, 0, false));
			connect(mesh.verts.back(), mesh.nodes.back());
		}
	}
	double dx = ((x01(0) - x00(0)) / (res(1) - 1)) / 2;
	double dy = ((x10(1) - x00(1)) / (res(0) - 1)) / 2;
	for (int i = 0; i < res(0) - 1; ++i) {
		for (int j = 0; j < res(1) - 1; ++j) {
			Vector3d x = x00;
			x(0) += (2 * j + 1) * dx;
			x(1) += (2 * i + 1) * dy;
			Vec3 ver(x[0], x[1], 0);
			Vec3 vel(0.0, 0.0, 0.0);
			mesh.add(new Vert(Vec3(x(0), x(1), 0.0), Vec3(0)));
			mesh.add(new Node(e2v(x), e2v(x), Vec3(0), 0, 0, false));
			connect(mesh.verts.back(), mesh.nodes.back());
		}
	}
	int center_index_cnt = 0;
	for (int i = 0; i < res(0) - 1; i++) {
		for (int j = 0; j < res(1) - 1; j++) {
			int k0 = (i * res(1)) + j; // upper right index
			int kc0 = ((res(0) * res(1)) + center_index_cnt); // center index
			center_index_cnt++;
			vector<Vert*> verts1;
			verts1.push_back(mesh.verts[k0]);
			verts1.push_back(mesh.verts[k0 + 1]);
			verts1.push_back(mesh.verts[kc0]);
			vector<Face*> faces1 = triangulateARC(verts1);
			for (int f = 0; f < faces1.size(); f++)
				mesh.add(faces1[f]);
			vector<Vert*> verts2;
			verts2.push_back(mesh.verts[k0 + 1]);
			verts2.push_back(mesh.verts[k0 + res(1) + 1]);
			verts2.push_back(mesh.verts[kc0]);
			vector<Face*> faces2 = triangulateARC(verts2);
			for (int f = 0; f < faces2.size(); f++)
				mesh.add(faces2[f]);
			vector<Vert*> verts3;
			verts3.push_back(mesh.verts[k0 + res(1) + 1]);
			verts3.push_back(mesh.verts[k0 + res(1)]);
			verts3.push_back(mesh.verts[kc0]);
			vector<Face*> faces3 = triangulateARC(verts3);
			for (int f = 0; f < faces3.size(); f++)
				mesh.add(faces3[f]);
			vector<Vert*> verts4;
			verts4.push_back(mesh.verts[k0 + res(1)]);
			verts4.push_back(mesh.verts[k0]);
			verts4.push_back(mesh.verts[kc0]);
			vector<Face*> faces4 = triangulateARC(verts4);
			for (int f = 0; f < faces4.size(); f++)
				mesh.add(faces4[f]);
		}
	}
	for (int i = 0; i < mesh.faces.size(); i++) {
		mesh.faces[i]->material = &material;
	}

	mark_nodes_to_preserve(mesh);
	compute_ms_data(mesh);

	v.resize(mesh.nodes.size() * 3);
	f.resize(mesh.nodes.size() * 3);

	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(mesh.nodes.size() * 3);
	norBuf.resize(mesh.nodes.size() * 3);
	texBuf.resize(mesh.nodes.size() * 2);
	eleBuf.resize(mesh.faces.size() * 3);

	updatePosNor();
}

void Cloth::updatePosNor()
{
	for (int i = 0; i < mesh.nodes.size(); i++) {
		Vec3 xm = mesh.nodes[i]->x;
		Vec3 Xm = mesh.nodes[i]->verts[0]->u;
		Vec3 nm = mesh.nodes[i]->n;
		Vec3 tm = mesh.nodes[i]->verts[0]->u;
		posBuf[3 * i + 0] = xm[0];
		posBuf[3 * i + 1] = xm[1];
		posBuf[3 * i + 2] = xm[2];
		norBuf[3 * i + 0] = nm[0];
		norBuf[3 * i + 1] = nm[1];
		norBuf[3 * i + 2] = nm[2];
		texBuf[2 * i + 0] = tm[0];
		texBuf[2 * i + 1] = tm[1];
	}
	for (int i = 0; i < mesh.faces.size(); i++) {
		eleBuf[3 * i + 0] = mesh.faces[i]->v[0]->index;
		eleBuf[3 * i + 1] = mesh.faces[i]->v[1]->index;
		eleBuf[3 * i + 2] = mesh.faces[i]->v[2]->index;
	}
}

void Cloth::updateBuffers()
{
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(mesh.nodes.size() * 3);
	norBuf.resize(mesh.nodes.size() * 3);
	texBuf.resize(mesh.nodes.size() * 2);
	eleBuf.resize(mesh.faces.size() * 3);

	// Update position and normal buffers
	updatePosNor();
}

void Cloth::step(double h)
{

}

#ifdef EOLC_ONLINE
void Cloth::init()
{
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size() * sizeof(float), &texBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	assert(glGetError() == GL_NO_ERROR);
}

void Cloth::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
	// Draw mesh
	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(p->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_DYNAMIC_DRAW);
	glDrawElements(GL_TRIANGLES, eleBuf.size(), GL_UNSIGNED_INT, (void*)0);
	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

void Cloth::drawSimple(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
	for (int i = 0; i < mesh.nodes.size(); i++) {
		if (!mesh.nodes[i]->EoL) {
			glColor3f(1.0f, 1.0f, 0.0f);
			glPointSize(10.0f);
			glBegin(GL_POINTS);
			glVertex3f(mesh.nodes[i]->x[0], mesh.nodes[i]->x[1], mesh.nodes[i]->x[2]);
			glEnd();
		}
	}

	for (int i = 0; i < mesh.edges.size(); i++) {
		Edge *e = mesh.edges[i];
		if (e->preserve) {
			glColor3f(0.0f, 1.0f, 0.0f);
			glLineWidth(3.0f);
			glBegin(GL_LINES);
			glVertex3f(e->n[0]->x[0], e->n[0]->x[1], e->n[0]->x[2]);
			glVertex3f(e->n[1]->x[0], e->n[1]->x[1], e->n[1]->x[2]);
			glEnd();
			glLineWidth(1.0f);
		}

		if (e->n[0]->EoL) {
			glColor3f(1.0f, 0.0f, 0.0f);
			glPointSize(10.0f);
			glBegin(GL_POINTS);
			glVertex3f(e->n[0]->x[0], e->n[0]->x[1], e->n[0]->x[2]);
			glEnd();
		}

		if (e->n[1]->EoL) {
			glColor3f(1.0f, 0.0f, 0.0f);
			glPointSize(10.0f);
			glBegin(GL_POINTS);
			glVertex3f(e->n[1]->x[0], e->n[1]->x[1], e->n[1]->x[2]);
			glEnd();
		}
	}
}
#endif // EOLC_ONLINE
