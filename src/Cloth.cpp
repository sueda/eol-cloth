#include "Cloth.h"
#include "conversions.h"
#include "Obstacles.h"
#include "FixedList.h"
#include "Constraints.h"
#include "Forces.h"
#include "GeneralizedSolver.h"
#include "UtilEOL.h"
#include "matlabOutputs.h"

#include "external\ArcSim\mesh.hpp"
#include "external\ArcSim\io.hpp"
#include "external\ArcSim\geometry.hpp"

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

Cloth::Cloth() :
	fsindex(0)
{
	consts = make_shared<Constraints>();
	myForces = make_shared<Forces>();
}

void Cloth::build(const Vector2i res,
	const Vector3d &x00,
	const Vector3d &x01,
	const Vector3d &x10,
	const Vector3d &x11)
{
	// Set boundary values
	// In this fixed 4 corner case, the outer boundary is made up of those 4 corners
	boundaries.resize(3, 4);
	boundaries.block<3, 1>(0, 0) = x00;
	boundaries.block<3, 1>(0, 1) = x01;
	boundaries.block<3, 1>(0, 2) = x11; // we flip this order to they are sequential of a continueous boundary path
	boundaries.block<3, 1>(0, 3) = x10;

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

void Cloth::updatePreviousMesh()
{
	delete_mesh(last_mesh);
	last_mesh = deep_copy(mesh);
}

void Cloth::velocityTransfer()
{
	v.resize(mesh.nodes.size() * 3 + mesh.EoL_Count * 2);
	v.setZero();

	// Loop through all of our nodes and update their velocities
	for (int n = 0; n < mesh.nodes.size(); n++) {
		Node* node = mesh.nodes[n];
		Vert* vert = node->verts[0];

		// Search through the unremeshed mesh to see if this node existed before, or was just introduced it this step
		// TODO:: Can more than one node fall within the close range? 
		bool found = false;
		double how_close = 1e-6;
		int closest = -1;
		for (int j = 0; j < last_mesh.nodes.size(); j++) {
			double dist = unsigned_vv_distance(vert->u, last_mesh.nodes[j]->verts[0]->u);
			if (unsigned_vv_distance(vert->u, last_mesh.nodes[j]->verts[0]->u) < how_close) {
				how_close = dist;
				closest = j;
				found = true;
				break;
			}
		}

		// If the node existed before do one of two things
		if (found) {
			// If the node was EOL but is no longer, we have to transfer it's EOL components back into fully LAG
			if (node->EoL_state == Node::WasEOL) {
				node->EoL_state = Node::IsLAG;
				MatrixXd F = MatrixXd::Zero(3, 2);
				Vector3d nodev = v2e(node->v);
				Vector2d nodeV = v322e(vert->v);
				for (int j = 0; j < vert->adjf.size(); j++) {
					F += incedent_angle(vert, vert->adjf[j]) * deform_grad(vert->adjf[j]);
				}
				is_seam_or_boundary(node) ? F *= (1 / M_PI) : F *= (1 / (2 * M_PI));
				Vector3d newvL = nodev - F*nodeV;
				v(3 * n) = newvL(0);
				v(3 * n + 1) = newvL(1);
				v(3 * n + 2) = newvL(2);
			}
			// If it already existed, just pull it's old information
			else {
				// TODO:: The node will have the data if it gets here correct?
				v(3 * n) = node->v[0];
				v(3 * n + 1) = node->v[1];
				v(3 * n + 2) = node->v[2];
				if (node->EoL) {
					v(mesh.nodes.size() * 3 + node->EoL_index * 2) = vert->v[0];
					v(mesh.nodes.size() * 3 + node->EoL_index * 2 + 1) = vert->v[1];
				}
			}
		}

		// If the node did not exist We have to do one of four things
		else {
			// If its a LAG point we can just use barycentric averaging
			if (!node->EoL) {
				Face* old_face = get_enclosing_face(last_mesh, Vec2(vert->u[0], vert->u[1]));
				Vec3 bary = get_barycentric_coords(Vec2(vert->u[0], vert->u[1]), old_face);
				Vector3d vwA = v2e(old_face->v[0]->node->v);
				Vector3d vwB = v2e(old_face->v[1]->node->v);
				Vector3d vwC = v2e(old_face->v[2]->node->v);
				// If any of its surrounding points are EOL, the Eulerian component velocity component must be taken into account
				if (old_face->v[0]->node->EoL) vwA += -deform_grad(old_face) * v322e(old_face->v[0]->v);
				if (old_face->v[1]->node->EoL) vwB += -deform_grad(old_face) * v322e(old_face->v[1]->v);
				if (old_face->v[2]->node->EoL) vwC += -deform_grad(old_face) * v322e(old_face->v[2]->v);
				Vector3d v_new_world = bary[0] * vwA + bary[1] * vwB + bary[2] * vwC;
				node->v = e2v(v_new_world);
				v(3 * n) = v_new_world(0);
				v(3 * n + 1) = v_new_world(1);
				v(3 * n + 2) = v_new_world(2);
			}
			else {
				// If the new EOL was the result of a conserved edge split, we average the components for efficiency 
				// For this to happen neither adjacent nodes can be NewEOLs, otherwise the averaging would be wrong
				if (node->EoL_state == Node::NewEOLFromSplit) {
					Vector3d newvL = Vector3d::Zero();
					Vector3d newvE = Vector3d::Zero();
					for (int e = 0; e < node->adje.size(); e++) {
						if (node->adje[e]->preserve) {
							Node* adjn = other_node(node->adje[e], node);
							newvL += v2e(adjn->v);
							newvE += v2e(adjn->verts[0]->v);
						}
					}
					newvL /= 2.0;
					newvE /= 2.0; // Should only ever be two here..?

					// We put in the new averaged velocities
					v(mesh.nodes.size() * 3 + node->EoL_index * 2) = newvE(0);
					v(mesh.nodes.size() * 3 + node->EoL_index * 2 + 1) = newvE(1);
					v(3 * n) = newvL(0);
					v(3 * n + 1) = newvL(1);
					v(3 * n + 2) = newvL(2);

				}
				// If this is brand new EOL with no previuosly known information, we calculate its world velocity and distribute it
				else {
					Face* old_face = get_enclosing_face(last_mesh, Vec2(vert->u[0], vert->u[1]));
					Matrix2d ftf = deform_grad(old_face).transpose() * deform_grad(old_face);
					Vector2d dtv = -deform_grad(old_face).transpose() * v2e(node->v);

					// We set up a simple KKT system with the purpose of putting as much world velocity as possible into the Eulerian component
					MatrixXd KKTl = MatrixXd::Zero(4, 4);
					KKTl.block(0, 0, 2, 2) = ftf;
					KKTl.block<2, 2>(0, 2) = Matrix2d::Identity();
					KKTl.block<2, 2>(2, 0) = Matrix2d::Identity();
					VectorXd KKTr(4);
					KKTr << dtv, VectorXd::Zero(2);
					ConjugateGradient<MatrixXd, Lower | Upper> cg;
					cg.compute(KKTl);
					VectorXd newvE = cg.solve(KKTr);

					// This gives us the Eulerian component
					vert->v = Vec3(newvE(0), newvE(1), 0.0);
					v(mesh.nodes.size() * 3 + node->EoL_index * 2) = newvE(0);
					v(mesh.nodes.size() * 3 + node->EoL_index * 2 + 1) = newvE(1);

					// We use the deformation gradient to extract the lagrangian component
					Vector3d newvL = v2e(node->v) + -deform_grad(old_face) * newvE.segment<2>(0);
					node->v = e2v(newvL);
					v(3 * n) = newvL(0);
					v(3 * n + 1) = newvL(1);
					v(3 * n + 2) = newvL(2);
				}
			}
		}
	}
}

void generateMatlab(const SparseMatrix<double>& M,
	const SparseMatrix<double>& MDK, const VectorXd& b,
	const SparseMatrix<double>& Aeq, const VectorXd& beq,
	const SparseMatrix<double>& Aineq, const VectorXd& bineq,
	VectorXd& v)
{
	mat_s2s_file(M, "M", "solver.m", false);
	mat_s2s_file(MDK, "MDK", "solver.m", false);
	vec_to_file(b, "b", "solver.m", false);
	mat_s2s_file(Aeq, "Aeq", "solver.m", false);
	vec_to_file(beq, "beq", "solver.m", false);
	mat_s2s_file(Aineq, "Aineq", "solver.m", false);
	vec_to_file(bineq, "bineq", "solver.m", false);
	vec_to_file(v, "v_input", "solver.m", false);
}

void Cloth::solve(shared_ptr<GeneralizedSolver> gs, double h)
{
	VectorXd b = -(myForces->M * v + h * myForces->f);
	generateMatlab(myForces->M,
		myForces->MDK, b,
		consts->Aeq, consts->beq,
		consts->Aineq, consts->bineq,
		v);
	bool success = gs->velocitySolve(consts->hasFixed, consts->hasCollisions,
		myForces->MDK, b,
		consts->Aeq, consts->beq,
		consts->Aineq, consts->bineq,
		v);
	vec_to_file(v, "v_solved", "solver.m", false);
}

void Cloth::step(shared_ptr<GeneralizedSolver> gs, shared_ptr<Obstacles> obs, const Vector3d& grav, double h, const bool& REMESHon, const bool& online)
{


	consts->fill(mesh, obs, fs[fsindex], h, online);
	if(REMESHon) velocityTransfer();
	myForces->fill(mesh, material, grav, h);
	double_to_file(h, "h", "solver.m", true);
	double_to_file(grav(2), "grav", "solver.m", false);
	double_to_file(material.density, "rho", "solver.m", false);
	double_to_file(material.e, "e", "solver.m", false);
	double_to_file(material.nu, "nu", "solver.m", false);
	MatrixXd x_X(mesh.nodes.size(), 5);
	VectorXi isEoL(mesh.nodes.size());
	isEoL.setZero();
	for (int i = 0; i < mesh.nodes.size(); i++) {
		if (mesh.nodes[i]->EoL) isEoL(i) = 1;
		x_X(i, 0) = mesh.nodes[i]->x[0];
		x_X(i, 1) = mesh.nodes[i]->x[1];
		x_X(i, 2) = mesh.nodes[i]->x[2];
		x_X(i, 3) = mesh.nodes[i]->verts[0]->u[0];
		x_X(i, 4) = mesh.nodes[i]->verts[0]->u[1];
	}
	mat_to_file(x_X, "x_X", "solver.m", false);
	vec_to_file(isEoL, "isEol","solver.m",false);
	MatrixXi faces2(3, mesh.faces.size());
	for (int i = 0; i < mesh.faces.size(); i++) {
		faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}
	VectorXi vvv(3);
	vvv << 1, 1, 1;
	mat_to_file(faces2.colwise() += vvv, "faces", "solver.m",false);
	vec_to_file(myForces->f, "f", "solver.m", false);
	solve(gs, h);

	for (int n = 0; n < mesh.nodes.size(); n++) {
		Node* node = mesh.nodes[n];
		Vert* vert = node->verts[0];
		node->v[0] = v(n * 3);
		node->v[1] = v(n * 3 + 1);
		node->v[2] = v(n * 3 + 2);
		node->x = node->x + h * node->v;
		if (node->EoL) {
			//cout << node->v << endl;
			vert->v[0] = v(mesh.nodes.size() * 3 + node->EoL_index * 2);
			vert->v[1] = v(mesh.nodes.size() * 3 + node->EoL_index * 2 + 1);
			//if (mesh.nodes[i]->verts[0]->u[0] != Xmin && mesh.nodes[i]->verts[0]->u[0] != Xmax) mesh.nodes[i]->verts[0]->u[0] = mesh.nodes[i]->verts[0]->u[0] + h * mesh.nodes[i]->verts[0]->v[0];
			//if (mesh.nodes[i]->verts[0]->u[1] != Ymin && mesh.nodes[i]->verts[0]->u[1] != Ymax) mesh.nodes[i]->verts[0]->u[1] = mesh.nodes[i]->verts[0]->u[1] + h * mesh.nodes[i]->verts[0]->v[1];
			vert->u[0] = vert->u[0] + h * vert->v[0];
			vert->u[1] = vert->u[1] + h * vert->v[1];
		}
	}

	updateBuffers();
}

void Cloth::updateFix(double t)
{
	if (fsindex + 1 >= fs.size()) return;
	if (t > fs[fsindex + 1]->when) {
		fsindex++;
	}
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

	myForces->drawSimple(mesh, MV, p);
}
#endif // EOLC_ONLINE

// Export
int Cloth::getBrenderCount() const
{
	return 2;
}

vector<string> Cloth::getBrenderNames() const
{
	vector<string> names;
	names.push_back("Cloth2D");
	names.push_back("Cloth3D");
	return names;
}

void Cloth::exportBrender(vector< shared_ptr< ofstream > > outfiles) const
{
	ofstream &outfile = *outfiles[0];

	for (int i = 0; i < mesh.nodes.size(); i++) {
		char vert[50];
		sprintf(vert, "v %f %f %f\n", mesh.nodes[i]->verts[0]->u[0], mesh.nodes[i]->verts[0]->u[1], 0);
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
		char norm[50];
		sprintf(norm, "vn %f %f %f\n", 0.0, 0.0, 1.0);
		outfile << norm;
	}
	//faces
	for (int i = 0; i < eleBuf.size(); i = i + 3) {
		char face[50];
		sprintf(face, "f %i/%i/%i %i/%i/%i %i/%i/%i\n", eleBuf[i] + 1, eleBuf[i] + 1, eleBuf[i] + 1, eleBuf[i + 1] + 1, eleBuf[i + 1] + 1, eleBuf[i + 1] + 1, eleBuf[i + 2] + 1, eleBuf[i + 2] + 1, eleBuf[i + 2] + 1);
		outfile << face;
	}

	ofstream &outfile2 = *outfiles[1];

	////vertex positions
	for (int i = 0; i < posBuf.size(); i = i + 3) {
		char vert[50];
		sprintf(vert, "v %f %f %f\n", posBuf[i], posBuf[i + 1], posBuf[i + 2]);
		outfile2 << vert;
	}
	//texture coordinates
	for (int i = 0; i < texBuf.size(); i = i + 2) {
		char vtex[50];
		sprintf(vtex, "vt %f %f\n", texBuf[i], texBuf[i + 1]);
		outfile2 << vtex;
	}
	//normal vectors
	for (int i = 0; i < norBuf.size(); i = i + 3) {
		char norm[50];
		sprintf(norm, "vn %f %f %f\n", norBuf[i], norBuf[i + 1], norBuf[i + 2]);
		outfile2 << norm;
	}
	//faces
	for (int i = 0; i < eleBuf.size(); i = i + 3) {
		char face[50];
		sprintf(face, "f %i/%i/%i %i/%i/%i %i/%i/%i\n", eleBuf[i] + 1, eleBuf[i] + 1, eleBuf[i] + 1, eleBuf[i + 1] + 1, eleBuf[i + 1] + 1, eleBuf[i + 1] + 1, eleBuf[i + 2] + 1, eleBuf[i + 2] + 1, eleBuf[i + 2] + 1);
		outfile2 << face;
	}
}