#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "glm/ext.hpp"

#include <math.h>

#include "Box.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Rigid.h"

// ArcSim
#include "mesh.hpp"
#include "vectors.hpp"
#include "util.hpp"
#include "io.hpp"

using namespace std;

Box::Box() :
	dim(1.0, 1.0, 1.0),
	x(0.0, 0.0, 0.0)
{
	//rigid = make_shared<Rigid>();
	v.resize(6);
}

Box::Box(const shared_ptr<Shape> s) :
	dim(1.0, 1.0, 1.0),
	x(0.0, 0.0, 0.0),
	E1(Eigen::Matrix4d::Identity()),
	E1inv(Eigen::Matrix4d::Identity()),
	adjoint(Eigen::Matrix<double, 6, 6>::Identity()),
	box(s)
{
	//rigid = make_shared<Rigid>();
	v.resize(6);
}

Box::~Box()
{
}

void Box::init()
{
	box->init();
}

void Box::step(double h)
{
	//adjoint = rigid->adjoint(E1inv);
	adjoint = Matrix6d::Identity();
	adjoint.block<3, 3>(0, 0) = E1.block<3, 3>(0, 0).transpose();
	adjoint.block<3, 3>(3, 3) = E1.block<3, 3>(0, 0).transpose();
	E1 = rigid->integrate(E1, adjoint*v, h);
	E1inv = E1.inverse();
}

void Box::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	MV->pushMatrix();
	//MV->print();
	//MV->translate(x(0), x(1), x(2));
	MV->translate(E1(0,3), E1(1,3), E1(2,3));
	//MV->scale(dim(0), dim(1), dim(2));
	//angleaxisd .angle
	//MV->rotate(rot(0), glm::vec3(1, 0, 0));
	//MV->rotate(rot(1), glm::vec3(0, 1, 0));
	//MV->rotate(rot(2), glm::vec3(0, 0, 1));
	//Eigen::Matrix<double, 4, 4, Eigen::ColMajor> E1row = E1;
	//E1row.row(1).swap(E1row.row(2));
	//glm::mat4 bbb = glm::make_mat4(E1row.data());

	Eigen::Matrix4d E1row = Eigen::Matrix4d::Identity();
	E1row.block<3, 3>(0, 0) = E1.block<3, 3>(0, 0);
	//E1row.row(1).swap(E1row.row(2));
	glm::mat4 bbb = glm::make_mat4(E1row.data());

	//cout << E1row << endl;
	//MV->print();
	MV->multMatrix(bbb);
	//MV->print();

	//MV->translate(E1(0, 3), E1(1, 3), E1(2, 3));
	MV->scale(dim(0), dim(1), dim(2));

	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	box->draw(prog);
	MV->popMatrix();
}

// Export
int Box::getBrenderCount() const
{
	return 1;
}

vector<string> Box::getBrenderNames() const
{
	vector<string> names;
	names.push_back("Box");
	return names;
}

void Box::exportBrender(vector< shared_ptr< ofstream > > outfiles) const
{
	ofstream &outfile = *outfiles[0];

	Eigen::Matrix4d T, S;
	S.setIdentity();
	S(0, 0) = dim(0);
	S(1, 1) = dim(1);
	S(2, 2) = dim(2);
	//T.setIdentity();
	//T.block<3, 1>(0, 3) = x;

	box->exportBrender(E1*S, outfile);


}