#include <iostream>

#ifdef EOLC_ONLINE
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "glm/ext.hpp"

#include "online/Program.h"
#include "online/MatrixStack.h"
#endif // EOLC_ONLINE

#include <math.h>

#include "Box.h"
#include "Shape.h"
#include "Rigid.h"

using namespace std;

Box::Box(const shared_ptr<Shape> s) :
	dim(1.0, 1.0, 1.0),
	x(0.0, 0.0, 0.0),
	E1(Eigen::Matrix4d::Identity()),
	E1inv(Eigen::Matrix4d::Identity()),
	adjoint(Eigen::Matrix<double, 6, 6>::Identity()),
	v(Eigen::Matrix<double, 6, 1>::Zero()),
	boxShape(s)
{

}

Box::~Box()
{
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

#ifdef EOLC_ONLINE

void Box::init()
{
	boxShape->init();
}

void Box::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	MV->pushMatrix();
	//MV->translate(E1(0,3), E1(1,3), E1(2,3));

	Eigen::Matrix4d E1row = Eigen::Matrix4d::Identity();
	E1row.block<3, 3>(0, 0) = E1.block<3, 3>(0, 0);
	glm::mat4 bbb = glm::make_mat4(E1.data());
	MV->multMatrix(bbb);

	MV->scale(dim(0), dim(1), dim(2));

	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	boxShape->draw(prog);
	MV->popMatrix();
}

#endif // EOLC_ONLINE

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

	boxShape->exportBrender(E1*S, outfile);


}