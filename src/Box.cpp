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

double faceNorm_data[] = {
	1.0,	0.0,	0.0,
	-1.0,	0.0,	0.0,
	0.0,	1.0,	0.0,
	0.0,	-1.0,	0.0,
	0.0,	0.0,	1.0,
	0.0,	0.0,	-1.0
};

int edgeFaces_data[] = {
	1,3,
	1,5,
	1,4,
	1,2,
	3,5,
	3,4,
	0,3,
	2,5,
	0,5,
	2,4,
	0,4,
	0,2
};

int edgeTan_data[] = {
	4,
	2,
	2,
	4,
	0,
	0,
	4,
	0,
	2,
	0,
	2,
	4
};

Box::Box(const shared_ptr<Shape> s) :
	num_points(8),
	num_edges(12),
	dim(1.0, 1.0, 1.0),
	x(0.0, 0.0, 0.0),
	E1(Eigen::Matrix4d::Identity()),
	E1inv(Eigen::Matrix4d::Identity()),
	adjoint(Eigen::Matrix<double, 6, 6>::Identity()),
	v(Eigen::Matrix<double, 6, 1>::Zero()),
	boxShape(s)
{
	Eigen::Map<Eigen::Matrix<double, 3, 6, Eigen::ColMajor> > faceNorms_(faceNorm_data);
	faceNorms = E1.block<3, 3>(0, 0) * faceNorms_;
	Eigen::Map<Eigen::Matrix<int, 2, 12, Eigen::ColMajor> > edgeFaces1_(edgeFaces_data);
	edgeFaces = edgeFaces1_;
	Eigen::Map<Eigen::Matrix<int, 1, 12> > edgeTan_(edgeTan_data);
	edgeTan = edgeTan_;
}

Box::~Box()
{
}



void Box::step(const double h)
{
	//adjoint = rigid->adjoint(E1inv);
	adjoint = Matrix6d::Identity();
	adjoint.block<3, 3>(0, 0) = E1.block<3, 3>(0, 0).transpose();
	adjoint.block<3, 3>(3, 3) = E1.block<3, 3>(0, 0).transpose();
	E1 = rigid->integrate(E1, adjoint*v, h);
	E1inv = E1.inverse();

	Eigen::Map<Eigen::Matrix<double, 3, 6, Eigen::ColMajor> > faceNorms_(faceNorm_data);
	faceNorms = E1.block<3, 3>(0, 0) * faceNorms_;
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

void Box::drawSimple(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	for (int i = 0; i < faceNorms.cols(); i++) {
		glColor3f(1.0f, 0.0f, 1.0f);
		glBegin(GL_LINES);
		glVertex3f(E1(0,3), E1(1, 3), E1(2, 3));
		glVertex3f(E1(0, 3) + faceNorms(0,i), E1(1, 3) + faceNorms(1, i), E1(2, 3) + faceNorms(2, i));
		glEnd();
	}
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