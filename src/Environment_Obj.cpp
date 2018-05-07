#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Environment_Obj.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"

#include "mesh.hpp"
#include "vectors.hpp"
#include "util.hpp"

using namespace std;

Env_obj::Env_obj() :
	dim(0.5, 0.5, 0.5),
	x(0.0, 0.0, 0.0)
{

}

Env_obj::Env_obj(const shared_ptr<Shape> s) :
	dim(0.5, 0.5, 0.5),
	x(0.0, 0.0, 0.0),
	env_obj(s)
{

}

Env_obj::~Env_obj()
{
}

void Env_obj::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	MV->pushMatrix();
	MV->translate(x(0), x(1), x(2));
	MV->scale(dim(0), dim(1), dim(2));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	env_obj->draw(prog);
	MV->popMatrix();
}
