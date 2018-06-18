#include "Points.h"

#ifdef EOLC_ONLINE

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

using namespace std;

void Points::drawSimple(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	for (int i = 0; i < pxyz.cols(); i++) {
		glColor3f(0.0f, 1.0f, 1.0f);
		glPointSize(10.0f);
		glBegin(GL_POINTS);
		glVertex3f(pxyz(0,i), pxyz(1, i), pxyz(2, i));
		glEnd();
	}
}

#endif // EOLC_ONLINE