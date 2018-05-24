#pragma once
#ifndef _MatrixStack_H_
#define _MatrixStack_H_

#include <stack>
#include <memory>
#include <glm/fwd.hpp>

class MatrixStack
{
public:
	MatrixStack();
	virtual ~MatrixStack();
	
	// glPushMatrix(): Copies the current matrix and adds it to the top of the stack
	void pushMatrix();
	// glPopMatrix(): Removes the top of the stack and sets the current matrix to be the matrix that is now on top
	void popMatrix();
	
	// glLoadIdentity(): Sets the top matrix to be the identity
	void loadIdentity();
	// glMultMatrix(): Right multiplies the top matrix
	void multMatrix(const glm::mat4 &matrix);
	
	// glTranslate(): Right multiplies the top matrix by a translation matrix
	void translate(const glm::vec3 &trans);
	void translate(float x, float y, float z);
	// glScale(): Right multiplies the top matrix by a scaling matrix
	void scale(const glm::vec3 &scale);
	void scale(float x, float y, float z);
	// glScale(): Right multiplies the top matrix by a scaling matrix
	void scale(float size);
	// glRotate(): Right multiplies the top matrix by a rotation matrix (angle in radians)
	void rotate(float angle, const glm::vec3 &axis);
	void rotate(float angle, float x, float y, float z);
	
	// glGet(GL_MODELVIEW_MATRIX): Gets the top matrix
	const glm::mat4 &topMatrix() const;
	
	// Prints out the specified matrix
	static void print(const glm::mat4 &mat, const char *name = 0);
	// Prints out the top matrix
	void print(const char *name = 0) const;
	
private:
	std::shared_ptr< std::stack<glm::mat4> > mstack;
	
};

#endif
