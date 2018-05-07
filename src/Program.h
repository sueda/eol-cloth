#pragma  once
#ifndef __Program__
#define __Program__

#include <map>
#include <string>

#define GLEW_STATIC
#include <GL/glew.h>

/**
 * An OpenGL Program (vertex and fragment shaders)
 */
class Program
{
public:
	Program();
	virtual ~Program();
	
	void setVerbose(bool v) { verbose = v; }
	bool isVerbose() const { return verbose; }
	
	void setShaderNames(const std::string &v, const std::string &f);
	virtual bool init();
	virtual void bind();
	virtual void unbind();

	void addAttribute(const std::string &name);
	void addUniform(const std::string &name);
	GLint getAttribute(const std::string &name) const;
	GLint getUniform(const std::string &name) const;
	
protected:
	std::string vShaderName;
	std::string fShaderName;
	
private:
	GLuint pid;
	std::map<std::string,GLint> attributes;
	std::map<std::string,GLint> uniforms;
	bool verbose;
};

#endif
