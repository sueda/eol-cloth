#include "runner.h"

#ifdef EOLC_ONLINE

#ifndef _GLIBCXX_USE_NANOSLEEP
#define _GLIBCXX_USE_NANOSLEEP
#endif
#include <thread>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "online\Camera.h"
#include "online\Program.h"
#include "online\MatrixStack.h"
#include "online\GLSL.h"

#endif // EOL_ONLINE

#include <Eigen\Core>

#include <memory>
#include <iostream>

#include "parseParams.h";
#include "genSet.h";
#include "Scene.h"

using namespace std;
using namespace Eigen;

shared_ptr<genSet> gs;
shared_ptr<Scene> scene;

#ifdef EOLC_ONLINE
bool keyToggles[256] = { false }; // only for English keyboards!

GLFWwindow *window; // Main application window

shared_ptr<Camera> camera;
shared_ptr<Program> progPhong;
shared_ptr<Program> progSimple;
#endif // EOLC_ONLINE

void init_offline(const shared_ptr<genSet> gs, const string &SIMSET_FILE)
{

}

void run_offline()
{
	while (true) {
		//step
	}
}

#ifdef EOLC_ONLINE

static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
}

static void char_callback(GLFWwindow *window, unsigned int key)
{
	keyToggles[key] = !keyToggles[key];
	switch (key) {
	case 'h':
		scene->step();
		break;
	case 'r':
		//scene->reset();
		break;
	case 'p':
		scene->partialStep();
		break;
	}
}

static void cursor_position_callback(GLFWwindow* window, double xmouse, double ymouse)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if (state == GLFW_PRESS) {
		camera->mouseMoved(xmouse, ymouse);
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the current mouse position.
	double xmouse, ymouse;
	glfwGetCursorPos(window, &xmouse, &ymouse);
	// Get current window size.
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	if (action == GLFW_PRESS) {
		bool shift = mods & GLFW_MOD_SHIFT;
		bool ctrl = mods & GLFW_MOD_CONTROL;
		bool alt = mods & GLFW_MOD_ALT;
		camera->mouseClicked(xmouse, ymouse, shift, ctrl, alt);
	}
}

void stepperFunc()
{
	while (true) {
		if (keyToggles[(unsigned)' ']) {
			scene->step();
		}
		this_thread::sleep_for(chrono::microseconds(1));
	}
}

bool init_online(const shared_ptr<genSet> gs, const string &SIMSET_FILE)
{
	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if (!glfwInit()) {
		return false;
	}
	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(1280, 720, "EOL Cloth", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return false;
	}
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return false;
	}

	glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	// Set vsync.
	glfwSwapInterval(1);
	// Set keyboard callback.
	glfwSetKeyCallback(window, key_callback);
	// Set char callback.
	glfwSetCharCallback(window, char_callback);
	// Set cursor position callback.
	glfwSetCursorPosCallback(window, cursor_position_callback);
	// Set mouse button callback.
	glfwSetMouseButtonCallback(window, mouse_button_callback);

	GLSL::checkVersion();

	// Set background color
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	// Enable z-buffer test
	glEnable(GL_DEPTH_TEST);
	// Enable alpha blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	progSimple = make_shared<Program>();
	progSimple->setShaderNames(gs->RESOURCE_DIR + "simple_vert.glsl", gs->RESOURCE_DIR + "simple_frag.glsl");
	progSimple->setVerbose(true); // Set this to true when debugging.
	progSimple->init();
	progSimple->addUniform("P");
	progSimple->addUniform("MV");
	//progSimple->setVerbose(false);

	progPhong = make_shared<Program>();
	progPhong->setVerbose(true); // Set this to true when debugging.
	progPhong->setShaderNames(gs->RESOURCE_DIR + "phong_vert.glsl", gs->RESOURCE_DIR + "phong_frag.glsl");
	progPhong->init();
	progPhong->addUniform("P");
	progPhong->addUniform("MV");
	progPhong->addUniform("kdFront");
	progPhong->addUniform("kdBack");
	progPhong->addAttribute("aPos");
	progPhong->addAttribute("aNor");
	//prog->setVerbose(false);

	camera = make_shared<Camera>();
	camera->setInitDistance(2.0f);

	scene = make_shared<Scene>();
	scene->load(gs->RESOURCE_DIR);
	load_simset(scene, SIMSET_FILE);
	scene->init(gs->online, gs->exportObjs);

	// If there were any OpenGL errors, this will print something.
	GLSL::checkError(GET_FILE_LINE);

	return true;
}

void render()
{
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);

	// Use the window size for camera.
	glfwGetWindowSize(window, &width, &height);
	camera->setAspect((float)width / (float)height);

	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if (keyToggles[(unsigned)'c']) {
		glEnable(GL_CULL_FACE);
	}
	else {
		glDisable(GL_CULL_FACE);
	}
	if (keyToggles[(unsigned)'l']) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	auto P = make_shared<MatrixStack>();
	auto MV = make_shared<MatrixStack>();

	// Apply camera transforms
	P->pushMatrix();
	camera->applyProjectionMatrix(P);
	MV->pushMatrix();
	camera->applyViewMatrix(MV);

	progSimple->bind();
	glUniformMatrix4fv(progSimple->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(progSimple->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	scene->drawSimple(MV, progSimple);

	progSimple->unbind();

	progPhong->bind();
	glUniformMatrix4fv(progPhong->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	scene->draw(MV, progPhong);
	MV->popMatrix();
	progPhong->unbind();

	MV->popMatrix();
	P->popMatrix();

	GLSL::checkError(GET_FILE_LINE);
}

void run_online()
{
	// Start simulation thread.
	thread stepperThread(stepperFunc);

	while (!glfwWindowShouldClose(window)) {
		render();
		// Swap front and back buffers.
		glfwSwapBuffers(window);
		// Poll for and process events.
		glfwPollEvents();
	}
	stepperThread.detach();
	glfwDestroyWindow(window);
	glfwTerminate();
}

#endif // EOL_ONLINE

void start_running(const string &GENSET_FILE, const string &SIMSET_FILE)
{
	gs = make_shared<genSet>();
	load_genset(gs, GENSET_FILE);

	if (gs->online) {
#ifdef EOLC_ONLINE
		init_online(gs, SIMSET_FILE);
		run_online();
#else
		cout << "ERROR: Attempting to run in online mode without building the online libraries." << endl;
#endif // EOL_ONLINE
	}
	else {
		init_offline(gs, SIMSET_FILE);
		run_offline();
	}
}


