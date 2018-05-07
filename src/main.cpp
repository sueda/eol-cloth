#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

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

#include "GLSL.h"
#include "Program.h"
#include "Camera.h"
#include "MatrixStack.h"
#include "Shape.h"
#include "Scene.h"
#include "Cloth.h"
#include "PreProcess.h"
#include "boxTriCollision.h"
#include "ChronoTimer.h"

// ImGui
#include "imgui.h"
#include "imgui_impl_glfw.h"

using namespace std;
using namespace Eigen;

bool keyToggles[256] = {false}; // only for English keyboards!

GLFWwindow *window; // Main application window
string RESOURCE_DIR = ""; // Where the resources are loaded from
string JSON_FILE = ""; // The json parameter file

shared_ptr<Camera> camera;
shared_ptr<Program> prog;
shared_ptr<Program> progSimple;
shared_ptr<Scene> scene;

ChronoTimer timer("timer");

bool showCollNorms;

static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
}

static void char_callback(GLFWwindow *window, unsigned int key)
{
	keyToggles[key] = !keyToggles[key];
	switch(key) {
		case 'h':
			scene->step();
			break;
		case 'r':
			scene->reset();
			break;
	}
}

static void cursor_position_callback(GLFWwindow* window, double xmouse, double ymouse)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if(state == GLFW_PRESS) {
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
	if(action == GLFW_PRESS) {
		bool shift = mods & GLFW_MOD_SHIFT;
		bool ctrl  = mods & GLFW_MOD_CONTROL;
		bool alt   = mods & GLFW_MOD_ALT;
		camera->mouseClicked(xmouse, ymouse, shift, ctrl, alt);
	}
}

static void init()
{
	GLSL::checkVersion();
	
	// Set background color
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	// Enable z-buffer test
	glEnable(GL_DEPTH_TEST);
	// Enable alpha blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	progSimple = make_shared<Program>();
	progSimple->setShaderNames(RESOURCE_DIR + "simple_vert.glsl", RESOURCE_DIR + "simple_frag.glsl");
	progSimple->setVerbose(true); // Set this to true when debugging.
	progSimple->init();
	progSimple->addUniform("P");
	progSimple->addUniform("MV");
	//progSimple->setVerbose(false);
	
	prog = make_shared<Program>();
	prog->setVerbose(true); // Set this to true when debugging.
	prog->setShaderNames(RESOURCE_DIR + "phong_vert.glsl", RESOURCE_DIR + "phong_frag.glsl");
	prog->init();
	prog->addUniform("P");
	prog->addUniform("MV");
	prog->addUniform("kdFront");
	prog->addUniform("kdBack");
	prog->addAttribute("aPos");
	prog->addAttribute("aNor");
	//prog->setVerbose(false);
	
	camera = make_shared<Camera>();
	camera->setInitDistance(2.0f);

	scene = make_shared<Scene>();
	scene->load(RESOURCE_DIR, JSON_FILE);
	scene->tare();
	scene->init();
	
	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

void render()
{
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);
	
	// Use the window size for camera.
	glfwGetWindowSize(window, &width, &height);
	camera->setAspect((float)width/(float)height);
	
	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(keyToggles[(unsigned)'c']) {
		glEnable(GL_CULL_FACE);
	} else {
		glDisable(GL_CULL_FACE);
	}
	if(keyToggles[(unsigned)'l']) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	
	auto P = make_shared<MatrixStack>();
	auto MV = make_shared<MatrixStack>();
	
	// Apply camera transforms
	P->pushMatrix();
	camera->applyProjectionMatrix(P);
	MV->pushMatrix();
	camera->applyViewMatrix(MV);

	// Draw grid
	progSimple->bind();
	glUniformMatrix4fv(progSimple->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(progSimple->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	/*glLineWidth(2.0f);
	float x0 = -0.5f;
	float x1 = 0.5f;
	float y0 = -0.5f;
	float y1 = 0.5f;
	int gridSize = 10;
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	for(int i = 1; i < gridSize; ++i) {
		if(i == gridSize/2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		} else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float x = x0 + i / (float)gridSize * (x1 - x0);
		glVertex3f(x, y0, 0.0f);
		glVertex3f(x, y1, 0.0f);
	}
	for(int i = 1; i < gridSize; ++i) {
		if(i == gridSize/2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		} else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float y = y0 + i / (float)gridSize * (y1 - y0);
		glVertex3f(x0, y, 0.0f);
		glVertex3f(x1, y, 0.0f);
	}
	glEnd();
	glColor3f(0.4f, 0.4f, 0.4f);
	glBegin(GL_LINE_LOOP);
	glVertex3f(x0, y0, 0.0);
	glVertex3f(x1, y0, 0.0);
	glVertex3f(x1, y1, 0.0);
	glVertex3f(x0, y1, 0.0);
	glEnd();*/


	// For drawing collision normals
	if (showCollNorms) {
		if (scene->EoLon) {
			for (int i = 0; i < scene->cloth->mesh.nodes.size(); i++) {
				if (scene->cloth->mesh.nodes[i]->EoL) {
					glColor3f(0.0f, 0.0f, 0.0f);
					glPointSize(5.0f);
					//glBegin(GL_POINTS);
					//glVertex3f(scene->cloth->mesh.nodes[i]->x[0], scene->cloth->mesh.nodes[i]->x[1], scene->cloth->mesh.nodes[i]->x[2]);
					//glEnd();
					//glLineWidth(2.0f);
					//glColor3f(0.0f, 1.0f, 0.0f);
					//glBegin(GL_LINES);
					//glVertex3f(scene->cloth->mesh.nodes[i]->x[0], scene->cloth->mesh.nodes[i]->x[1], scene->cloth->mesh.nodes[i]->x[2]);
					//glVertex3f(scene->cloth->mesh.nodes[i]->x[0] + scene->cloth->mesh.nodes[i]->nor_ave[0]*.3, scene->cloth->mesh.nodes[i]->x[1] + scene->cloth->mesh.nodes[i]->nor_ave[1] * .3, scene->cloth->mesh.nodes[i]->x[2] + scene->cloth->mesh.nodes[i]->nor_ave[2] * .3);
					//glEnd();
					//glColor3f(0.0f, 0.0f, 1.0f);
					//glBegin(GL_LINES);
					//glVertex3f(scene->cloth->mesh.nodes[i]->x[0], scene->cloth->mesh.nodes[i]->x[1], scene->cloth->mesh.nodes[i]->x[2]);
					//glVertex3f(scene->cloth->mesh.nodes[i]->x[0] + scene->cloth->mesh.nodes[i]->b1[0] * .3, scene->cloth->mesh.nodes[i]->x[1] + scene->cloth->mesh.nodes[i]->b1[1] * .3, scene->cloth->mesh.nodes[i]->x[2] + scene->cloth->mesh.nodes[i]->b1[2] * .3);
					//glEnd();
					glColor3f(1.0f, 0.0f, 1.0f);
					glBegin(GL_LINES);
					glVertex3f(scene->cloth->mesh.nodes[i]->x[0], scene->cloth->mesh.nodes[i]->x[1], scene->cloth->mesh.nodes[i]->x[2]);
					glVertex3f(scene->cloth->mesh.nodes[i]->x[0] + scene->cloth->mesh.nodes[i]->v[0] * 10, scene->cloth->mesh.nodes[i]->x[1] + scene->cloth->mesh.nodes[i]->v[1] * 10, scene->cloth->mesh.nodes[i]->x[2] + scene->cloth->mesh.nodes[i]->v[2] * 10);
					glEnd();

				}
				//else {
				//	glColor3f(0.0f, 1.0f, 1.0f);
				//	glBegin(GL_LINES);
				//	glVertex3f(scene->cloth->mesh.nodes[i]->x[0], scene->cloth->mesh.nodes[i]->x[1], scene->cloth->mesh.nodes[i]->x[2]);
				//	glVertex3f(scene->cloth->mesh.nodes[i]->x[0] + scene->cloth->mesh.nodes[i]->n[0], scene->cloth->mesh.nodes[i]->x[1] + scene->cloth->mesh.nodes[i]->n[1], scene->cloth->mesh.nodes[i]->x[2] + scene->cloth->mesh.nodes[i]->n[2]);
				//	glEnd();
				//}
			}
			for (int i = 0; i < scene->cloth->collisions.size(); i++) {
				if (scene->cloth->collisions[i]->count1 == 3 && scene->cloth->mesh.nodes[scene->cloth->collisions[i]->verts2(0)]->EoL) {
					glColor3f(1.0f, 0.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex3f(scene->cloth->collisions[i]->pos2(0), scene->cloth->collisions[i]->pos2(1), scene->cloth->collisions[i]->pos2(2));
					glVertex3f(scene->cloth->collisions[i]->pos2(0) + scene->cloth->collisions[i]->nor1(0)*.3, scene->cloth->collisions[i]->pos2(1) + scene->cloth->collisions[i]->nor1(1)*.3, scene->cloth->collisions[i]->pos2(2) + scene->cloth->collisions[i]->nor1(2)*.3);
					glEnd();
					glColor3f(0.0f, 1.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex3f(scene->cloth->collisions[i]->pos2(0), scene->cloth->collisions[i]->pos2(1), scene->cloth->collisions[i]->pos2(2));
					glVertex3f(scene->cloth->collisions[i]->pos2(0) + scene->cloth->collisions[i]->tan1(0)*.3, scene->cloth->collisions[i]->pos2(1) + scene->cloth->collisions[i]->tan1(1)*.3, scene->cloth->collisions[i]->pos2(2) + scene->cloth->collisions[i]->tan1(2)*.3);
					glEnd();
					glColor3f(0.0f, 0.0f, 1.0f);
					glBegin(GL_LINES);
					glVertex3f(scene->cloth->collisions[i]->pos2(0), scene->cloth->collisions[i]->pos2(1), scene->cloth->collisions[i]->pos2(2));
					glVertex3f(scene->cloth->collisions[i]->pos2(0) + scene->cloth->collisions[i]->bin1(0)*.3, scene->cloth->collisions[i]->pos2(1) + scene->cloth->collisions[i]->bin1(1)*.3, scene->cloth->collisions[i]->pos2(2) + scene->cloth->collisions[i]->bin1(2)*.3);
					glEnd();
				}
			}
		}
		else {
			for (int i = 0; i < scene->cloth->collisions.size(); i++) {
				if (scene->cloth->collisions[i]->count1 == 3 && scene->cloth->collisions[i]->count2 == 1) {
					glColor3f(1.0f, 0.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex3f(scene->cloth->collisions[i]->pos2(0), scene->cloth->collisions[i]->pos2(1), scene->cloth->collisions[i]->pos2(2));
					glVertex3f(scene->cloth->collisions[i]->pos2(0) + scene->cloth->collisions[i]->nor1(0)*.3, scene->cloth->collisions[i]->pos2(1) + scene->cloth->collisions[i]->nor1(1)*.3, scene->cloth->collisions[i]->pos2(2) + scene->cloth->collisions[i]->nor1(2)*.3);
					glEnd();
				}
				if (scene->cloth->collisions[i]->count1 == 2 && scene->cloth->collisions[i]->count2 == 2) {
					glColor3f(0.0f, 1.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex3f(scene->cloth->collisions[i]->pos2(0), scene->cloth->collisions[i]->pos2(1), scene->cloth->collisions[i]->pos2(2));
					glVertex3f(scene->cloth->collisions[i]->pos2(0) + scene->cloth->collisions[i]->nor1(0)*.3, scene->cloth->collisions[i]->pos2(1) + scene->cloth->collisions[i]->nor1(1)*.3, scene->cloth->collisions[i]->pos2(2) + scene->cloth->collisions[i]->nor1(2)*.3);
					glEnd();
				}
				if (scene->cloth->collisions[i]->count1 == 1 && (scene->cloth->collisions[i]->count2 == 1 || scene->cloth->collisions[i]->count2 == 2 || scene->cloth->collisions[i]->count2 == 3)) {
					glColor3f(0.0f, 0.0f, 1.0f);
					glBegin(GL_LINES);
					glVertex3f(scene->cloth->collisions[i]->pos1_(0), scene->cloth->collisions[i]->pos1_(1), scene->cloth->collisions[i]->pos1_(2));
					glVertex3f(scene->cloth->collisions[i]->pos1_(0) + scene->cloth->collisions[i]->nor1(0)*.3, scene->cloth->collisions[i]->pos1_(1) + scene->cloth->collisions[i]->nor1(1)*.3, scene->cloth->collisions[i]->pos1_(2) + scene->cloth->collisions[i]->nor1(2)*.3);
					glEnd();
				}
				glColor3f(1.0f, 0.0f, 1.0f);
				glBegin(GL_LINES);
				glVertex3f(scene->cloth->collisions[i]->pos1_(0), scene->cloth->collisions[i]->pos1_(1), scene->cloth->collisions[i]->pos1_(2));
				glVertex3f(scene->cloth->collisions[i]->pos1_(0) + scene->cloth->mesh.nodes[scene->cloth->collisions[i]->verts1(0)]->v[0], scene->cloth->collisions[i]->pos1_(1) + scene->cloth->mesh.nodes[scene->cloth->collisions[i]->verts1(0)]->v[1], scene->cloth->collisions[i]->pos1_(2) + scene->cloth->mesh.nodes[scene->cloth->collisions[i]->verts1(0)]->v[2]);
				glEnd();
			}
		}
	}

	progSimple->unbind();

	// Draw scene
	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	scene->draw(MV, prog);
	MV->popMatrix();
	prog->unbind();
	
	//////////////////////////////////////////////////////
	// Cleanup
	//////////////////////////////////////////////////////
	
	// Pop stacks
	MV->popMatrix();
	P->popMatrix();
	
	GLSL::checkError(GET_FILE_LINE);
}

void stepperFunc()
{
	while(true) {
		if(keyToggles[(unsigned)' ']) {
			//timer.tic();
			scene->step();
			//timer.toc();
			//timer.print();
		}
		this_thread::sleep_for(chrono::microseconds(1));
	}
}

int main(int argc, char **argv)
{
	if(argc < 3) {
		cout << "Please specify the resource directory and json file." << endl;
		return 0;
	}
	RESOURCE_DIR = argv[1] + string("/");
	JSON_FILE = argv[1] + string("/") + argv[2];
	
	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if(!glfwInit()) {
		return -1;
	}
	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(1280, 720, "Nick Weidner", NULL, NULL);
	if(!window) {
		glfwTerminate();
		return -1;
	}
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if(glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}
	// ImGui
	ImGui_ImplGlfw_Init(window, true);
	ImVec4 clear_color = ImColor(114, 144, 154);
	showCollNorms = false;

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
	// Initialize scene.
	init();
	// Start simulation thread.
	thread stepperThread(stepperFunc);
	// Loop until the user closes the window.
	while(!glfwWindowShouldClose(window)) {
		ImGui_ImplGlfw_NewFrame();
		// 1. Show a simple window
		// Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets appears in a window automatically called "Debug"
		{
			ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Always);
			ImGui::Begin("Stats", NULL, ImGuiWindowFlags_NoMove);
			static float f = 0.0f;
			ImGui::Text("Vert-Face: Blue");
			ImGui::Text("Edge-Edge: Green");
			ImGui::Text("Face-Vert: Red");
			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
			ImGui::Text("Node resolution: %d number of points", scene->cloth->mesh.nodes.size());
			ImGui::Text("Face resolution: %d number of triangles", scene->cloth->mesh.faces.size());
			ImGui::Checkbox("Show collision normals", &showCollNorms);
			ImGui::Checkbox("Output Matlab physics debug file", &scene->cloth->matlab_debug_physics);
			ImGui::Checkbox("Output Matlab collisoin debug file", &scene->pps->matlab_debug_collision);
			ImGui::Checkbox("Output 2D postscript", &scene->pps->export_postscript);
			ImGui::End();
		}
		
		// Render scene.
		render();
		ImGui::Render();
		// Swap front and back buffers.
		glfwSwapBuffers(window);
		// Poll for and process events.
		glfwPollEvents();

		

	}
	// Quit program.
	stepperThread.detach();
	glfwDestroyWindow(window);
	ImGui_ImplGlfw_Shutdown();
	glfwTerminate();
	return 0;
}
