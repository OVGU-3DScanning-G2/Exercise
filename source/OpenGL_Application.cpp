// OpenGL_Application.cpp : Defines the entry point for the console application.
//

#define USE_GLFW

#include <string>     //we want to process text + filenames
#include <iostream>   //for making output text to the console with "cout"
#include <vector>     //for std::vector class functions
#include <algorithm>  //for std::sort
#include <stdio.h>
#include <ctime>

#define GLFW_INCLUDE_GLU
#include "GLFW/glfw3.h" //inlcude the function definition

#include <GL/glu.h>

#include "GLcamera.h"
#include "Point3d.h"
#include "KDTree.h"

//Normally compiler & linker options are set in the project file and not in the source code
//Its just here to show what dependencies are needed
#pragma comment(lib, "opengl32.lib")     //link against standard MS Windows OpenGL-Library
#pragma comment(lib, "glu32.lib")        //link to some standard OpenGL convenience functions
#pragma comment(lib, "GLFW/glfw3.lib")   //link against the the GLFW OpenGL SDK

//using namespace std; //everything what is in the "Standard" C++ namespace, so the "std::" prefix can be avoided

void updateScene(const std::vector<Point3d>& points);
void loadFileXYZ(const char* filename, std::vector<Point3d>& points); //declare our own read-function, implementation is done below the main(...)-function

void drawPoints(std::vector<Point3d>& points, double size, GLbyte color_r, GLbyte color_g, GLbyte color_b);
void drawBox();
void drawCircle();
void drawCoordinateAxes();
void drawBackground();

//Variables to control our virtual camera movements
GLcamera m_camera;
Point3d  m_sceneCenter;   //scene center and rotation point
double   m_sceneRadius = 0; //scene radius
Point3d  m_bbmin, m_bbmax;//bounding box
//Own Variables
//-----------------------------------------
std::vector<Point3d> res;
std::vector<Point3d> points;
std::vector<Point3d> abfrage;
KDTree data;
double abfrageLaenge;
int startDim = 0;
int pointSize = 2;
//-----------------------------------------

int m_windowWidth = 0;
int m_windowHeight = 0;

double lastposX, lastposY; //last mouse positions

void initializeGL()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

//This function is called, wenn the framebuffer size (e.g. if the window size changes)
void resizeGL(GLFWwindow* window, int width, int height)
{
	m_windowWidth = width;
	m_windowHeight = height;

	glViewport(0, 0, width, height);  //set our viewport

	m_camera.setWindowSize(width, height);
	m_camera.updateProjection(); //adjust projection to new window size
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
	{
		//KDTree - Abfrage
		//----------------------------------------------------------------------------
		abfrage.clear();
		int random = (std::rand() % (points.size()));
		abfrage.emplace_back(points[random]);

		Point3d S = m_bbmax - m_bbmin;
		abfrageLaenge = S.x * 0.25;

		clock_t begin = clock();
		res = data.getRange(abfrageLaenge, abfrage[0], startDim);
		clock_t end = clock();

		std::cout << "Time needed to calculate range query: " << double(end - begin) / CLOCKS_PER_SEC << "s\r";
		//----------------------------------------------------------------------------
	}

	if (key == GLFW_KEY_N && action == GLFW_RELEASE)
	{
		abfrage.clear();
		int random = (std::rand() % (points.size()));
		abfrage.emplace_back(points[random]);

		res.clear();

		clock_t begin = clock();
		res.emplace_back(data.getNN(abfrage[0]));
		clock_t end = clock();

		std::cout << "Time needed to calculate NN: " << double(end - begin) / CLOCKS_PER_SEC << "s\r";
	}
}

//This function is called, wenn a mouse button is pressed
//We use this function to get the initial mouse position for our camera rotation
void mouse_press_event(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) //wenn the left mouse button is pressed
	{
		glfwGetCursorPos(window, &lastposX, &lastposY); //get the current mouse coursor position
	}
}

//This function is called when the mouse is moved
//The mouse pos is (x,y) and we take the sphere equation (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2
//we solve for z to get the 3rd coordinate
//The sphere center at(x0,y0,z0) and x0=windowWidth/2, y0=windowHeight/2, z0=0;
//This results in z= sqrt ( r^2 - (x-x0)^2 - (y-y0)^2 )
void mouse_move_event(GLFWwindow* window, double currposX, double currposY)
{
	if (lastposX == currposX && lastposY == currposY) return;

	//check the mouse button state
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	//compute object rotation wenn the left mouse button is pressed
	if (state == GLFW_PRESS) //if the left mouse button is pressed while it is moved then we rotate
	{
		m_camera.rotate((int)lastposX, (int)lastposY, (int)currposX, (int)currposY);

		//the current mouse position is the last mouse position for the next interaction
		lastposX = currposX;
		lastposY = currposY;
	}
}

//The callback function receives two - dimensional scroll offsets.
void mouse_wheel_event(GLFWwindow* window, double xoffset, double yoffset)
{
	const double factor = (yoffset < 0) ? 1.1 : 0.9;
	m_camera.zoom(factor);
}

int main(int argc, char* argv[]) //this function is called, wenn ou double-click an ".EXE" file
{
	std::cout << "Hello world! \n This is my console window" << std::endl;

	//try to load point cloud data from file
	clock_t begin = clock();

	//loadFileXYZ("data/Stanford Dragon.xyz", points);
	loadFileXYZ("data/cone.xyz", points);

	clock_t end = clock();
	std::cout << "Time needed to load data: " << double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;

	//OK, we now compute the min and max coordinates for our bounding box
	updateScene(points);

	//Load KD-Tree
	//----------------------------------------------------------------------------
	begin = clock();

	data = KDTree(points, startDim);

	end = clock();
	std::cout << "Time needed to build kdTree: " << double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;
	//----------------------------------------------------------------------------

	//Create an OpenGL window with GLFW
	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		getc(stdin);
		return -1;
	}

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(1280, 720, "My OpenGL Window", NULL, NULL);
	if (!window)
	{
		fprintf(stderr, "Failed glfwCreateWindow\n");
		glfwTerminate();
		return -1;
	}

	//setting up events/callbacks
	glfwSetKeyCallback(window, key_callback);
	glfwSetFramebufferSizeCallback(window, resizeGL);
	glfwSetCursorPosCallback(window, mouse_move_event);
	glfwSetMouseButtonCallback(window, mouse_press_event);
	glfwSetScrollCallback(window, mouse_wheel_event);

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	//Prepare our virtual camera
	glfwGetFramebufferSize(window, &m_windowWidth, &m_windowHeight);

	//Initialize Camera
	m_camera.setWindowSize(m_windowWidth, m_windowHeight);    //setup window parameters
	m_camera.initializeCamera(m_sceneCenter, m_sceneRadius);  //set the camera outside the scene
	m_camera.updateProjection();                              //adjust projection to window size

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		/* Render here */
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //clear buffers
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);               //clear background color
		glClearDepth(1.0f);                                 //clear depth buffer

		//draws the scene background
		drawBackground();

		//draw points
		//----------------------------------------------------------------------------
		drawPoints(points, pointSize, 0, 255, 0);
		//----------------------------------------------------------------------------

		//draw abfrage
		//----------------------------------------------------------------------------
		drawPoints(abfrage, pointSize + 10, 255, 165, 0);
		//----------------------------------------------------------------------------

		//draw res-points
		//----------------------------------------------------------------------------
		drawPoints(res, pointSize + 10, 255, 0, 0);
		//----------------------------------------------------------------------------

		//draw bounding box for Abfrage
		//----------------------------------------------------------------------------
		if (!res.empty())
		{
			glPushMatrix();
			glPushAttrib(GL_POLYGON_BIT);
			glColor3ub(255, 255, 255);
			glTranslated(abfrage[0].x - abfrageLaenge, abfrage[0].y - abfrageLaenge, abfrage[0].z - abfrageLaenge);
			glScaled(abfrageLaenge * 2, abfrageLaenge * 2, abfrageLaenge * 2);

			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); //draw wire frame instead of filled quads
			drawBox();
			glPopAttrib();
			glPopMatrix();
		}
		//----------------------------------------------------------------------------

		//draw coordinate axes
		drawCoordinateAxes();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();

	return 0;
}

void drawPoints(std::vector<Point3d>& points, double size, GLbyte color_r, GLbyte color_g, GLbyte color_b)
{
	glPointSize(size);

	glEnable(GL_DEPTH_TEST);

	if (!points.empty())
	{ /* Drawing Points with VertexArrays */
		glColor3ub(color_r, color_g, color_b);
		glEnableClientState(GL_VERTEX_ARRAY); //enable data upload to GPU
		glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &points[0]);

		//draw point cloud
		glDrawArrays(GL_POINTS, 0, (unsigned int)points.size());
		glDisableClientState(GL_VERTEX_ARRAY);  //disable data upload to GPU
	}
}

//computes the bounding box of a 3d point cloud
void updateScene(const std::vector<Point3d>& points)
{
	//if there are no points we return an empty bounding box
	if (points.empty())
	{
		m_bbmin.x = 0;  m_bbmin.y = 0;  m_bbmin.z = 0;
		m_bbmax.x = 0;  m_bbmax.y = 0;  m_bbmax.z = 0;
		return;
	}

	//We now compute the min and max coordinates for our bounding box
	m_bbmin = points.front(); //initialize min with the first point
	m_bbmax = points.front(); //initialize max with the first point

	for (unsigned int i = 0; i < points.size(); ++i)
	{
		const Point3d& pt = points[i]; //do not copy but get a reference to the i-th point in the vector
		if (pt.x < m_bbmin.x) m_bbmin.x = pt.x;
		else if (pt.x > m_bbmax.x) m_bbmax.x = pt.x;

		if (pt.y < m_bbmin.y) m_bbmin.y = pt.y;
		else if (pt.y > m_bbmax.y) m_bbmax.y = pt.y;

		if (pt.z < m_bbmin.z) m_bbmin.z = pt.z;
		else if (pt.z > m_bbmax.z) m_bbmax.z = pt.z;
	}

	//check how many points we have read from file
	std::cout << "point vector now contains: " << points.size() << " points" << std::endl;

	if (points.empty())
	{
		std::cout << "ERROR: no points to show...(press enter to exit)" << std::endl;
		getc(stdin);
		return;
	}

	m_sceneCenter = (m_bbmax + m_bbmin) * 0.5;
	m_sceneRadius = distance3d(m_sceneCenter, m_bbmax);

	std::cout << "\nBounding Box was computed:\n";
	std::cout << "minPoint is: " << m_bbmin.x << "," << m_bbmin.y << "," << m_bbmin.z << std::endl;
	std::cout << "maxPoint is: " << m_bbmax.x << "," << m_bbmax.y << "," << m_bbmax.z << std::endl;
}

//Here is the implementation of our file reader
void loadFileXYZ(const char* filename, std::vector<Point3d>& points)
{
	FILE* file = 0;
	int error = fopen_s(&file, filename, "rt"); //r= read, t=text
	if (error != 0)
	{
		std::cout << "file " << filename << " could not be opened!" << std::endl;
		return; //nothing can be done else -> end function
	}

	std::cout << "reading file: " << filename << std::endl;

	while (!feof(file)) //as long we have not reached the end-of-file
	{
		Point3d point;
		int items = fscanf_s(file, "%lf %lf %lf\n", &point.x, &point.y, &point.z);

		if (items != 3) //we ecpected that 3 values have been read (except we are already at the end of file)
		{
			std::cout << "file format error" << std::endl;
			break; //abort while loop
		}
		else
		{
			points.emplace_back(point); //add the current point to our point vector
		}
	}

	//dont forget to close to file
	fclose(file);

	unsigned int numberOfPoints = points.size();

	std::cout << "reading finished: " << numberOfPoints << " points have be read" << std::endl;
}

/** draws a unit box.
origin is at (0,0,0) and maximum at (1,1,1).
use glTranslate to move the origin to the minimum coordinates
afterwards use glScale(sx,sy,sz) resize the box.
*/
void drawBox()
{
	glBegin(GL_QUADS);
	glVertex3d(0, 0, 0); glVertex3d(1, 0, 0); glVertex3d(1, 1, 0); glVertex3d(0, 1, 0); //front
	glVertex3d(0, 0, 1); glVertex3d(1, 0, 1); glVertex3d(1, 1, 1); glVertex3d(0, 1, 1); //back
	glVertex3d(0, 0, 0); glVertex3d(0, 0, 1); glVertex3d(0, 1, 1); glVertex3d(0, 1, 0); //left
	glVertex3d(1, 0, 0); glVertex3d(1, 0, 1); glVertex3d(1, 1, 1); glVertex3d(1, 1, 0); //right
	glVertex3d(0, 0, 0); glVertex3d(1, 0, 0); glVertex3d(1, 0, 1); glVertex3d(0, 0, 1); //bottom
	glVertex3d(0, 1, 0); glVertex3d(1, 1, 0); glVertex3d(1, 1, 1); glVertex3d(0, 1, 1); //top
	glEnd();
}

//draws a unit circle with radius 1 at the origin(0,0,0)
//use glTranslate, glScale and glRotate to resize and reposition the circle
void drawCircle()
{
	const int segments = 180;

	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < segments; ++i)
	{
		const double theta = 2.0 * 3.1415926 * double(i) / double(segments);
		glVertex2d(cos(theta), sin(theta));
	}
	glEnd();
}

//draws the coordinate system
void drawCoordinateAxes()
{
	//draw coordinate frame
	glBegin(GL_LINES);
	//draw line for X-Axis
	glColor3ub(255, 0, 0);
	glVertex3d(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
	glVertex3d(m_sceneCenter.x + m_sceneRadius, m_sceneCenter.y, m_sceneCenter.z);
	//draw line for Y-Axis
	glColor3ub(0, 255, 0);
	glVertex3d(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
	glVertex3d(m_sceneCenter.x, m_sceneCenter.y + m_sceneRadius, m_sceneCenter.z);
	//draw line for Z-Axis
	glColor3ub(0, 0, 255);
	glVertex3d(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
	glVertex3d(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z + m_sceneRadius);
	glEnd();

	//draw center point as a sphere
	glPushMatrix();
	glColor3ub(255, 255, 0);
	glTranslated(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
	GLUquadric* quad = gluNewQuadric();
	gluSphere(quad, 0.01, 30, 30);
	gluDeleteQuadric(quad);
	glPopMatrix();

	glPushMatrix();
	glColor3ub(64, 164, 164);
	glTranslated(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
	glScaled(m_sceneRadius, m_sceneRadius, m_sceneRadius);
	drawCircle();
	//draw another circle 90 degree rotated
	glRotated(90, 1, 0, 0);
	drawCircle();
	glPopMatrix();

	//draw bounding box
	glPushMatrix();
	glPushAttrib(GL_POLYGON_BIT);
	glColor3ub(255, 255, 255);
	Point3d S = m_bbmax - m_bbmin;
	glTranslated(m_bbmin.x, m_bbmin.y, m_bbmin.z);
	glScaled(S.x, S.y, S.z);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); //draw wire frame instead of filled quads
	drawBox();
	glPopAttrib();
	glPopMatrix();
}

//------------------------------------------------------------------------------
/** draws a static permanent 2d color gradient.
Drawing static background means that we draw a 2D rectangle, which is not
influenced by scene rotation and by the camera projection.
Basically we just want to draw a 2D texture/image.
*/
//------------------------------------------------------------------------------
void drawBackground()
{
	const float winWidth((float)m_windowWidth);
	const float winHeight((float)m_windowHeight);

	glPushAttrib(GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT | GL_CURRENT_BIT);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_PROJECTION);           //At first we select the Projection matrix
	glPushMatrix();                        //save Projektion matrix to restore it at the end
	glLoadIdentity();                      //initialize projection matrix
	gluOrtho2D(0, winWidth, 0, winHeight); //select orthographic projecton

	glMatrixMode(GL_MODELVIEW);            //now change to Modelview because we want to draw something
	glPushMatrix();
	glLoadIdentity();

	glBegin(GL_QUADS);
	glColor3ub(100, 100, 100); //color bottom
	glVertex2f(0.0f, 0.0f);  glVertex2f(winWidth, 0.0f);
	glColor3ub(38, 38, 38);  //color top
	glVertex2f(winWidth, winHeight);  glVertex2f(0.0f, winHeight);
	glEnd();

	glMatrixMode(GL_PROJECTION); //select to Projektionmatrix
	glPopMatrix();  //reset 2D-projection

	glMatrixMode(GL_MODELVIEW); //select MODELVIEW
	glPopMatrix();  //reset ModelviewMatrix to last state

	glPopAttrib();  //restore last attributes
}
