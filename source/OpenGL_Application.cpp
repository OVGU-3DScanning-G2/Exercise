// OpenGL_Application.cpp : Defines the entry point for the console application.
//

#define USE_GLFW

#include <string>     //we want to process text + filenames
#include <iostream>   //for making output text to the console with "cout"
#include <vector>     //for std::vector class functions
#include <algorithm>  //for std::sort
#include <stdio.h>
#include <ctime>
#ifdef __linux__
#include <fstream>
#endif

#define GLFW_INCLUDE_GLU
#include "GLFW/glfw3.h" //inlcude the function definition

#include <GL/glu.h>

#include "../include/GLcamera.h"
#include "../include/Point3d.h"
#include "../include/KDTree.h"
#include "../include/Algorithms.h"

//Normally compiler & linker options are set in the project file and not in the source code
//Its just here to show what dependencies are needed
#pragma comment(lib, "opengl32.lib")     //link against standard MS Windows OpenGL-Library
#pragma comment(lib, "glu32.lib")        //link to some standard OpenGL convenience functions
#pragma comment(lib, "GLFW/glfw3.lib")   //link against the the GLFW OpenGL SDK


std::string filename = "../data/cone.xyz";

//using namespace std; //everything what is in the "Standard" C++ namespace, so the "std::" prefix can be avoided

void updateScene(const std::vector<Point3d>& points);
void loadFileXYZ(const char* filename, std::vector<Point3d>& points); //declare our own read-function, implementation is done below the main(...)-function

void drawPoints(std::vector<Point3d>& points, std::vector<Point3d>& pointColors, double pointSize);
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
std::vector<Point3d> resColors;
std::vector<Point3d*> ptrRes;
std::vector<Point3d> points;
std::vector<Point3d> pointsColors;
std::vector<Point3d> abfrage;
std::vector<Point3d> abfrageColors;
std::vector<Point3d> oldPoints;
std::vector<Point3d> oldPointsColors;
std::vector<Point3d> cornerPointsLine, cornerPointsPlane;
bool drawBestFitLine = false, drawBestFitPlane = false, doComputeBestFitLine = true,drawBestFitSphere=false;
double bestFitSphereRadius = 0;
Point3d bestFitSphereCenter;
KDTree data;
double abfrageLaenge;
int numNeighborhood = 2;
int startDim = 0;
int pointSize = 4;
//-----------------------------------------

int m_windowWidth = 0;
int m_windowHeight = 0;

double lastposX, lastposY; //last mouse positions

void initializeGL()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

Point3d colorFromGradientHSV(double index)
{
	if(index < 0) index = 0;
	else if(index > 1) index = 1;

	const double H(240.0*(1.0 - index));
	const double hi(std::floor(H / 60.0));
	const double f(H / 60 - hi);
	const double V(1.0);
	const double S(1.0);
	const double p(V*(1.0 - S));
	const double q(V*(1.0 - S*f));
	const double t(V*(1.0 - S*(1 - f)));

	if(hi == 1) return Point3d((255 * q), (255 * V), (255 * p));
	else if(hi == 2) return Point3d((255 * p), (255 * V), (255 * t));
	else if(hi == 3) return Point3d((255 * p), (255 * q), (255 * V));
	else if(hi == 4) return Point3d((255 * t), (255 * p), (255 * V));
	else if(hi == 5) return Point3d((255 * V), (255 * p), (255 * q));
	else return Point3d((255 * V), (255 * t), (255 * p));
}

void setDefaultPointColors(std::vector<Point3d>& colors, int size, Point3d color)
{
	colors.clear();

	for (int i = 0; i < size; i++)
	{
		colors.emplace_back(color);
	}
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

void action_loadFile(GLFWwindow* window) {
		//try to load point cloud data from file
		clock_t begin = clock();

		//loadFileXYZ("data/Stanford Dragon.xyz", points);
		points.clear();
		pointsColors.clear();
		cornerPointsLine.clear();
		cornerPointsPlane.clear();
		drawBestFitLine = false;
		drawBestFitPlane = false;
        drawBestFitSphere = false;
		loadFileXYZ(filename.c_str(), points); // FILENAME MOVED TO LINE 32

		//Tipp -> #pragma omp parallel for

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

		/* Make the window's context current */
		glfwMakeContextCurrent(window);

		//Prepare our virtual camera
		glfwGetFramebufferSize(window, &m_windowWidth, &m_windowHeight);

		//Initialize Camera
		m_camera.setWindowSize(m_windowWidth, m_windowHeight);    //setup window parameters
		m_camera.initializeCamera(m_sceneCenter, m_sceneRadius);  //set the camera outside the scene
		m_camera.updateProjection();                              //adjust projection to window size

																  //Farben der Punkte festlegen
		setDefaultPointColors(pointsColors, points.size(), Point3d(0, 255, 0));
}

void action_bestFitSphere(){
	/*
	*
	* Insert Spherefitting
	*
	*/
	// points.clear();
	// pointsColors.clear();

	// sum x
	//        std::vector<Point3d>
	//        points
	if (points.empty())
	{
		std::cout << "ERROR: cant execute bestFitSphere because points empty." <<  std::endl; 
	} else {
		/*
		 *
		 * Insert Spherefitting
		 *
		 */
		// points.clear();
		// pointsColors.clear();
		// sum x
		//        std::vector<Point3d>
		//        points
		clock_t begin = clock();
		double sum_x = 0, sum_y = 0, sum_z = 0, sum_xx = 0, sum_xy = 0, sum_xz = 0,
				sum_yy = 0, sum_yz = 0, sum_zz = 0, sum_xxx = 0, sum_xxy = 0,
				sum_xxz = 0, sum_xyy = 0, sum_xyz = 0, sum_xzz = 0, sum_yyy = 0,
				sum_yyz = 0, sum_yzz = 0, sum_zzz = 0;
		#pragma omp parallel for reduction(+: sum_x, sum_y, sum_z,sum_xx, sum_xy, sum_xz,sum_yy, sum_yz,sum_zz,sum_xxx, sum_xxy, sum_xxz, sum_xyy, sum_xyz, sum_xzz,sum_yyy, sum_yyz, sum_yzz, sum_zzz)
		for (int i = 0; i < points.size(); ++i) {
			double X = points[i].x;
			double XX = X * X;
			double Y = points[i].y;
			double YY = Y * Y;
			double Z = points[i].z;
			double ZZ = Z * Z;
			sum_x += X;
			sum_xx += XX;
			sum_xxx += XX * X;
			sum_y += Y;
			sum_yy += YY;
			sum_yyy += YY * Y;
			sum_z += Z;
			sum_zz += ZZ;
			sum_zzz += ZZ * Z;
			sum_xy += X * Y;
			sum_yz += Y * Z;
			sum_xz += X * Z;
			sum_xxy += XX * Y;
			sum_xxz += XX * Z;
			sum_xyy += YY * X;
			sum_yyz += YY * Z;
			sum_xzz += ZZ * X;
			sum_yzz += ZZ * Y;
		}
		double A1 = sum_xx + sum_yy + sum_zz;
		int n = points.size();
		double a = 2 * (sum_x * sum_x - n * sum_xx);
		double b = 2 * (sum_x * sum_y - n * sum_xy);
		double c = 2 * (sum_x * sum_z - n * sum_xz);
		double d = -n * (sum_xxx + sum_xyy + sum_xzz) + A1 * sum_x;
		double e = 2 * (sum_x * sum_y - n * sum_xy);
		double f = 2 * (sum_y * sum_y - n * sum_yy);
		double g = 2 * (sum_y * sum_z - n * sum_yz);
		double h = -n * (sum_xxy + sum_yyy + sum_yzz) + A1 * sum_y;
		double j = 2 * (sum_x * sum_z - n * sum_xz);
		double k = 2 * (sum_y * sum_z - n * sum_yz);
		double l = 2 * (sum_z * sum_z - n * sum_zz);
		double m = -n * (sum_xxz + sum_yyz + sum_zzz) + A1 * sum_z;
		double delta = a * (f * l - g * k) - e * (b * l - c * k)
				+ j * (b * g - c * f);
		double center_x = (d * (f * l - g * k) - h * (b * l - c * k)
				+ m * (b * g - c * f)) / delta;
		double center_y = (a * (h * l - m * g) - e * (d * l - m * c)
				+ j * (d * g - h * c)) / delta;
		double center_z = (a * (f * m - h * k) - e * (b * m - d * k)
				+ j * (b * h - d * f)) / delta;
		bestFitSphereCenter = Point3d(center_x, center_y, center_z);
		bestFitSphereRadius = sqrt(
				pow(center_x, 2) + pow(center_y, 2) + pow(center_z, 2)
						+ (A1
								- 2
										* (center_x * sum_x + center_y * sum_y
												+ center_z * sum_z)) / n);
		clock_t end = clock();
		std::cout << "Time needed to calculate Best Fit Sphere: "
				<< double(end - begin) / CLOCKS_PER_SEC << "s" << std::endl;
		std::vector<double> diffs;
		double maxValue = DBL_MIN, minValue = DBL_MAX;
		for (int i = 0; i < points.size(); i++) {
			diffs.emplace_back(
					abs(
							distance3d(bestFitSphereCenter, points[i])
									- bestFitSphereRadius));
			if (diffs[i] > maxValue || i == 0)
				maxValue = diffs[i];

			if (diffs[i] < minValue || i == 0)
				minValue = diffs[i];
		}
		maxValue = (double) (1) / (maxValue - minValue);
		//Normalisieren der Werte
		for (int i = 0; i < diffs.size(); i++) {
			diffs[i] = (diffs[i] - minValue) * maxValue;
		}
		//Stadnardabweichung/Mittelwert/Varianz
		double meanDiffs = 0, varianceDiffs = 0, standardDeviationDiffs = 0;
		for (int i = 0; i < diffs.size(); i++) {
			meanDiffs += diffs[i];
		}
		meanDiffs /= diffs.size();
		for (int i = 0; i < diffs.size(); i++) {
			varianceDiffs += pow(diffs[i] - meanDiffs, 2);
		}
		varianceDiffs /= diffs.size();
		standardDeviationDiffs = sqrt(varianceDiffs);
		std::cout << "Mittelwert: " << meanDiffs << std::endl;
		std::cout << "Varianz: " << varianceDiffs << std::endl;
		std::cout << "Standardabweichung: " << standardDeviationDiffs << std::endl;
		std::cout << "Radius: " << bestFitSphereRadius << std::endl;
		std::cout << "Center: " << bestFitSphereCenter.x << ", "
				<< bestFitSphereCenter.y << ", " << bestFitSphereCenter.z
				<< std::endl;
		for (int i = 0; i < points.size(); i++) {
			pointsColors[i] = colorFromGradientHSV(diffs[i]) * (1.0 / 255);
		}
		drawBestFitSphere = true;
	}
}

void action_rangeRequest() {
		//KDTree - RangeAbfrage
	if (points.empty())
	{
		std::cout << "ERROR: cant execute RangeReuest because points empty." <<  std::endl; 
	} else {
		//----------------------------------------------------------------------------
		abfrage.clear();
		int random = (std::rand() % (points.size()));
		abfrage.emplace_back(points[random]);

		Point3d S = m_bbmax - m_bbmin;
		abfrageLaenge = S.x * 0.25;

		clock_t begin = clock();
		ptrRes = data.getRange(abfrageLaenge, abfrage[0], startDim);

		res.clear();
		for(Point3d* point : ptrRes)
		{
			res.emplace_back(*point);
		}

		clock_t end = clock();

		setDefaultPointColors(resColors, res.size(), Point3d(255, 0, 0));
		setDefaultPointColors(abfrageColors, abfrage.size(), Point3d(125, 125, 0));

		std::cout << "Time needed to calculate range query: " << double(end - begin) / CLOCKS_PER_SEC << "s\r";
		//----------------------------------------------------------------------------
 }
}

void action_nearestNeighbor() {
	if (points.empty())
	{
		std::cout << "ERROR: cant execute NearestNeighbor because points empty." <<  std::endl; 
	} else {
		abfrage.clear();
		int random = (std::rand() % (points.size()));
		abfrage.emplace_back(points[random]);

		res.clear();

		clock_t begin = clock();
		ptrRes = data.getKNN(abfrage[0], numNeighborhood);

		res.clear();
		for(Point3d* point : ptrRes)
		{
			res.emplace_back(*point);
		}
		clock_t end = clock();

		setDefaultPointColors(resColors, res.size(), Point3d(255, 0, 0));
		setDefaultPointColors(abfrageColors, abfrage.size(), Point3d(125, 125, 0));

		std::cout << "Time needed to calculate NN: " << double(end - begin) / CLOCKS_PER_SEC << "s\r";
	}
}

void action_incNeighborhood() {
	//! check if Nearest Neighbour was already searched
	if (points.empty()){
		std::cout << "ERROR: cant execute incNeighborhood because points empty." <<  std::endl; 
	} else {
		if(abfrage.size()<1){
			//! if not, select random point
			abfrage.clear();
			int random = (std::rand() % (points.size()));
			abfrage.emplace_back(points[random]);
		}

		numNeighborhood++;

		ptrRes = data.getKNN(abfrage[0], numNeighborhood);

		res.clear();
		for(Point3d* point : ptrRes)
		{
			res.emplace_back(*point);
		}

		setDefaultPointColors(resColors, res.size(), Point3d(255, 0, 0));
		setDefaultPointColors(abfrageColors, abfrage.size(), Point3d(125, 125, 0));
	}
}

void action_decNeighborhood() {
	//! check if Nearest Neighbour was already searched
	if (points.empty()){
		std::cout << "ERROR: cant execute decNeighborhood because points empty." <<  std::endl; 
	} else {
		if(abfrage.size()<1){
			//! if not, select random point
			abfrage.clear();
			int random = (std::rand() % (points.size()));
			abfrage.emplace_back(points[random]);
		}

		if (numNeighborhood > 1)
		{
			numNeighborhood--;
		}
		ptrRes = data.getKNN(abfrage[0], numNeighborhood);

		res.clear();
		for(Point3d* point : ptrRes)
		{
			res.emplace_back(*point);
		}

		setDefaultPointColors(resColors, res.size(), Point3d(255, 0, 0));
		setDefaultPointColors(abfrageColors, abfrage.size(), Point3d(125, 125, 0));
	}
}

void action_smooth() {
	if (points.empty()){
		std::cout << "ERROR: cant execute smooth because points empty." <<  std::endl; 
	} else {
		clock_t begin = clock();
		oldPoints = points;
		setDefaultPointColors(oldPointsColors, oldPoints.size(), Point3d(255, 255, 255));
		points = data.smooth(points, numNeighborhood);

		//data = KDTree(points, startDim);
		clock_t end = clock();

		std::cout << "Time needed to smooth: " << double(end - begin) / CLOCKS_PER_SEC << "s\r";

		//Einfärben der Differenz zum Vorg�ngermodel
		std::vector<double> diffs;
		double maxValue = 0, minValue = 0;

		for (int i = 0; i < oldPoints.size(); i++)
		{
			diffs.emplace_back(sqrt(pow(points[i].x - oldPoints[i].x, 2)
				+ pow(points[i].y - oldPoints[i].y, 2)
				+ pow(points[i].z - oldPoints[i].z, 2)));

			if (diffs[i] > maxValue || i == 0)
				maxValue = diffs[i];

			if (diffs[i] < minValue || i == 0)
				minValue = diffs[i];
		}

		maxValue = (double)1 / (maxValue - minValue);

		//Normalisieren der Werte
		for (int i = 0; i < diffs.size(); i++)
		{
			diffs[i] = (diffs[i] - minValue) * maxValue;
		}

		//Übernahme der Werte
		for (int i = 0; i < oldPoints.size(); i++)
		{
			pointsColors[i] = colorFromGradientHSV(diffs[i]) * (1.0 / 255);
		}
	}
}

void action_thinning(){
	if (points.empty()){
		std::cout << "ERROR: cant execute thinning because points empty." <<  std::endl; 
	} else {
		clock_t begin = clock();

		data.thinning(numNeighborhood);
		points = data.getNotThinnedPoints();
		clock_t end = clock();

		data = KDTree(points, startDim);

		std::cout << points.size() << " points left." << std::endl;
		std::cout << "Time needed to thinn: " << double(end - begin) / CLOCKS_PER_SEC << "s\r";
	}
}

void action_bestFitLine(){
	if (points.empty()){
		std::cout << "ERROR: cant execute bestFitLine because points empty." <<  std::endl; 
	} else {
		std::vector<double> diffs;
		
		computeBestFitLine(points, cornerPointsLine);
		drawBestFitLine = true;
		drawBestFitPlane = false;

		Point3d pointOnLine = cornerPointsLine[0];
		Point3d lineDirection(cornerPointsLine[0].x - cornerPointsLine[1].x,
			cornerPointsLine[0].y - cornerPointsLine[1].y,
			cornerPointsLine[0].z - cornerPointsLine[1].z);

		double maxValue = DBL_MIN, minValue = DBL_MAX;

		for (int i = 0; i < points.size(); i++)
		{
			diffs.emplace_back(distancePt2Line(points[i], pointOnLine, lineDirection));

			if (diffs[i] > maxValue || i == 0)
				maxValue = diffs[i];

			if (diffs[i] < minValue || i == 0)
				minValue = diffs[i];
		}

		maxValue = (double)1 / (maxValue - minValue);

		//Normalisieren der Werte
		for (int i = 0; i < diffs.size(); i++)
		{
			diffs[i] = (diffs[i] - minValue) * maxValue;
		}
		
		// Colour:
		double meanDiffs = 0, varianceDiffs = 0, standardDeviationDiffs = 0;

		for (int i = 0; i < diffs.size(); i++)
		{
			meanDiffs += diffs[i];
		}

		meanDiffs /= diffs.size();

		for (int i = 0; i < diffs.size(); i++)
		{
			varianceDiffs += pow(diffs[i] - meanDiffs, 2);
		}

		varianceDiffs /= diffs.size();

		standardDeviationDiffs = sqrt(varianceDiffs);

		std::cout << "Mittelwert: " << meanDiffs << std::endl;
		std::cout << "Varianz: " << varianceDiffs << std::endl;
		std::cout << "Standardabweichung: " << standardDeviationDiffs << std::endl;

		for (int i = 0; i < points.size(); i++)
		{
			pointsColors[i] = colorFromGradientHSV(diffs[i]) * (1.0 / 255);
		}
	}
}

void action_bestFitPlane(){
	if (points.empty()){
		std::cout << "ERROR: cant execute bestFitPlane because points empty." <<  std::endl; 
	} else {
		std::vector<double> diffs;
		
		computeBestFitPlane(points, cornerPointsPlane);
		drawBestFitLine = false;
		drawBestFitPlane = true;

		Point3d pointOnPlane = cornerPointsPlane[0];
		Point3d direction1(cornerPointsPlane[0].x - cornerPointsPlane[1].x,
			cornerPointsPlane[0].y - cornerPointsPlane[1].y,
			cornerPointsPlane[0].z - cornerPointsPlane[1].z);
		Point3d direction2(cornerPointsPlane[0].x - cornerPointsPlane[2].x,
			cornerPointsPlane[0].y - cornerPointsPlane[2].y,
			cornerPointsPlane[0].z - cornerPointsPlane[2].z);
		Point3d planeDirection(direction1.y * direction2.z - direction2.y * direction1.z,
			direction1.z * direction2.x - direction2.z * direction1.x,
			direction1.x * direction2.y - direction2.x * direction1.y);

		double maxValue = DBL_MIN, minValue = DBL_MAX;

		for (int i = 0; i < points.size(); i++)
		{
			diffs.emplace_back(abs(distancePt2Plane(points[i], pointOnPlane, planeDirection)));

			if (diffs[i] > maxValue || i == 0)
				maxValue = diffs[i];

			if (diffs[i] < minValue || i == 0)
				minValue = diffs[i];
		}

		maxValue = (double)1 / (maxValue - minValue);

		//Normalisieren der Werte
		for (int i = 0; i < diffs.size(); i++)
		{
			diffs[i] = (diffs[i] - minValue) * maxValue;
		}
		
		//Colour:
		double meanDiffs = 0, varianceDiffs = 0, standardDeviationDiffs = 0;

		for (int i = 0; i < diffs.size(); i++)
		{
			meanDiffs += diffs[i];
		}

		meanDiffs /= diffs.size();

		for (int i = 0; i < diffs.size(); i++)
		{
			varianceDiffs += pow(diffs[i] - meanDiffs, 2);
		}

		varianceDiffs /= diffs.size();

		standardDeviationDiffs = sqrt(varianceDiffs);

		std::cout << "Mittelwert: " << meanDiffs << std::endl;
		std::cout << "Varianz: " << varianceDiffs << std::endl;
		std::cout << "Standardabweichung: " << standardDeviationDiffs << std::endl;

		for (int i = 0; i < points.size(); i++)
		{
			pointsColors[i] = colorFromGradientHSV(diffs[i]) * (1.0 / 255);
		}
	}
}


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	
	if (key == GLFW_KEY_1 && action == GLFW_RELEASE)
	{
		filename = "data/cone.xyz";
	}

	if (key == GLFW_KEY_2 && action == GLFW_RELEASE)
	{
		filename = "data/cap.xyz";
	}

	if (key == GLFW_KEY_3 && action == GLFW_RELEASE)
	{
		filename = "data/Stanford Horse.xyz";
	}

	if (key == GLFW_KEY_4 && action == GLFW_RELEASE)
	{
		filename = "data/Stanford Bunny.xyz";
	}

	if (key == GLFW_KEY_5 && action == GLFW_RELEASE)
	{
		filename = "data/Stanford Skeleton Hand.xyz";
	}

	if (key == GLFW_KEY_6 && action == GLFW_RELEASE)
	{
		filename = "data/Stanford Dragon.xyz";
	}

	if (key == GLFW_KEY_7 && action == GLFW_RELEASE)
	{
		filename = "data/Stanford Happy Buddha.xyz";
	}

	if (key == GLFW_KEY_8 && action == GLFW_RELEASE)
	{
		filename = "data/sphere.xyz";
	}

	if (key == GLFW_KEY_9 && action == GLFW_RELEASE)
	{
		filename = "data/sphere2.xyz";
	}

	if (key == GLFW_KEY_R && action == GLFW_RELEASE)
	{
		action_loadFile(window);
	}

	if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
	{
		action_rangeRequest();
	}

    if (key == GLFW_KEY_F && action == GLFW_RELEASE)
    {
		action_bestFitSphere();
	}

	if (key == GLFW_KEY_N && action == GLFW_RELEASE)
	{
		action_nearestNeighbor();
	}
	
	if (!points.empty())
	{
		if (key == GLFW_KEY_UP && action == GLFW_RELEASE)
		{
			action_incNeighborhood();
		}

		if (key == GLFW_KEY_DOWN && action == GLFW_RELEASE)
		{
			action_decNeighborhood();
		}

		if (key == GLFW_KEY_S && action == GLFW_RELEASE)
		{
			action_smooth();
		}

		if (key == GLFW_KEY_T && action == GLFW_RELEASE)
		{
			action_thinning();
		}

		if (key == GLFW_KEY_B && action == GLFW_RELEASE)
		{
			//Einfärben
			if (doComputeBestFitLine)
			{
				action_bestFitLine();
				doComputeBestFitLine = false;
			}
			else
			{
				action_bestFitPlane();
				doComputeBestFitLine = true;
			}
		}
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
	std::cout 	<< "usage 3ds_01_1 [ -nn <x> <y> <z> | -range <x> <y> <z> ] [-f <filename>]" << std::endl \
				<< std::endl \
				<< "	THIS DOESNT DO ANYTHING YET:" << std::endl \
				<< "	--nn <x> <y> <z>           -   find Nearest Neighbour for <x;y;z>" << std::endl \
				<< "	--range <r> <x> <y> <z>    -   highlight range <r> around <x;y;z>" << std::endl \
				<< "	-f <filename>              -   load specific xyz-File" << std::endl\
				<< "    --pointSizes <small> <big> -   set sizes for drawing points" << std::endl\
				<< std::endl;

	std::string job = "";
	double x = 0, y = 0, z = 0, r = 0;

	// parse arguments
	for(int i = 1; i<argc; i++){
		if(argv[i][0] == '-'){
			std::string option = argv[i];
			if( option == "-f"){
				i++;
				filename = argv[i];
			} else if( option == "--range") {
				job = "range";
				i++;
				r = std::stod(std::string(argv[i]));
				i++;
				x= std::stod(std::string(argv[i]));
				i++;
				y = std::stod(std::string(argv[i]));
				i++;
				z = std::stod(std::string(argv[i]));
			} else if( option == "--nn") {
				job = "nn";
				i++;
				x = std::stod(std::string(argv[i]));
				i++;
				y = std::stod(std::string(argv[i]));
				i++;
				z = std::stod(std::string(argv[i]));
			} else if ( option == "--pointSizes" ){
				i++;
				pointSize = std::stod(std::string(argv[i]));
				// i++;
				// pointSize_big = std::stod(std::string(argv[i]));
			}
		}
	}

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
		drawPoints(points, pointsColors, pointSize);
		//----------------------------------------------------------------------------

		//draw abfrage
		//----------------------------------------------------------------------------
		drawPoints(abfrage, abfrageColors, pointSize + 10);
		//----------------------------------------------------------------------------

		//draw res-points
		//----------------------------------------------------------------------------
		drawPoints(res, resColors, pointSize + 10);
		//----------------------------------------------------------------------------

		//draw old-points
		//----------------------------------------------------------------------------
		//drawPoints(oldPoints, oldPointsColors, pointSize);
		//----------------------------------------------------------------------------

		glPushAttrib(GL_LINE_BIT);
		glLineWidth(2);
		if (drawBestFitLine)
		{
			glColor3f(1.0f, 1.0f, 0.0f);
			glBegin(GL_LINES);
			glVertex3d(cornerPointsLine[0].x, cornerPointsLine[0].y, cornerPointsLine[0].z);
			glVertex3d(cornerPointsLine[1].x, cornerPointsLine[1].y, cornerPointsLine[1].z);
			glEnd();
		}
		if (drawBestFitPlane)
		{
			glColor3f(0.0f, 1.0f, 0.0f);
			glBegin(GL_LINE_LOOP);
			for (size_t i = 0; i < 4; ++i) {
				glVertex3d(cornerPointsPlane[i].x, cornerPointsPlane[i].y, cornerPointsPlane[i].z);
			}
			glEnd();
			glPointSize(8);
			glBegin(GL_POINTS);
			for (size_t i = 0; i < 4; ++i) {
				glVertex3d(cornerPointsPlane[i].x, cornerPointsPlane[i].y, cornerPointsPlane[i].z);
			}
			glEnd();
		}
        if (drawBestFitSphere)
        {
			//draw center point as a sphere
			glPushMatrix();
			glEnable (GL_BLEND); glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glColor4f(0, 0.4, 0, 0.7);
			glTranslated(bestFitSphereCenter.x, bestFitSphereCenter.y, bestFitSphereCenter.z);
			GLUquadric* quad = gluNewQuadric();
			gluSphere(quad, bestFitSphereRadius, 30, 30);
			gluDeleteQuadric(quad);
			glPopMatrix();
        }
		glPopAttrib();

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

void drawPoints(std::vector<Point3d>& points, std::vector<Point3d>& pointColors, double pointSize)
{
	glPointSize(pointSize);

	glEnable(GL_DEPTH_TEST);

	if (!points.empty())
	{ /* Drawing Points with VertexArrays */
		glEnableClientState(GL_COLOR_ARRAY);
		glEnableClientState(GL_VERTEX_ARRAY); //enable data upload to GPU
		glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &points[0]);

		if (!pointColors.empty())
			glColorPointer(3, GL_DOUBLE, sizeof(Point3d), &pointColors[0]);
		else
			glColor3ub(0, 0, 0);

		//draw point cloud
		glDrawArrays(GL_POINTS, 0, (unsigned int)points.size());
		glDisableClientState(GL_COLOR_ARRAY);  //disable data upload to GPU
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
#ifdef __linux__
    points.clear();
    std::ifstream file(filename);
    if(!file)
    {
        std::cout << "File not found\t " << filename << " is unavailable." << std::endl;
        return;
    }

    std::cout << "reading file: " << filename << std::endl;
    double a,b,c;
    while(file >> a >> b >> c)
//        points.push_back(Point3d(a,b,c));
        points.emplace_back(a,b,c);// equal to the cmd above

    size_t numberOfPoints = points.size();
    std::cout << "reading finished: " << numberOfPoints << " points have be read" << std::endl;

// compile for Windows if not Linux
#else
	FILE* file = 0;
	int error = fopen_s(&file, filename, "rt"); //r= read, t=text
	if (error != 0)
	{
		std::cout << "file " << filename << " could not be opened!" << std::endl;
		return; //nothing can be done else -> end function
	}

	std::cout << "reading file: " << filename << std::endl;

	points.clear();

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
#endif
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
