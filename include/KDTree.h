#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <iostream>

#include "include/Point3d.h"

class KDTree
{
public:
	KDTree();
	KDTree(std::vector<Point3d>& points, int dim = 0);
	std::vector<Point3d> getRange(double laenge, Point3d& point, int dim);
	std::vector<Point3d> getKNN(Point3d& point, int k);
	Point3d median;
	KDTree* left;
	KDTree* right;
protected:
	// dim: dimension to spit at
private:
	Point3d getNN(std::vector<Point3d>& point, int dim);
};
