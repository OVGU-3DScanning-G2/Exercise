#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <iostream>

#include "Point3d.h"

class KDTree
{
public:
	KDTree();
	KDTree(std::vector<Point3d>& points, int dim = 0);
	std::vector<Point3d*> getRange(double laenge, Point3d& point, int dim);
	std::vector<Point3d*> getKNN(Point3d& point, int k);
	std::vector<Point3d> smooth(std::vector<Point3d>& points, int strength);
	void thinning(int strength);
	std::vector<Point3d> getNotThinnedPoints();
	Point3d* median;
	KDTree* left;
	KDTree* right;
protected:
	// dim: dimension to spit at
private:
	void getNN(Point3d& point, std::vector<Point3d*>& neighbours, int k, int dim);
};
