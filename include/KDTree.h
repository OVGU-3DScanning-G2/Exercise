#include <vector>
#include "include/Point3d.h"

class KDTree
{
public:
	KDTree();
	KDTree(std::vector<Point3d>& points, int dim = 0);
	std::vector<Point3d> getRange(double laenge, Point3d& point, int dim);
	Point3d getNN(Point3d& point);
	Point3d median;
	KDTree* left;
	KDTree* right;
protected:
	// dim: dimension to spit at
private:
	Point3d getNN(Point3d& point, double minDist);
};
