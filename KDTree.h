#include <vector>
#include "Point3d.h"

class KDTree
{
public:
	KDTree(std::vector<Point3d>& points, int dim = 0);
	std::vector<Point3d> abfrage(int laenge, Point3d& point);
	Point3d value;
	KDTree* left;
	KDTree* right;
protected:
	// dim: dimension to spit at
};