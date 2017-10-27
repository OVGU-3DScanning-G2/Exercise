#include "KDTree.h"


KDTree::KDTree(std::vector<Point3d>& points, int dim = 0){
	if (points.size() == 1)
	{
		this.value = points[0];
		this.left = NULL;
		this.right = NULL;
	}
	else{
		switch(dim){
			case 0:
				std::sort( points.begin(),points.end(),sortByXvalue);
				break;
			case 1:
				std::sort( points.begin(),points.end(),sortByYvalue);
				break;
			case 2:
				std::sort( points.begin(),points.end(),sortByZvalue);
				break;
		}
		std::size_t const half_size = lines.size() / 2;
		std::vector<Point3d> split_lo(points.begin(), points.begin() + half_size);
		std::vector<Point3d> split_hi(points.begin() + half_size+1, points.end());
		this.left = KDTree(split_lo, (dim + 1 )%3);
		this.right = KDTree(split_hi, (dim + 1 )%3);
		this.value = points[half_size];
		std::cout << "Value left: " <<  this.left.value << "\nValue right: " <<  this.right.value << "\n--------" <<std::endl; 
	}
}

std::vector<Point3d> KDTree::abfrage(int laenge, Point3d& point)
{
	//Code für die Abfrage der Punkte die innerhalb des Rechtecks liegen hier rein
}
bool sortByXvalue ( const Point3d& p1, const Point3d& p2)
{
	if (p1.x < p2.x)
		return true;
	else
		return false;
}

bool sortByYvalue ( const Point3d& p1, const Point3d& p2)
{
	if (p1.y < p2.y)
		return true;
	else
		return false;
}

bool sortByZvalue ( const Point3d& p1, const Point3d& p2)
{
	if (p1.z < p2.z)
		return true;
	else
		return false;
}