/*!
   \file "KDTree.cpp"
   \brief "Implementation of KDTree object."
*/

#include "../include/KDTree.h"
#include <random>
static Point3d searchPoint;

bool sortByXvalue(const Point3d& p1, const Point3d& p2)
{
	return p1.x < p2.x;
}

bool sortByYvalue(const Point3d& p1, const Point3d& p2)
{
	return p1.y < p2.y;
}

bool sortByZvalue(const Point3d& p1, const Point3d& p2)
{
	return p1.z < p2.z;
}

KDTree::KDTree()
{
	median = NULL;
	left = NULL;
	right = NULL;
}

KDTree::KDTree(std::vector<Point3d>& points, int dim){
	if (points.size() == 1)
	{
		median = &points[0];
		left = NULL;
		right = NULL;
	}
	else {
		// sort points
		int median_index = points.size() / 2;
		switch (dim) {
		case 0:
			std::nth_element(points.begin(), points.begin() + median_index,
				points.end(), sortByXvalue);
			break;
		case 1:
			std::nth_element(points.begin(), points.begin() + median_index,
				points.end(), sortByYvalue);
			break;
		case 2:
			std::nth_element(points.begin(), points.begin() + median_index,
				points.end(), sortByZvalue);
			break;
		}

		// find median point
		//  should be point in center of array after sorting
		if (points.size() == 2)
		{
			median = &points[0];
			left = NULL;
			right = new KDTree(*new std::vector<Point3d>(points.begin() + 1, points.end()), (dim + 1) % 3);
		}
		else
		{
			std::size_t half_size = points.size() / 2;

			if (points.size() % 2 == 0)
				half_size = half_size - 1;

			std::vector<Point3d>* split_lo = new std::vector<Point3d>(points.begin(), points.begin() + half_size);
			std::vector<Point3d>* split_hi = new std::vector<Point3d>(points.begin() + half_size + 1, points.end());

			median = &points[half_size];
			left = new KDTree(*split_lo, (dim + 1) % 3);
			right = new KDTree(*split_hi, (dim + 1) % 3);
		}
	}
}

std::vector<Point3d*> KDTree::getRange(double laenge, Point3d& point, int dim)
{
	std::vector<Point3d*> res = std::vector<Point3d*>();

	bool checkLeft = false;
	bool checkRight = false;

	switch (dim)
	{
	case 0:
		if (median->x >= point.x - laenge)
			checkLeft = true;
		if (median->x <= point.x + laenge)
			checkRight = true;
		break;
	case 1:
		if (median->y >= point.y - laenge)
			checkLeft = true;
		if (median->y <= point.y + laenge)
			checkRight = true;
		break;
	case 2:
		if (median->z >= point.z - laenge)
			checkLeft = true;
		if (median->z <= point.z + laenge)
			checkRight = true;
		break;
	}

	if (checkLeft && left != NULL)
	{
		std::vector<Point3d*> res2 = left->getRange(laenge, point, (dim + 1) % 3);
		res.insert(res.end(), res2.begin(), res2.end());
	}

	if (checkRight && right != NULL)
	{
		std::vector<Point3d*> res2 = right->getRange(laenge, point, (dim + 1) % 3);
		res.insert(res.end(), res2.begin(), res2.end());
	}

	if ((left == NULL && right == NULL) || (checkLeft && checkRight))
	{
		bool push = false;

		if (median->x >= point.x - laenge && median->x <= point.x + laenge && median->y >= point.y - laenge && median->y <= point.y + laenge
			&& median->z >= point.z - laenge && median->z <= point.z + laenge)
			push = true;

		if(push)
			res.emplace_back(this->median);
	}

	return res;
}


static double euclid(Point3d& p1, Point3d& p2)
{
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

bool samePointInVector(Point3d& p, std::vector<Point3d>& points)
{
	for(Point3d point : points)
	{
		if (p == point)
			return true;
	}

	return false;
}


std::vector<Point3d*> KDTree::getKNN(Point3d& point, int k)
{
	if (k < 1)
	{
		k = 1;
		std::cout << "Automatically set k to 1, because it was smaller than 1." << std::endl;
	}

	searchPoint = point;

	std::vector<Point3d*> neighbours{new Point3d(DBL_MAX, DBL_MAX, DBL_MAX)};

	getNN(point, neighbours, k, 0);

	return neighbours;
}

bool sortWithSearchPoint(Point3d* p1, Point3d* p2)
{
	return (euclid(*p1, searchPoint) < euclid(*p2, searchPoint));
}

void insertNeighbour(Point3d& insertPoint, std::vector<Point3d*>& neighbours, int k)
{
	neighbours.push_back(&insertPoint);

	//Sortieren
	std::sort(neighbours.begin(), neighbours.end(), sortWithSearchPoint);

	if (neighbours.size() > k)
		neighbours.pop_back();
}

void KDTree::getNN(Point3d& point, std::vector<Point3d*>& neighbours, int k, int dim)
{
	//Rekursion nach unten
	//----------------------
	if (this == NULL)
		return; //R�ckgabe von Unendlich beim Ankommen vom Ende des KDTree

	if (left == NULL && right == NULL)
	{
		if (point == *median)
		{
			return;													  //Rückgabe von Unendlich beim Ankommen vom
																	  //Blatt des KDTree (Punkt ist angefragter Punkt)
		}
		else
		{
			insertNeighbour(*median, neighbours, k);
			return; //Rückgabe des Punktes, wenn bei Blatt angekommen (nicht der angefragter Punkt)
		}
	}

	bool goLeft = false;

	switch (dim) //�berpr�fung, ob nach links gegangen werden muss oder nach rechts
	{
	case 0:
		if (point.x <= median->x) //X
			goLeft = true;
		break;
	case 1:
		if (point.y <= median->y) //Y
			goLeft = true;
		break;
	case 2:
		if (point.z <= median->z) //Z
			goLeft = true;
		break;
	}

	//Point3d actual;

	if (goLeft)
	{
		left->getNN(point, neighbours, k, (dim + 1) % 3); //Ermittlung des n�hesten Punktes im linken Teil
	}
	else
	{
		right->getNN(point, neighbours, k, (dim + 1) % 3); //Ermittlung des n�hesten Punktes im rechten Teil
	}
	//----------------------

	//Recursion nach oben
	//----------------------
	double maxMinDist = euclid(*(neighbours.back()), point); //Distanz zwischen ermitteltem Punkt und angefragten Punkt
	bool oldGoLeft = goLeft;

	switch(dim) //�berpr�fung ob es Punkte im anderen Teilbaum gibt, die n�her seien k�nnten
	{
	case 0:
		if (abs(median->x - point.x) <= maxMinDist) //X
		{
			goLeft = !goLeft;
		}
		break;
	case 1:
		if (abs(median->y - point.y) <= maxMinDist) //Y
		{
			goLeft = !goLeft;
		}
		break;
	case 2:
		if (abs(median->z - point.z) <= maxMinDist) //Z
		{
			goLeft = !goLeft;
		}
		break;
	}

	if (oldGoLeft != goLeft) //Ermittlung des Punkte im anderen Teilbaums der dem angefragten Punkt am n�chsten ist
	{
		//Point3d otherBranchPoint;

		if (goLeft)
		{
			left->getNN(point, neighbours, k, (dim + 1) % 3); //Ermittlung des Punktes f�r den linken Teilbaum
		}
		else
		{
			right->getNN(point, neighbours, k, (dim + 1) % 3); //Ermittlung des Punktes f�r den rechten Teilbaum
		}

		//if (euclid(otherBranchPoint, point) < euclid(actual, point) && !samePointInVector(otherBranchPoint, neighbours))
		//	actual = otherBranchPoint; //�bernahme des Punktes, wenn er n�her ist als der andere
	}

	if (euclid(*median, point) < maxMinDist && !(*median == point))
		insertNeighbour(*median, neighbours, k); //�bernhame des Medians, sollte dieser noch n�her dran sein

	//inserNeighbour(actual, neighbours, k);
	return; //�bergabe des ermittelten Wertes
	//----------------------
}


std::vector<Point3d> KDTree::smooth(std::vector<Point3d>& points, int strength)
{
	std::vector<Point3d> newPoints;
	std::vector<Point3d*> neighborHood;
	double newX = 0, newY = 0, newZ = 0;

	for(int i = 0; (unsigned)i < (unsigned)points.size(); i++)
	{
		neighborHood = getKNN(points[i], strength);

		double maxdist = euclid(points[i], *neighborHood.back());
		double gewicht = 0;
		double sumGewichte = 0;

		for(int j = 0; (unsigned)j < (unsigned)neighborHood.size(); j++)
		{
			gewicht = exp(-1 * euclid(*neighborHood[j], points[i]) / maxdist);
			newX += neighborHood[j]->x * gewicht;
			newY += neighborHood[j]->y * gewicht;
			newZ += neighborHood[j]->z * gewicht;

			sumGewichte += gewicht;
		}

		//newX += points[i].x;
		//newY += points[i].y;
		//newZ += points[i].z;

		//sumGewichte++;

		Point3d newPoint = Point3d(newX / sumGewichte, newY / sumGewichte, newZ / sumGewichte);
		newPoints.emplace_back(newPoint);

		newX = 0;
		newY = 0;
		newZ = 0;
	}

	return newPoints;
}

void KDTree::thinning(int strength)
{
	//Diese Methode ist fehlerhaft (Springt f�r jeden Punkt rein, bis auf einen???, obwohl diese im KDTree richtig markiert werden)
	/*for (int i = 0; i < points.size(); i++)
	{
		if (points[i].thinned == false)
		{
			std::vector<Point3d*> neighbours = getKNN(points[i], strength);

			for (int j = 0; j < neighbours.size(); j++)
			{
				neighbours[j]->thinned = true;
			}
		}
		else
		{
			std::cout << "Point already thinned" << std::endl;
		}
	}*/

	if (this == NULL)
		return;

	if (median->thinned == false)
	{
		std::vector<Point3d*> neighbours = getKNN(*median, strength);

		for(Point3d* point : neighbours)
		{
			point->thinned = true;
		}
	}

	left->thinning(strength);
	right->thinning(strength);
}

std::vector<Point3d> KDTree::getNotThinnedPoints()
{
	std::vector<Point3d> res;

	if (this == NULL)
		return res;

	if (median->thinned == false)
		res.emplace_back(*median);

	std::vector<Point3d> res2 = right->getNotThinnedPoints();
	res.insert(res.end(), res2.begin(), res2.end());

	res2 = left->getNotThinnedPoints();
	res.insert(res.end(), res2.begin(), res2.end());

	return res;
}
