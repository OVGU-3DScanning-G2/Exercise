#include <random>
#include "../include/KDTree.h"

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

std::vector<Point3d> KDTree::getRange(double laenge, Point3d& point, int dim)
{
	std::vector<Point3d> res = std::vector<Point3d>();

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
		std::vector<Point3d> res2 = left->getRange(laenge, point, (dim + 1) % 3);
		res.insert(res.end(), res2.begin(), res2.end());
	}

	if (checkRight && right != NULL)
	{
		std::vector<Point3d> res2 = right->getRange(laenge, point, (dim + 1) % 3);
		res.insert(res.end(), res2.begin(), res2.end());
	}

	if ((left == NULL && right == NULL) || (checkLeft && checkRight))
	{
		bool push = false;

		if (median->x >= point.x - laenge && median->x <= point.x + laenge && median->y >= point.y - laenge && median->y <= point.y + laenge
			&& median->z >= point.z - laenge && median->z <= point.z + laenge)
			push = true;

		if(push)
			res.emplace_back(*(this->median));
	}

	return res;
}


double euclid(Point3d& p1, Point3d& p2)
{
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

bool samePoints(Point3d& p1, Point3d& p2)
{
	if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z)
		return true;
	else
		return false;
}

bool samePointInVector(Point3d& p, std::vector<Point3d>& points)
{
    for(Point3d point : points)
	{
		if (samePoints(p, point))
			return true;
	}

	return false;
}


std::vector<Point3d> KDTree::getKNN(Point3d& point, int k)
{
	if (k < 1)
	{
		k = 1;
		std::cout << "Automatically set k to 1, because it was smaller than 1." << std::endl;
	}

	std::vector<Point3d> points{point};

	for (int i = 0; i < k; i++)
	{
		points.emplace_back(getNN(points, 0));
	}

	points.erase(points.begin());

	return points;
}


Point3d KDTree::getNN(std::vector<Point3d>& points, int dim)
{
	//Rekursion nach unten
	//----------------------

	// go back if below leaf
	if (this == NULL)
		return Point3d(DBL_MAX, DBL_MAX, DBL_MAX); //Rückgabe von Unendlich beim Ankommen vom Ende des KDTree

	// go back if current leaf == point
	if (left == NULL && right == NULL)
	{
		if (samePointInVector(*median, points))
			return Point3d(DBL_MAX, DBL_MAX, DBL_MAX); //Rückgabe von Unendlich beim Ankommen vom
														//Blatt des KDTree (Punkt ist angefragter Punkt)
		else
			return *median; //Rückgabe des Punktes, wenn bei Blatt angekommen (nicht der angefragter Punkt)
	}

	bool goLeft = false;

	switch (dim) //Überprüfung, ob nach links gegangen werden muss oder nach rechts
	{
	case 0:
		if (points.front().x <= median->x) //X
			goLeft = true;
		break;
	case 1:
		if (points.front().y <= median->y) //Y
			goLeft = true;
		break;
	case 2:
		if (points.front().z <= median->z) //Z
			goLeft = true;
		break;
	}

	Point3d actual;

	if (goLeft)
	{
		actual = left->getNN(points, (dim + 1) % 3); //Ermittlung des nähesten Punktes im linken Teil
	}
	else
	{
		actual = right->getNN(points, (dim + 1) % 3); //Ermittlung des nähesten Punktes im rechten Teil
	}
	//----------------------

	//Recursion nach oben
	//----------------------
	double minDist = euclid(actual, points.front()); //Distanz zwischen ermitteltem Punkt und angefragten Punkt
	bool oldGoLeft = goLeft;

	switch(dim) //Überprüfung ob es Punkte im anderen Teilbaum gibt, die näher seien könnten
	{
	case 0:
		if (abs(median->x - points.front().x) <= minDist) //X
		{
			goLeft = !goLeft;
		}
		break;
	case 1:
		if (abs(median->y - points.front().y) <= minDist) //Y
		{
			goLeft = !goLeft;
		}
		break;
	case 2:
		if (abs(median->z - points.front().z) <= minDist) //Z
		{
			goLeft = !goLeft;
		}
		break;
	}

	if (oldGoLeft != goLeft) //Ermittlung des Punkte im anderen Teilbaums der dem angefragten Punkt am nächsten ist
	{
		Point3d otherBranchPoint;

		if (goLeft)
		{
			otherBranchPoint = left->getNN(points, (dim + 1) % 3); //Ermittlung des Punktes für den linken Teilbaum
		}
		else
		{
			otherBranchPoint = right->getNN(points, (dim + 1) % 3); //Ermittlung des Punktes für den rechten Teilbaum
		}

		if (euclid(otherBranchPoint, points.front()) < euclid(actual, points.front()) && !samePointInVector(otherBranchPoint, points))
			actual = otherBranchPoint; //Übernahme des Punktes, wenn er näher ist als der andere
	}

	if (euclid(*median, points.front()) < euclid(actual, points.front()) && !samePointInVector(*median, points))
		actual = *median; //Übernhame des Medians, sollte dieser noch näher dran sein

	return actual; //Übergabe des ermittelten Wertes
	//----------------------
}


std::vector<Point3d> KDTree::smooth(std::vector<Point3d>& points, int strength)
{
	std::vector<Point3d> newPoints;
	std::vector<Point3d> neighborHood;
	double newX = 0, newY = 0, newZ = 0;

	for(int i = 0; (unsigned)i < (unsigned)points.size(); i++)
	{
		neighborHood = getKNN(points[i], strength);

		double maxdist = euclid(points[i], neighborHood.back());
		double gewicht = 0;
		double sumGewichte = 0;

		for(int j = 0; (unsigned)j < (unsigned)neighborHood.size(); j++)
		{
			gewicht = exp(-1 * euclid(neighborHood[j], points[i]) / maxdist);
			newX += neighborHood[j].x * gewicht;
			newY += neighborHood[j].y * gewicht;
			newZ += neighborHood[j].z * gewicht;

			sumGewichte += gewicht;
		}

		newX += points[i].x;
		newY += points[i].y;
		newZ += points[i].z;

		sumGewichte++;

		Point3d newPoint = Point3d(newX / sumGewichte, newY / sumGewichte, newZ / sumGewichte);
		newPoints.emplace_back(newPoint);

		newX = 0;
		newY = 0;
		newZ = 0;
	}

	return newPoints;
}

std::vector<Point3d> KDTree::thinning(std::vector<Point3d>& points, int strength)
{
	std::vector<Point3d> newPoints;



	return newPoints;
}
