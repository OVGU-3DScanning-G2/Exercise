#include "include/KDTree.h"
#include <iostream>
#include <algorithm>

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
	left = NULL;
	right = NULL;
}

KDTree::KDTree(std::vector<Point3d>& points, int dim){
	if (points.size() == 1)
	{
		median = points[0];
		left = NULL;
		right = NULL;
	}
	else {
		// sort points
		switch (dim) {
		case 0:
			std::sort(points.begin(), points.end(), sortByXvalue);
			break;
		case 1:
			std::sort(points.begin(), points.end(), sortByYvalue);
			break;
		case 2:
			std::sort(points.begin(), points.end(), sortByZvalue);
			break;
		}

		// find median point
		//  should be point in center of array after sorting
		if (points.size() == 2)
		{
			median = points[0];
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

			median = points[half_size];
			left = new KDTree(*split_lo, (dim + 1) % 3);
			right = new KDTree(*split_hi, (dim + 1) % 3);
		}

		//std::cout << this->value.x << " " << this->value.y << " " << this->value.z << std::endl;
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
		if (median.x >= point.x - laenge)
			checkLeft = true;
		if (median.x <= point.x + laenge)
			checkRight = true;
		break;
	case 1:
		if (median.y >= point.y - laenge)
			checkLeft = true;
		if (median.y <= point.y + laenge)
			checkRight = true;
		break;
	case 2:
		if (median.z >= point.z - laenge)
			checkLeft = true;
		if (median.z <= point.z + laenge)
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

		if (median.x >= point.x - laenge && median.x <= point.x + laenge && median.y >= point.y - laenge && median.y <= point.y + laenge
			&& median.z >= point.z - laenge && median.z <= point.z + laenge)
			push = true;

		if(push)
			res.emplace_back(this->median);
	}

	return res;
}

Point3d KDTree::getNN(Point3d& point)
{
	double startMinDist = pow(median.x - point.x, 2) + pow(median.y - point.y, 2) + pow(median.z - point.z, 2);

	return getNN(point, startMinDist);
}

Point3d KDTree::getNN(Point3d& point, double minDist)
{
	double currLeftDist = -1;
	if(left != NULL && point.x != left->median.x && point.y != left->median.y && point.z != left->median.z)
		currLeftDist = pow(left->median.x - point.x, 2) + pow(left->median.y - point.y, 2) + pow(left->median.z - point.z, 2);
	double currRightDist = -1; 
	if(right != NULL && point.x != right->median.x && point.y != right->median.y && point.z != right->median.z)
		currRightDist = pow(right->median.x - point.x, 2) + pow(right->median.y - point.y, 2) + pow(right->median.z - point.z, 2);

	double currMinDist = 0;

	if (currLeftDist > currRightDist && currRightDist > 0)
		currMinDist = currRightDist;
	else
	{
		if (currLeftDist > 0)
			currMinDist = currLeftDist;
		else
			currMinDist = currRightDist;
	}

	if (minDist <= currMinDist || (currLeftDist < 0 && currRightDist < 0))
		return median;
	else
	{
		if (currMinDist == currLeftDist)
			return left->getNN(point, currMinDist);
		else
			return right->getNN(point, currMinDist);
	}
}
