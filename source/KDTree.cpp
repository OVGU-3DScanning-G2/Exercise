#include "KDTree.h"
#include <iostream>
#include <algorithm>

bool sortByXvalue(const Point3d& p1, const Point3d& p2)
{
	if (p1.x < p2.x)
		return true;
	else
		return false;
}

bool sortByYvalue(const Point3d& p1, const Point3d& p2)
{
	if (p1.y < p2.y)
		return true;
	else
		return false;
}

bool sortByZvalue(const Point3d& p1, const Point3d& p2)
{
	if (p1.z < p2.z)
		return true;
	else
		return false;
}

KDTree::KDTree(std::vector<Point3d>& points, int dim){
	if (points.size() == 1)
	{
		value = points[0];
		left = NULL;
		right = NULL;
	}
	else {
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

		if (points.size() == 2)
		{
			value = points[0];
			left = new KDTree(*new std::vector<Point3d>(points.begin(), points.begin() + 1), (dim + 1) % 3);
			right = new KDTree(*new std::vector<Point3d>(points.begin() + 1, points.end()), (dim + 1) % 3);
		}
		else
		{
			std::size_t half_size = points.size() / 2;

			if (points.size() % 2 == 0)
				half_size = half_size - 1;

			std::vector<Point3d>* split_lo = new std::vector<Point3d>(points.begin(), points.begin() + half_size + 1);
			std::vector<Point3d>* split_hi = new std::vector<Point3d>(points.begin() + half_size + 1, points.end());

			value = points[half_size];
			left = new KDTree(*split_lo, (dim + 1) % 3);
			right = new KDTree(*split_hi, (dim + 1) % 3);
		}

		//std::cout << this->value.x << " " << this->value.y << " " << this->value.z << std::endl;
	}
}

std::vector<Point3d> KDTree::abfrage(double laenge, Point3d& point, int dim)
{
	std::vector<Point3d> res = std::vector<Point3d>();

	bool checkLeft = false;
	bool checkRight = false;

	switch (dim)
	{
	case 0:
		if (value.x >= point.x - laenge)
			checkLeft = true;
		if (value.x <= point.x + laenge)
			checkRight = true;
		break;
	case 1:
		if (value.y >= point.y - laenge)
			checkLeft = true;
		if (value.y <= point.y + laenge)
			checkRight = true;
		break;
	case 2:
		if (value.z >= point.z - laenge)
			checkLeft = true;
		if (value.z <= point.z + laenge)
			checkRight = true;
		break;
	}

	if (checkLeft && left != NULL)
	{
		std::vector<Point3d> res2 = left->abfrage(laenge, point, (dim + 1) % 3);
		res.insert(res.end(), res2.begin(), res2.end());
	}

	if (checkRight && right != NULL)
	{
		std::vector<Point3d> res2 = right->abfrage(laenge, point, (dim + 1) % 3);
		res.insert(res.end(), res2.begin(), res2.end());
	}

	if (left == NULL && right == NULL)
	{
		bool push = false;

		if (value.x >= point.x - laenge && value.x <= point.x + laenge && value.y >= point.y - laenge && value.y <= point.y + laenge
			&& value.z >= point.z - laenge && value.z <= point.z + laenge)
			push = true;

		if(push)
			res.push_back(this->value);
	}
		
	return res;
}