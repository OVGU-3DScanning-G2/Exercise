#ifndef MY_ALGORITHMS_H
#define MY_ALGORITHMS_H

#include "../include/Point3d.h"
#include "../include/Matrix.h"
#include <iostream>
#include <vector>

Point3d computeCenter(const std::vector<Point3d>& points);                  ///< Computes and returns the center of the point cloud
void computeCovarianceMatrix3x3(const std::vector<Point3d>& points, Matrix& M); ///< Coputes the 3x3 covariance matrix
void computeBestFitLine (const std::vector<Point3d>& points, std::vector<Point3d>& corners);  ///< Computes best-fit line
void computeBestFitPlane(const std::vector<Point3d>& points, std::vector<Point3d>& corners);  ///< Computes best-fit plane

double distancePt2Line (const Point3d& point, const Point3d& pointOnLine , const Point3d& lineDirection );  ///< distance point-to-line (3d)
double distancePt2Plane(const Point3d& point, const Point3d& pointOnPlane, const Point3d& planeDirection);  ///< distance point-to-plane

#endif //MY_ALGORITHMS_H
