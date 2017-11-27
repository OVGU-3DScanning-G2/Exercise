/*!
   \file "KDTree.cpp"
   \brief "Implementation of KDTree object."
*/

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <iostream>

#include "Point3d.h"


class KDTree
{
public:
	KDTree(); /*!< Erzeugt einen leeren KDTree mit NullPtr als median, left und right.*/
	KDTree(std::vector<Point3d>& points, int dim = 0); /*!< Erzeugt einen KDTree für die Punktmenge points. */
	std::vector<Point3d*> getRange(double laenge, Point3d& point, int dim); /*!< Bestimmt die Punkte die sich innerhalb eines Quaders mit den Maßen
																			länge x länge x länge und dem Zentrum point befinden.*/
	std::vector<Point3d*> getKNN(Point3d& point, int k); /*!< Berechnet die k nächsten Nachbarn des Punktes point mit zuhilfenahme der Funktion getNN. */
	std::vector<Point3d> smooth(std::vector<Point3d>& points, int strength); /*!< Smoothed die Punkte des KDTree anhand des strength Kriteriums,
																			 wobei die k=strength Nachbarn für die Kalkulation herangezugen werden.*/
	void thinning(int strength); /*!< Dünnt den KDTree mit der Stärke strength (es werden die k=strength Nachbarn entfernt) aus. */
	std::vector<Point3d> getNotThinnedPoints(); /*!< Gibt alle nicht markierten Punkte des KDTree zurück. */
	Point3d* median; /*!< Speichert den Median der Menge für die beiden Kindknoten. */
	KDTree* left; /*!< Linker Kindknoten des KDTree. */
	KDTree* right; /*!< Rechter Kindknoten des KDTree. */
protected:
	// dim: dimension to spit at
private:
	void getNN(Point3d& point, std::vector<Point3d*>& neighbours, int k, int dim); /*!< Berechnet die k nächsten Nachbarn des Punktes point. */
};
/*! \class KDTree KDTree.h "include/KDTree.h"
	\brief KDTree

	This is the KD-Tree.
 */
