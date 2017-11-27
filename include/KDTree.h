#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <iostream>

#include "include/Point3d.h"

class KDTree
{
public:
	KDTree(); /*!< Erzeugt einen leeren KDTree mit NullPtr als median, left und right. */
	KDTree(std::vector<Point3d>& points, int dim = 0); /*!< Erzeugt einen KDTree fuer die Punktmenge points. */
	std::vector<Point3d*> getRange(double laenge, Point3d& point, int dim); /*!< Bestimmt die Punkte die sich innerhalb eines Quaders mit den Massen
																			leange x laenge x laenge und dem Zentrum point befinden.*/
	std::vector<Point3d*> getKNN(Point3d& point, int k); /*!< Berechnet die k naechsten Nachbarn des Punktes point mit zuhilfenahme der Funktion getNN. */
	std::vector<Point3d> smooth(std::vector<Point3d>& points, int strength); /*!< Smoothed die Punkte des KDTree anhand des strength Kriteriums,
																			 wobei die k=strength Nachbarn fuer die Kalkulation herangezugen werden.*/
	void thinning(int strength); /*!< Duennt den KDTree mit der Staerke strength (es werden die k=strength Nachbarn entfernt) aus. */
	std::vector<Point3d> getNotThinnedPoints(); /*!< Gibt alle nicht markierten Punkte des KDTree zurueck. */
	Point3d* median; /*!< Speichert den Median der Menge fuer die beiden Kindknoten. */
	KDTree* left; /*!< Linker Kindknoten des KDTree. */
	KDTree* right; /*!< Rechter Kindknoten des KDTree. */
protected:
	// dim: dimension to spit at
private:
	void getNN(Point3d& point, std::vector<Point3d*>& neighbours, int k, int dim); /*!< Berechnet die k naechsten Nachbarn des Punktes point. */
};
