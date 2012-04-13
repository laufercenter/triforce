/* Copyright 2011, Nils J. D. Drechsel & Jordi Villà-Freixa
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef TESSELLATION_H_
#define TESSELLATION_H_


#include <armadillo>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <map>
#include <vector>
#include <list>
#include <set>

#include "molecule.h"



using namespace std;
using namespace arma;



#define CONCAVE 0
#define CONVEX 1
#define ANGLE_TRESHOLD 0.001
#define THRESHOLD_IP 0.05
#define THRESHOLD_NUMERICAL 0.0001
#define EPSILON 0.01
#define ORDER_CLOCKWISE 0
#define ORDER_COUNTERCLOCKWISE 1




typedef vec Vector;




typedef struct
{
	double matrix[3][3];
}
Matrix3D;


typedef struct
{
	Vector vector;
	int with;
	int from;
	bool visited;
	bool flagged;
	int id;
}
IntersectionPoint;


typedef struct 
{
	int id;
	int index;
	Vector vector;
	Vector normal;
	double openingAngle;
	double g;
	double a;
	double sphereRadius;
//	double radius;
	std::list<IntersectionPoint> forwardIntersections;
	int form;
	bool intersect;
		
	
}
CircularRegion;

typedef struct 
{
	Vector k_j;
	Vector j_k;	
}
IntersectionPair;





class Tesselation{
	
public:
	Tesselation(Molecule *m);
	void build();
	
	
private:
	
	Molecule *molecule;
	vector<Vector> atoms;
	


	double cot(double a);
	double csc(double a);
	double getAngleBetweenNormals(Vector &a, Vector &b);
	double getAngle(Vector &a, Vector &b);
	bool isZero(double v);
	void determineProjection(Vector &origin, double radius, CircularRegion &circle);
	IntersectionPair determineIntersectionPoints(double radius, CircularRegion &K, CircularRegion &J);
	void makeCircularRegions(Vector &origin, double radius, std::vector<vec> &atoms, std::vector<int> &indizee, std::vector<double> &radii, std::vector<CircularRegion> &circles);
	void filterCircularRegions(double radius, std::vector<CircularRegion> &circles);
	void filterEmptyCircularRegions(std::vector<CircularRegion> &circles);
	void reindexCircularRegions(std::vector<CircularRegion> &circles);
	void outputCircularRegions(std::vector<CircularRegion> &circles);
	void filterIntersectionPoints(std::vector<CircularRegion> &circles, int except);
	void clearFlags(std::vector<CircularRegion> &circles);
	std::vector<CircularRegion>* deepCopy(std::vector<CircularRegion> &circles);
	void outputGaussBonnetPath(std::list<IntersectionPoint*> &points);
	void prepareCircularRegions(std::vector<CircularRegion> &circles, std::vector<CircularRegion> **newCircles);
	void insertFakeIntersectionPoints(std::vector<CircularRegion> &circles);
	void buildGaussBonnetPath(Vector &origin, double radius, std::vector<vec> &atoms, std::vector<int> &indizee, std::vector<double> &radii, std::vector<CircularRegion> &circles);
	void harvestIntersectionPoints(std::vector<CircularRegion> &circles, std::vector<vec> &intersections);
	bool hasUnflaggedIntersectionPoints(CircularRegion &circle, IntersectionPoint **ip);
	std::list<IntersectionPoint*>* retrieveIntersections(CircularRegion &circle);
	void showIntersections(std::list<IntersectionPoint*> &intersections);
	std::vector<std::list<IntersectionPoint*>*>*  harvestGaussBonnetPaths(std::vector<CircularRegion> &circles);


	
	
	
};



#endif // TESSELLATION_H_
