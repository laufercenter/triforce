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
#include <stdio.h>

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

enum OcclusionState{
	OCCLUDED,
	UNOBSTRUCTED,
	UNDEFINED
};



typedef struct
{
	double d;
	bool blocked;
	map<int,OcclusionState> tertiaryIntersections;
	bool sasa;
	bool visited;
	
	
}
CircularIntersection;

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
	list<IntersectionPoint> forwardIntersections;
	map<int,CircularIntersection> circularIntersections;
	int form;
	bool intersect;
	bool flagged;
		
	
}
CircularRegion;

typedef struct 
{
	Vector k_j;
	Vector j_k;	
}
IntersectionPair;


enum OcclusionType{
	BOTTOM,
	TOP,
	BOTH,
	INSIDE,
	OUTSIDE,
	TRIFORCE,
	NA
};


class Tessellation{
	
public:
	Tessellation(Molecule &m);
	void build();
	vector<vector<list<IntersectionPoint*>*>*>* intersectionPoints();
	vector<vector<CircularRegion>* >* circularRegions();
	void setOutputFile(string filename);
	void outputGaussBonnetData(string filename);
	
	
private:
	
	Molecule molecule;
	vector<Vector> atoms;
	vector<double> *radii;
	//#atoms #circularregions
	vector<vector<CircularRegion>* > circleSet;
	//#atoms #sasas #circularregions
	vector<vector<list<IntersectionPoint*>*>*> intersectionSet;


	void emptyIntersections();
	void emptyCircularRegions();
	double vsign(double v);
	double cot(double a);
	double csc(double a);
	double getAngleBetweenNormals(Vector &a, Vector &b);
	double getAngle(Vector &a, Vector &b);
	bool isZero(double v);
	void determineProjection(Vector &origin, double radius, CircularRegion &circle);
	IntersectionPair determineIntersectionPoints(double radius, CircularRegion &K, CircularRegion &J);
	void makeCircularRegions(Vector &origin, double radius, vector<vec> &atoms, vector<double> &radii, vector<CircularRegion> &circles);
	int filterCircularRegions(double radius, vector<CircularRegion> &circles);
	void filterEmptyCircularRegions(vector<CircularRegion> &circles);
	void reindexCircularRegions(vector<CircularRegion> &circles);
	void outputCircularRegions(vector<CircularRegion> &circles);
	void filterIntersectionPoints(vector<CircularRegion> &circles, int except);
	void clearFlags(vector<CircularRegion> &circles);
	vector<CircularRegion>* deepCopy(vector<CircularRegion> &circles);
	void outputGaussBonnetPath(list<IntersectionPoint*> &points);
	void prepareCircularRegions(vector<CircularRegion> &circles, vector<CircularRegion> **newCircles);
	void insertFakeIntersectionPoints(vector<CircularRegion> &circles);
	void buildGaussBonnetPath(Vector &origin, double radius, vector<Vector> &atoms, vector<double> &radii, vector<CircularRegion> **circles, vector<list<IntersectionPoint*>*>** intersections);
	void harvestIntersectionPoints(vector<CircularRegion> &circles, vector<vec> &intersections);
	bool hasUnflaggedIntersectionPoints(CircularRegion &circle, IntersectionPoint **ip);
	list<IntersectionPoint*>* retrieveIntersections(CircularRegion &circle);
	void showIntersections(list<IntersectionPoint*> &intersections);
	vector<list<IntersectionPoint*>*>*  harvestGaussBonnetPaths(vector<CircularRegion> &circles);
	void determineCircularIntersections(vector<CircularRegion> &circles);
	double complLongAngle(Vector &vi, Vector &vj, Vector &vk);
	OcclusionType occludesForwardIntersectionPoint(CircularRegion &I, CircularRegion &J, CircularRegion &K, double dij, double dik);
	void buildIntersectionGraph(double radius, vector<CircularRegion> &circles, vector<list<IntersectionPoint*>*> &intersections);
	int sgn(double d);
	int checkIntegrityAndBlock(CircularRegion &I, CircularRegion &J);
	void setTertiaryIntersection(CircularRegion &I, CircularRegion &J, CircularRegion &K, OcclusionState state);
	int whoOccludes(CircularRegion *I, CircularRegion *J, vector<CircularRegion> &circles);
	void printSasa(list<int> &sasa);
	double intersectionDistance(CircularRegion I, CircularRegion J, CircularRegion &K, double dij, double djk);
	int selectNextIntersectingCluster(CircularRegion &I, CircularRegion &J, vector<CircularRegion> &circles);
	bool mergeCutList(int x, list<int>::iterator start, list<int> &sasa);
	
	
	
};



#endif // TESSELLATION_H_
