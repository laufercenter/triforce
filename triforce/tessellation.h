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
#define THRESHOLD_NUMERICAL 0.0000001
#define EPSILON 0.01
#define ORDER_CLOCKWISE 0
#define ORDER_COUNTERCLOCKWISE 1




typedef vec Vector;

enum OcclusionState{
	OCCLUDED,
	UNOBSTRUCTED,
	UNDEFINED
};


enum Direction{
	IN,
	OUT
};





typedef struct
{
	double d;
	bool visited;
	
	
}
CircularIntersection;

typedef struct
{
	double out;
	Vector vectorOut;
	double in;
	Vector vectorIn;
}
Interfaces;

typedef struct
{
	int id0;
	int id1;
	
	int pointsTo0;
	int pointsTo1;
	
	Vector vector;
	double angle0;
	double angle1;
	
	bool visited;
}
IntersectionNode;


typedef struct
{
	int id0;
	int id1;
}
IntersectionAddress;




struct InteractionNodeComparator: public std::binary_function<IntersectionAddress, IntersectionAddress, bool>
{
	bool operator()(const IntersectionAddress& lhs, const IntersectionAddress& rhs) const
	{
		if(lhs.id0 == rhs.id0){
			return lhs.id1 < rhs.id1;
		}
		else return lhs.id0 < rhs.id0;
	}
};







struct IntersectionBranch;

typedef multimap<double, IntersectionBranch> IntersectionBranches;


typedef struct IntersectionBranch
{
	IntersectionNode* node;
	IntersectionBranches::iterator it;
	bool visited;
	Direction direction;
	IntersectionBranches* body;
	int id;
}
IntersectionBranch;





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
	map<int,CircularIntersection> circularIntersections;
	IntersectionBranches intersectionBranches;
	
	int form;
	bool intersect;
	bool flagged;
		
	
}
CircularRegion;

typedef vector<CircularRegion> CircularRegionsPerAtom;


typedef struct
{
	int id0;
	int id1;
	Vector vector;
	double angle0;
	double angle1;
	Vector normalForCircularRegion;
	double lambda;
}
SASANode;

typedef vector<SASANode> SASANodeList;

typedef struct{
	Vector tessellationOrigin;
	SASANodeList sasa;
	double radius;
	Vector vector;
}
SASA;

typedef vector<SASA> SASAs;


typedef map<IntersectionAddress, IntersectionNode, InteractionNodeComparator> IntersectionGraph;










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
	SASAs &sasas();
	
	
private:
	
	Molecule molecule;
	vector<Vector> atoms;
	vector<double> *radii;
	//#atoms #circularregions
	//#atoms #sasas #circularregions
	SASAs sasasForMolecule;


	CircularRegionsPerAtom coverHemisphere(Vector tessellationOrigin, double radius, CircularRegionsPerAtom circles);
	void buildGaussBonnetPath(Vector &origin, double radius, vector<Vector> &atoms, vector<double> &radii, SASAs &sasas);
	double vsign(double v);
	double cot(double a);
	double csc(double a);
	double getAngleBetweenNormals(Vector &a, Vector &b);
	double getAngle(Vector &a, Vector &b);
	bool isZero(double v);
	bool isInPositiveEpsilonRange(double v, double eps);
	void determineProjection(Vector &origin, double radius, CircularRegion &circle);
	IntersectionPair determineIntersectionPoints(double radius, CircularRegion &K, CircularRegion &J);
	void makeCircularRegions(Vector &origin, double radius, vector<vec> &atoms, vector<double> &radii, vector<CircularRegion> &circles);
	int filterCircularRegions(double radius, vector<CircularRegion> &circles);
	void outputGaussBonnetPath(SASA &points);
	void reindexCircularRegions(CircularRegionsPerAtom &circles);
	void insertFakeIntersectionPoints(vector<CircularRegion> &circles);
	int sgn(double d);
	void determineCircularIntersections(CircularRegionsPerAtom &circles);
	double complLongAngle(Vector &vi, Vector &vj, Vector &vk);
	double angularInterface(Vector &x0, Vector &v, Vector &p0, Vector &p1);
	void measurementPoints(Vector &p0, Vector &p1, Vector &tessellationOrigin, CircularRegion &I);
	Interfaces angularInterfaces(Vector &x0, Vector &x1, Vector &tessellationOrigin, CircularRegion &I);
	Interfaces retrieveInterfaces(Vector tessellationOrigin, CircularRegion &I, CircularRegion J, double dij, double radius);
	IntersectionBranches::iterator increaseBranchInterator(IntersectionBranches::iterator it);
	IntersectionBranches::iterator decreaseBranchInterator(IntersectionBranches::iterator it);
	IntersectionBranches::iterator increaseBranchInterator(multimap<double, IntersectionBranch>::iterator it, int ignore);
	IntersectionBranches::iterator decreaseBranchInterator(IntersectionBranches::iterator it, int ignore);
	void disconnectIntersectionPoint(IntersectionNode &a);
	void connectIntersectionPoints(IntersectionNode &a, IntersectionNode &b);
	void createIntersectionNode(IntersectionAddress &address, IntersectionGraph &intersectionGraph);
	void createIntersectionNode(int id0, int id1, IntersectionGraph &intersectionGraph);
	void createIntersectionBranch(IntersectionAddress &address, Interfaces interfacesI, Interfaces interfacesJ, CircularRegion &I, CircularRegion &J, IntersectionGraph &intersectionGraph);
	void printBranch(const char* s, multimap<double, IntersectionBranch>::iterator &it);
	void printIntersectionGraph(IntersectionGraph &g);
	void buildIntersectionGraph(double radius, Vector &tessellationOrigin, CircularRegionsPerAtom &circles, SASAs &sasas, string filename);
	void outputGaussBonnetData(string filename, CircularRegionsPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph);

	
	
};



#endif // TESSELLATION_H_
