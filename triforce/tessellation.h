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
#include "interpolation.h"



using namespace std;
using namespace arma;



#define ANGLE_TRESHOLD 0.001
#define THRESHOLD_IP 0.05
#define THRESHOLD_NUMERICAL 0.0000001
#define EPSILON 0.01
#define ORDER_CLOCKWISE 0
#define ORDER_COUNTERCLOCKWISE 1

#define FD 0.0000001
#define FDT 0.1

#define MINISCULE 0.00001



typedef vec Vector;
typedef Col<int> VectorInt;

enum OcclusionState{
	OCCLUDED,
	UNOBSTRUCTED,
	UNDEFINED
};


enum CircularInterfaceForm{
	CONCAVE,
	CONVEX,
	SPLITTER
};



enum Location{
	INTERNAL,
	EXTERNAL
};

enum Direction{
	IN,
	OUT
};

enum Hemisphere{
	FRONTHEMISPHERE,
	BACKHEMISPHERE
};


typedef struct
{
	double rotation;
	Vector drotation_dxi;
	Vector drotation_dxj;
	Vector drotation_dxl;
	Vector vector;
}
Rotation;


typedef struct
{
	Rotation out;
	Rotation in;
}
PHIContainer;



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
	
	int prev0;
	int prev1;
	
	Vector vector;
	Rotation rotation0;
	Rotation rotation1;
	
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
	int visited;
	Direction direction;
	Location forward;
	Location backward;
	IntersectionBranches* body;
	int id;
	bool flagged;
}
IntersectionBranch;






struct IteratorComparator: public std::binary_function<IntersectionBranches::iterator, IntersectionBranches::iterator, bool>
{
	bool operator()(const IntersectionBranches::iterator& lhs, const IntersectionBranches::iterator& rhs) const
	{
		if(lhs->second.node->id0 == rhs->second.node->id0){
			return lhs->second.node->id1 < rhs->second.node->id1;
		}
		else return lhs->second.node->id0 < rhs->second.node->id0;
	}
};




typedef struct 
{
	int id;
	int index;
	Vector vector;
	Vector normal;
	Rotation lambda;
	Rotation psi;
	double g;
	double g_normalised;
	double a;
	double sphereRadius;
	double d;
//	double radius;
	Matrix dmu_dx;
	map<int,CircularIntersection> circularIntersections;
	IntersectionBranches intersectionBranches;
	
	CircularInterfaceForm form;
	bool intersect;
	bool flagged;
	bool valid;
		
	
}
CircularInterface;

typedef vector<CircularInterface> CircularInterfacesPerAtom;


typedef struct
{
	int id0;
	int id1;
	int index0;
	int index1;
	Vector vector;
	Rotation rotation0;
	Rotation rotation1;
	Vector normalForCircularInterface;
	Rotation lambda;
	Rotation psi;
	CircularInterfaceForm form;
}
SASANode;

typedef vector<SASANode> SASANodeList;


typedef struct{
	Vector tessellationOrigin;
	Hemisphere hemisphere;
	SASANodeList sasa;
}
SASA;

typedef vector<SASA> SASAs;

typedef struct{
	SASAs sasas;
	double radius;
	Vector vector;
}
SASAsForAtom;

typedef vector<SASAsForAtom> SASAsForMolecule;


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
	Tessellation(Molecule &m, Interpolation *dataPHI);
	
	void build(bool split);
	SASAsForMolecule &sasas();
	
	
private:
	
	Molecule molecule;
	vector<Vector> atoms;
	vector<double> radii;
	Vector torigin;
	double tradius;
	//#atoms #circularregions
	//#atoms #sasas #circularregions
	SASAsForMolecule sasasForMolecule;

	CircularInterfacesPerAtom coverHemisphere(Vector tessellationOrigin, double radius, CircularInterfacesPerAtom circles, CircularInterfaceForm form);
	void buildGaussBonnetPath(int i, vector<Vector> &atoms, vector<double> &radii, SASAsForMolecule &sasas, bool split);
	double vsign(double v);
	double cot(double a);
	double csc(double a);
	double getAngleBetweenNormals(Vector &a, Vector &b);
	double getAngle(Vector &a, Vector &b);
	bool isZero(double v);
	bool isInPositiveEpsilonRange(double v, double eps);
	void determineProjection(Vector &origin, double radius, CircularInterface &circle);
	IntersectionPair determineIntersectionPoints(double radius, CircularInterface &K, CircularInterface &J);
	void makeCircularInterfaces(int i,Vector &origin, double radius, vector<vec> &atoms, vector<double> &radii, vector<CircularInterface> &circles);
	int filterCircularInterfaces(Vector tessellationOrigin, double radius, vector<CircularInterface> &circles);
	void outputGaussBonnetPath(SASA &points);
	void reindexCircularInterfaces(CircularInterfacesPerAtom &circles);
	void insertArtificialIntersectionPoints(CircularInterface &I, IntersectionGraph &intersectionGraph, Vector &tessellationOrigin);
	int sgn(double d);
	void determineCircularIntersections(CircularInterfacesPerAtom &circles);
	double complLongAngle(Vector &vi, Vector &vj, Vector &vk);
	double angularInterface(Vector &x0, Vector &v, Vector &p0, Vector &p1);
	void measurementPoints(Vector &p0, Vector &p1, Vector &tessellationOrigin, CircularInterface &I);
	Interfaces angularInterfaces(Vector &x0, Vector &x1, Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J);
	PHIContainer retrieveInterfaces(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J, double dij, double radius);
	IntersectionBranches::iterator increaseBranchInterator(IntersectionBranches::iterator it, CircularInterfacesPerAtom &circles);
	IntersectionBranches::iterator decreaseBranchInterator(IntersectionBranches::iterator it, CircularInterfacesPerAtom &circles);
	IntersectionBranches::iterator increaseBranchInterator(multimap<double, IntersectionBranch>::iterator it, int ignore, CircularInterfacesPerAtom &circles);
	IntersectionBranches::iterator decreaseBranchInterator(IntersectionBranches::iterator it, int ignore, CircularInterfacesPerAtom &circles);
	void disconnectIntersectionPoint(IntersectionNode &a);
	void connectIntersectionPoints(IntersectionNode &a, IntersectionNode &b, IntersectionGraph &intersectionGraph);
	void createIntersectionNode(IntersectionAddress &address, IntersectionGraph &intersectionGraph);
	void createIntersectionNode(int id0, int id1, IntersectionGraph &intersectionGraph);
	void createIntersectionBranch(IntersectionAddress &address, PHIContainer &PHII, PHIContainer &PHIJ, CircularInterface &I, CircularInterface &J, IntersectionGraph &intersectionGraph);
	void printBranch(const char* s, multimap<double, IntersectionBranch>::iterator &it);
	void printIntersectionGraph(IntersectionGraph &g);
	void buildIntersectionGraph(double radius, Vector &tessellationOrigin, CircularInterfacesPerAtom &circles, SASAs &sasas, Hemisphere hemisphere, string filename);
	void outputGaussBonnetData(string filename, double radius, CircularInterfacesPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph);
	void deleteIntersectionPoint(IntersectionBranches::iterator &it,IntersectionGraph &intersectionGraph, CircularInterfacesPerAtom &circles);
	void depleteCircularInterfaces(Vector tessellationOrigin, double radius, vector<CircularInterface> &circles);
	double rotationalAngle(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J);
	bool isWithinNumericalLimits(double x, double l);
	Rotation calculateOmega(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J);
	Rotation calculateEta(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J);
	PHIContainer calculatePHI(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J, double dij, double radius);
	void determinePsiRotations(Vector &tessellationOrigin, CircularInterfacesPerAtom &circles);
	Rotation calculatePsi(Vector tessellationOrigin, CircularInterface &circle);
	Matrix matrixCross(Matrix &m, Vector &v);
	bool addToEraseList(map<IntersectionBranches::iterator,bool,IteratorComparator> &masterEraseList, map<IntersectionBranches::iterator,bool,IteratorComparator> &eraseList, IntersectionBranches::iterator &it, int limit);
	bool addToEraseListCascade(map<IntersectionBranches::iterator,bool,IteratorComparator> &masterEraseList, map<IntersectionBranches::iterator,bool,IteratorComparator> &eraseList, IntersectionBranches::iterator &it, int limit, CircularInterfacesPerAtom &circles);


	
	
};



#endif // TESSELLATION_H_
