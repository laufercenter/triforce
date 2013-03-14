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
#define THRESHOLD_NUMERICAL 0.000001
#define THRESHOLD_STRONG_NUMERICAL 0.000000000001
#define EPSILON 0.01
#define ORDER_CLOCKWISE 0
#define ORDER_COUNTERCLOCKWISE 1

#define FD 0.000001
#define FDT 10.0025

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

enum Orientation{
	FORWARD,
	BACKWARD
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
}
Rotation;


typedef struct
{
	double rotation;
	double di;
	double dij;
	Vector ni;
	Vector nij;
	double dot_ni_nij;
	double s0;
	double s2;
	
	int id_i;
	int id_j;
}
OmegaRotation;

typedef struct
{
	double rotation;
	double sig0;
	double sig1;
	double sig2;
	double dot_IJ;
	double rho;

	int id_i;
	int id_j;
	
}
EtaRotation;


typedef struct
{
	double rotation;
	
	OmegaRotation omega;
	EtaRotation eta;
	
	int s;
	
	
}
PHIRotation;



typedef struct
{
	PHIRotation out;
	PHIRotation in;
}
PHIContainer;



typedef struct
{
	double rotation;
	double d_i;
	double r_l;
	double r_i;
	double g;
}
LambdaRotation;



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
	PHIRotation PHI;
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
	Vector normal;
	LambdaRotation lambdaRotation;
	Rotation lambda;
	Rotation psi;
	double g;
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
	bool hasDerivatives;
		
	
}
CircularInterface;

typedef vector<CircularInterface> CircularInterfacesPerAtom;


typedef struct
{
	int id0;
	int id1;
	int id2;
	int index0;
	int index1;
	int index2;
	CircularInterfaceForm form0;
	CircularInterfaceForm form1;
	CircularInterfaceForm form2;
	Rotation rotation0;
	Rotation rotation1;
	Vector normalForCircularInterface;
	Rotation lambda;
	Rotation psi;
	Vector tessellationAxis;
	Hemisphere hemisphere;
}
SASASegment;

typedef vector<SASASegment> SASASegmentList;



typedef vector<SASASegmentList> SASAs;


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
	
	void build(bool split);
	SASAs &sasas();
	
	
	Vector calculateInterfaceNormal(const Vector &v_l, const Vector &v_i);
	Vector calculateInterfaceNormal(const Vector &v_l, const Vector &v_i, double &d);
	LambdaRotation calculateLambda(int index_i, double d_i, double r_l, double r_i, Vector &mu_i, CircularInterfaceForm &form);
	Rotation calculateLambdaDerivatives(LambdaRotation &r, CircularInterface &circle);
	LambdaRotation calculateLambda(double d_i, double r_l, double r_i, Vector &mu_i);
	Rotation calculatePsi(Vector &tessellationAxis, Vector &mu_i, Matrix &dmu_dx, CircularInterfaceForm form, int index);
	Rotation calculatePsi(Vector &tessellationAxis, Vector &mu_i);
	EtaRotation calculateEta(Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j);
	EtaRotation calculateEta(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j);
	Rotation calculateEtaDerivatives(EtaRotation &r, CircularInterfacesPerAtom &circles);
	OmegaRotation calculateOmega(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j);
	OmegaRotation calculateOmega(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j);
	Rotation calculateOmegaDerivatives(OmegaRotation &r, CircularInterfacesPerAtom &circles, Vector &tessellationAxis);
	PHIContainer calculatePHI(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j);
	PHIContainer calculatePHI(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j);
	Rotation calculatePHIDerivatives(PHIRotation &r, CircularInterfacesPerAtom &circles, Vector &tessellationAxis);
	double calculateRho(Vector &mu_i, Vector &mu_j);
	double calculateRho(Vector &mu_i, Vector &mu_j, bool derivatives);
	double sacos(Vector &a, Vector &b);
	double sdot(Vector &a, Vector &b);
	double l(Vector &a);
	void calculateProjectionAndDerivatives(Vector &tessellationAxis, CircularInterface &circle);

	
	
private:
	
	Molecule molecule;
	vector<Vector> atoms;
	vector<double> radii;
	Vector torigin;
	double tradius;
	int ti;
	//#atoms #circularregions
	//#atoms #sasas #circularregions
	SASAs sasasForMolecule;

	CircularInterfacesPerAtom coverHemisphere(Vector tessellationAxis, double radius, CircularInterfacesPerAtom circles, CircularInterfaceForm form);
	void buildGaussBonnetPath(int i, vector<Vector> &atoms, vector<double> &radii, SASAs &sasas, bool split, vector<int> &neighbourlist);
	double vsign(double v);
	double cot(double a);
	double csc(double a);
	double getAngleBetweenNormals(Vector &a, Vector &b);
	double getAngle(Vector &a, Vector &b);
	bool isZero(double v);
	bool isInPositiveEpsilonRange(double v, double eps);
	void determineProjection(Vector &origin, double radius, CircularInterface &circle);
	IntersectionPair determineIntersectionPoints(double radius, CircularInterface &K, CircularInterface &J);
	bool makeCircularInterfaces(int i,Vector &origin, double radius, vector<vec> &atoms, vector<double> &radii, vector<CircularInterface> &circles, vector<int> &neighbourlist);
	int filterCircularInterfaces(Vector tessellationAxis, double radius, vector<CircularInterface> &circles);
	void outputGaussBonnetPath(SASAs &points);
	void reindexCircularInterfaces(CircularInterfacesPerAtom &circles);
	void insertArtificialIntersectionPoints(CircularInterface &I, Vector &tessellationAxis, Hemisphere hemisphere, SASASegmentList &sasa);
	int sgn(double d);
	void determineCircularIntersections(CircularInterfacesPerAtom &circles);
	double complLongAngle(Vector &vi, Vector &vj, Vector &vk);
	IntersectionBranches::iterator increaseBranchInterator(IntersectionBranches::iterator it, CircularInterface &circle);
	IntersectionBranches::iterator decreaseBranchInterator(IntersectionBranches::iterator it, CircularInterface &circle);
	void createIntersectionBranch(PHIContainer &PHII, PHIContainer &PHIJ, CircularInterface &I, CircularInterface &J);
	void printBranch(const char* s, multimap<double, IntersectionBranch>::iterator &it);
	void printIntersectionGraph(IntersectionGraph &g, CircularInterfacesPerAtom &circles);
	void buildIntersectionGraph(double radius, Vector &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, string filename);
	void outputGaussBonnetData(string filename, double radius, CircularInterfacesPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph);
	void depleteCircularInterfaces(Vector tessellationAxis, double radius, vector<CircularInterface> &circles);
	bool isWithinNumericalLimits(double x, double l);
	bool isWithinStrongNumericalLimits(double x, double l);
	OmegaRotation calculateOmega(Vector &tessellationAxis, CircularInterface &I, CircularInterface &J);
	EtaRotation calculateEta(Vector &tessellationAxis, CircularInterface &I, CircularInterface &J);
	PHIContainer calculatePHI(Vector &tessellationAxis, CircularInterface &I, CircularInterface &J, double radius);
	void determinePsiRotations(Vector &tessellationAxis, CircularInterfacesPerAtom &circles);
	Rotation calculatePsi(Vector &tessellationAxis, CircularInterface &circle);
	Matrix matrixCross(Matrix &m, Vector &v);
	Vector normalise(Vector x);
	Vector normalise(Vector x, double &l);


	
	
};



#endif // TESSELLATION_H_
