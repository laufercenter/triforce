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
#include "depth3d.h"
#include "multiLayeredDepthBuffer.h"
#include "benchmark.h"


using namespace std;
using namespace arma;



#define ANGLE_TRESHOLD 0.001
#define THRESHOLD_IP 0.05
#define THRESHOLD_NUMERICAL 0.000001
#define THRESHOLD_STRONG_NUMERICAL 0.000000000001
#define EPSILON 0.01
#define ORDER_CLOCKWISE 0
#define ORDER_COUNTERCLOCKWISE 1
#define THRESHOLD_INTERFACE 0.0001
#define THRESHOLD_GENERAL_POSITION 0.025


#define FD 0.001
#define FDT 0.01

#define MINISCULE 0.0000001



typedef fvec Vector;
typedef Col<unsigned int> VectorInt;

enum DerivativeMode{
	FIXED_POSITION,
	GENERAL_POSITION,
	ALIGNED
};

enum BranchMode{
	BRANCH_ALL,
	BRANCH_OUT,
	BRANCH_IN
};

enum OcclusionState{
	OCCLUDED,
	UNOBSTRUCTED,
	UNDEFINED
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


enum CircularInterfaceForm{
	CONCAVE,
	CONVEX,
	SPLITTER
};




struct SegmentInfo;

typedef struct{
	int i0;
	int i1;
}
PartialSegmentID;

typedef struct{
	int i0;
	int i1;
	int i2;
}
FullSegmentID;

typedef struct{
	int index;
	Vector v;
	Vector auxiliary;
	Vector planeNormal;
	Matrix dchi_dx;
	Hemisphere hemisphere;
	DerivativeMode mode;
}
TessellationAxis;


struct SegmentComparator: public std::binary_function<PartialSegmentID, PartialSegmentID, bool>
{
	bool operator()(const PartialSegmentID& lhs, const PartialSegmentID& rhs) const
	{
		if(lhs.i0 == rhs.i0){
			return lhs.i1 < rhs.i1;
		}
		else return lhs.i0 < rhs.i0;
	}
};

struct SegmentSetComparator: public std::binary_function<FullSegmentID, FullSegmentID, bool>
{
	bool operator()(const FullSegmentID& lhs, const FullSegmentID& rhs) const
	{
		if(lhs.i0 == rhs.i0){
			if(lhs.i1==rhs.i1){
				return lhs.i2 < rhs.i2;
			}
			else{
				return lhs.i1 < rhs.i1;
			}
		}
		else return lhs.i0 < rhs.i0;
	}
};

typedef struct{
	int id;
	Direction direction;
}
SplitterIntersection;

struct SplitterIntersectionComparator: public std::binary_function<SplitterIntersection, SplitterIntersection, bool>
{
	bool operator()(const SplitterIntersection& lhs, const SplitterIntersection& rhs) const
	{
		if(lhs.id == rhs.id){
			return lhs.direction < rhs.direction;
		}
		else return lhs.id < rhs.id;
	}
};



typedef struct
{
	float rotation;
	Vector drotation_dxi;
	Vector drotation_dxj;
	Vector drotation_dxl;
	Vector drotation_dxt;
}
Rotation;


typedef struct
{
	float rotation;
	float di;
	float dij;
	Vector ni;
	Vector nij;
	Vector vi;
	Vector vij;
	float dot_ni_nij;
	float s0;
	float s2;
	
	int id_i;
	int id_j;
}
OmegaRotation;

typedef struct
{
	float rotation;
	float sig0;
	float sig1;
	float sig2;
	float dot_IJ;
	float rho;

	int id_i;
	int id_j;
	
}
EtaRotation;


typedef struct
{
	float rotation;
	
	OmegaRotation omega;
	EtaRotation eta;
	
	int s;
	
	
}
PHIRotation;



struct IntersectionBranch;


typedef vector<SegmentInfo> SegmentList;

typedef map<PartialSegmentID, SegmentList::iterator, SegmentComparator> SegmentGraph;


typedef struct SegmentInfo{
	PartialSegmentID backw;
	PartialSegmentID id0;
	PartialSegmentID id1;
	PartialSegmentID forw;
	SegmentList::iterator it_backw;
	SegmentList::iterator it_forw;
	bool hasForward;
	bool hasBackward;
	unsigned int visited;
	float dist;
	PHIRotation PHI0;
	PHIRotation PHI1;
	bool extended;
	Vector v0;
	Vector v1;
	float weight;
	float weight1;
	IntersectionBranch* source;
	int i;
} SegmentInfo;




typedef struct
{
	PHIRotation out;
	PHIRotation in;
}
PHIContainer;



typedef struct

{
	float rotation;
	float d_i;
	float r_l;
	float r_i;
	float g;
}
LambdaRotation;



typedef struct
{
	float d;
	bool visited;
	
	
}
CircularIntersection;

typedef struct
{
	float out;
	Vector vectorOut;
	float in;
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





typedef struct
{
	int id;
	float rho;
	float dot_mui_muj;
}
RhoContainer;




typedef multimap<float, IntersectionBranch*> IntersectionBranches;


typedef struct IntersectionBranch{
	IntersectionNode* node;
	IntersectionBranches::iterator it;
	int visited;
	int visited1;
	Direction direction;
	Location forward;
	Location backward;
	IntersectionBranches* body;
	int id;
	bool flagged;
	PHIRotation PHI;
	RhoContainer rho;
	Vector v0;
	Vector v1;
	float weight;
	float weight1;
	int i;
}
IntersectionBranch;






struct IteratorComparator: public std::binary_function<IntersectionBranches::iterator, IntersectionBranches::iterator, bool>
{
	bool operator()(const IntersectionBranches::iterator& lhs, const IntersectionBranches::iterator& rhs) const
	{
		if(lhs->second->node->id0 == rhs->second->node->id0){
			return lhs->second->node->id1 < rhs->second->node->id1;
		}
		else return lhs->second->node->id0 < rhs->second->node->id0;
	}
};

typedef map<int,RhoContainer> CircularIntersections;

typedef struct CircularInterface
{
	int id;
	int index;
	Vector normal;
	LambdaRotation lambdaRotation;
	Rotation lambda;
	Rotation psi;
	float kappa[2];
	float psi2[2];
	float g;
	float sphereRadius;
	float d;
	Matrix dmu_dx;
	CircularIntersections circularIntersections;
	IntersectionBranches intersectionBranches;
	
	CircularInterfaceForm form;
	bool intersect;
	bool flagged;
	bool valid;
	bool hasDerivatives;
	bool extended;
	vector<Vector> exposedVectors;
	vector<float> exposedPHI[2];
	Vector base[2];
	Vector center;
	bool hasBase[2];
	bool hasCenter;
	bool erased;
	Hemisphere hemisphere;
	unsigned int link;
}
CircularInterface;



typedef vector<CircularInterface> CircularInterfacesPerAtom;


typedef struct
{
	int i;
	int i2;
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
	TessellationAxis tessellationAxis;
	Hemisphere hemisphere;
	float depthBufferEstimatedArea;
	float radius;
	Vector v0;
	Vector v1;
	float weight;
	float weight1;
	float kappa;
	unsigned int link;
	int htype;
	PHIRotation PHI0;
	PHIRotation PHI1;
	bool artificial;
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
	Tessellation(Molecule &m, unsigned int numbBuffer, Depth3D &depthData, Data1D &occludedDistribution, Data1D &exposedDistribution);
	
	void build(bool useDepthBuffer, bool split=true);
	SASAs &sasas();
	void update();
	
	
	Vector calculateInterfaceNormal(const Vector &v_l, const Vector &v_i);
	Vector calculateInterfaceNormal(const Vector &v_l, const Vector &v_i, float &d);
	LambdaRotation calculateLambda(int index_i, float d_i, float r_l, float r_i, Vector &mu_i, CircularInterfaceForm &form);
	Rotation calculateLambdaDerivatives(LambdaRotation &r, CircularInterface &circle);
	LambdaRotation calculateLambda(float d_i, float r_l, float r_i, Vector &mu_i);
	Rotation calculatePsi(TessellationAxis &tessellationAxis, Vector &mu_i, Matrix &dmu_dx, CircularInterfaceForm form, int index);
	Rotation calculatePsi(TessellationAxis &tessellationAxis, Vector &mu_i);
	EtaRotation calculateEta(Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, RhoContainer &rhoContainer);
	EtaRotation calculateEta(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer);
	Rotation calculateEtaDerivatives(EtaRotation &r, CircularInterfacesPerAtom &circles);
	OmegaRotation calculateOmega(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, RhoContainer &rhoContainer, int index_i);
	OmegaRotation calculateOmega(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer, int index_i);
	Rotation calculateOmegaDerivatives(OmegaRotation &r, CircularInterfacesPerAtom &circles, TessellationAxis &tessellationAxis);
	PHIContainer calculatePHI(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, RhoContainer &rhoContainer, int index_i);
	PHIContainer calculatePHI(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer, int index_i);
	Rotation calculatePHIDerivatives(PHIRotation &r, CircularInterfacesPerAtom &circles, TessellationAxis &tessellationAxis);
	float calculateRho(Vector &mu_i, Vector &mu_j);
	float calculateRho(Vector &mu_i, Vector &mu_j, bool derivatives);
	float sacos(Vector &a, Vector &b);
	float sdot(Vector &a, Vector &b);
	float l(Vector &a);
	void calculateProjectionAndDerivatives(TessellationAxis &tessellationAxis, CircularInterface &circle);
	Benchmark getBenchmark();
	void outputTessellation(string filename);
// 	Vector generalPosition(CircularInterfacesPerAtom &circles);
	
	float V2PHI(TessellationAxis &tessellationAxis, Hemisphere hemisphere, Vector v, Vector &normal, float g);
	float cot(float a);
	float csc(float a);
	float getAngleBetweenNormals(Vector &a, Vector &b);
	float getAngle(Vector &a, Vector &b);
	bool isZero(float v);
	float vsign(float v);
	bool isInPositiveEpsilonRange(float v, float eps);
	bool isWithinNumericalLimits(float x, float l);
	bool isWithinStrongNumericalLimits(float x, float l);
	int sgn(float d);

	

	Vector generalPosition(CircularInterfacesPerAtom &circles);
	void print(FILE* outputfile);
	float calculateKappa(TessellationAxis tessellationAxis, Vector &v, CircularInterfaceForm form);
	fmat33 rotz(float theta);
	fmat33 roty(float theta);
	fmat33 rotx(float theta);
	Vector PHI2V(TessellationAxis tessellationAxis, float PHI, float psi, float lambda, float kappa, CircularInterfaceForm form);

	
private:
	
	Depth3D depthData;
	Data1D occludedDistribution;
	Data1D exposedDistribution;
	Molecule molecule;
	vector<Vector> atoms;
	vector<float> radii;
	Vector torigin;
	float tradius;
	int ti;
	CircularInterfacesPerAtom tcircles;
	Vector ttessellationAxis;
	//#atoms #circularregions
	//#atoms #sasas #circularregions
	SASAs sasasForMolecule;
	SASAs prevSasasForMolecule;
	bool hasDepthBuffer;
	unsigned int numbBuffer;
	Benchmark benchmark;
	float depthBufferEstimatedArea;
	int totalInterfaces;
	int survivedInterfaces;


	void setupTessellationAxis(TessellationAxis &tessellationAxis, Hemisphere hemisphere, int closestNeighbour, vector<Vector> &atoms, vector<float> &radii, Vector &origin, float radius, DerivativeMode dmode, CircularInterfacesPerAtom &circles);
	CircularInterfacesPerAtom coverHemisphere(TessellationAxis &tessellationAxis, CircularInterfacesPerAtom circles, CircularInterfaceForm form);
	void buildGaussBonnetPath(int i, vector<Vector> &atoms, vector<float> &radii, SASAs &sasas, bool split, vector<int> &neighbourlist, int closestNeighbour, bool useDepthBuffer);
	void determineProjection(Vector &origin, float radius, CircularInterface &circle);
	void makeCircularInterfaces(int i,Vector &origin, float radius, vector<fvec> &atoms, vector<float> &radii, vector<CircularInterface> &circles, vector<int> &neighbourlist, bool &isInsideAnotherSphere);
	int filterCircularInterfaces(vector<CircularInterface> &circles, bool splitterOnly);
	void outputGaussBonnetPath(SASAs &points);
	void reindexCircularInterfaces(CircularInterfacesPerAtom &circles);
	void insertArtificialIntersectionPoints(CircularInterface &I, TessellationAxis &tessellationAxis, Hemisphere hemisphere, SASASegmentList &sasa, unsigned int &globalSegmentCounter);
	void determineCircularIntersections(CircularInterfacesPerAtom &circles, bool splitterOnly);
	IntersectionBranches::iterator increaseBranchInterator(IntersectionBranches::iterator it, CircularInterface &circle);
	IntersectionBranches::iterator decreaseBranchInterator(IntersectionBranches::iterator it, CircularInterface &circle);
	void createIntersectionBranch(PHIContainer &PHII, CircularInterface &I, CircularInterface &J, RhoContainer &rho);
	void createIntersectionBranch(PHIContainer &PHII, CircularInterface &I, CircularInterface &J, RhoContainer &rho, IntersectionBranch &b);
	

	void printBranch(const char* s, multimap<float, IntersectionBranch>::iterator &it);
	void printIntersectionGraph(IntersectionGraph &g, CircularInterfacesPerAtom &circles);
	
	
	//bool buildIntersectionGraph(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, string filename, MultiLayeredDepthBuffer &buffer0, MultiLayeredDepthBuffer &buffer1, bool useDepthBuffer, bool split, unsigned int &globalSegmentCounter, bool derivatives);
	
	void buildIntersectionGraphFirstPass(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles);
	void buildIntersectionGraphSplitterPass(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles);
	void splitterSanityCheck(CircularInterfacesPerAtom &circles);
	void copyIntersectionGraph(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, CircularInterfacesPerAtom &newCircles);
	void buildIntersectionGraphArtificialPointsPass(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, unsigned int &globalSegmentCounter);
	bool buildIntersectionGraphCollectionPass(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, string filename, MultiLayeredDepthBuffer &buffer0, MultiLayeredDepthBuffer &buffer1, bool useDepthBuffer, bool split, unsigned int &globalSegmentCounter, bool derivatives);
	
	
	void outputGaussBonnetData(string filename, float radius, CircularInterfacesPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph);
	OmegaRotation calculateOmega(TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J, RhoContainer &rhoContainer);
	EtaRotation calculateEta(TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J, RhoContainer &rhoContainer);
	PHIContainer calculatePHI(TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J, float radius, RhoContainer &rhoContainer);
	void determinePsiRotations(TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles);
	Rotation calculatePsi(TessellationAxis &tessellationAxis, CircularInterface &circle);
	Matrix matrixCross(Matrix &m, Vector &v);
	Vector normalise(Vector x);
	Vector normalise(Vector x, float &l);
	void convertExposedVectors2PHIValues(TessellationAxis &tessellationAxis,Hemisphere hemisphere, CircularInterface &circle);
	float exposition(Hemisphere hemisphere, IntersectionBranches::iterator it0, IntersectionBranches::iterator it1, CircularInterface &circle);
	//void splitTessellation(SASASegmentList &sasa, TessellationAxis &frontTessellationAxis, TessellationAxis &backTessellationAxis, CircularInterfacesPerAtom &circles);
	unsigned int coverHemisphere2(TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, CircularInterfaceForm form);
	
	void reindexCircularInterfaces2(CircularInterfacesPerAtom &circles);
	void reindexCircularInterfaces3(CircularInterfacesPerAtom &circles);
	void addPsiAndLambda(TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles);
	void sortGaussBonnetPaths(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, string filename, MultiLayeredDepthBuffer &buffer0, MultiLayeredDepthBuffer &buffer1, bool useDepthBuffer, bool split, unsigned int &globalSegmentCounter, bool derivatives);
	
	void cleanCircularIntersections(CircularInterfacesPerAtom &circles);
	
	
	
	
	
	
	
	
	
	
	
	


	
	
};



#endif // TESSELLATION_H_
