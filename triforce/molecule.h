/* Copyright 2012, Nils J. D. Drechsel & Jordi Villà-Freixa
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


	

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <string>
#include <vector>
#include <map>

#include <armadillo>

#include "topology.h"
#include "neighbourList.h"


#define N_SPECIES 6


using namespace std;
using namespace arma;

typedef fvec Vector;


enum ForceField {
	Amber99SBildn,
	Amber99SB,
	Amber99,
	Amber96
};

enum Source {
	EXTERNAL_SOURCE,
	INTERNAL_SOURCE
};





#ifndef DOUBLE_COORDINATE_POINTERS
	typedef float CoordinatesDT;
	typedef fvec VectorCoordinates;
#else
	typedef double CoordinatesDT; 
	typedef vec VectorCoordinates;
#endif
	
#ifndef DOUBLE_AREA_POINTERS
	typedef float AreasDT;
#else
	typedef double AreasDT; 
#endif

#ifndef DOUBLE_COORDINATE_POINTERS
	typedef float ForcesDT;
#else
	typedef double ForcesDT; 
#endif

	
typedef struct
{
	CoordinatesDT* x;
	CoordinatesDT* y;
	CoordinatesDT* z;
}
AtomicPointers; 



class Molecule{
	
public:
	Molecule();
	Molecule(Topology topology);
	void init();
	
	void addAtom(CoordinatesDT* x, CoordinatesDT* y, CoordinatesDT* z, AreasDT* area, ForcesDT* forceX, ForcesDT* forceY, ForcesDT* forceZ, string name, string type, int i=-1);
	void addAtom(CoordinatesDT* x, CoordinatesDT* y, CoordinatesDT* z, AreasDT* area, ForcesDT* forceX, ForcesDT* forceY, ForcesDT* forceZ, float radius, string name, int i=-1);
	void addAtom(CoordinatesDT* x, CoordinatesDT* y, CoordinatesDT* z, AreasDT* area, ForcesDT* forceX, ForcesDT* forceY, ForcesDT* forceZ, float eps, float sig, string name, int i=-1);
	void addAtom(CoordinatesDT* x, CoordinatesDT* y, CoordinatesDT* z, AreasDT* area, ForcesDT* forceX, ForcesDT* forceY, ForcesDT* forceZ, float radius, float eps, float sig, string name, int i=-1);
	void update();
	void updateAtomicPositions();
	void addInternallyStoredAtom(float x, float y, float z, string name, string type, int i=-1);
	void addInternallyStoredAtom(float x, float y, float z, float radius, string name, int i=-1);
	void addInternallyStoredAtom(float x, float y, float z, float eps, float sig, string name, int i=-1);
	void addInternallyStoredAtom(float x, float y, float z, float radius, float eps, float sig, string name, int i=-1);
	Vector getInternallyStoredAtomCoordinates(int i);
	void setInternallyStoredAtomCoordinates(int i, Vector &v);
	void perturbInternallyStoredAtomCoordinates(int i, Vector p);
	void jiggleInternallyStoredAtomCoordinates(float mag);
	Vector getCenter();
	
	void print(FILE* outputfile);
	void printxyz(FILE* outputfile);
	void printxyzr(FILE* outputfile);
	
	void printDifference(Molecule *mol);
	
	vector<Vector> &fetchCoordinates();
	vector<float> &fetchRadii();
	vector<AreasDT*> &fetchAreaPointers();
	vector<vector<ForcesDT*> > &fetchForcePointers();
	vector<float> &fetchEpsilons();
	vector<float> &fetchSigmas();
	vector<float*> &fetchNonpolarPointers();
	vector<unsigned int> &fetchSpecies();
	void generateNeighbourList();
	vector<int> getNeighborListFor(int i);
	vector<int> &fetchClosestNeighbours();
	void calculateClosestNeighbours();
	void refreshInternalPointers();
	
	
private:
	
	string string2UpperCase(string s);
	void constructAtoms(unsigned int end);
	

	
	Topology topology;
	vector<AtomicPointers> atomicPointers;
	vector<Vector> atoms;
	vector<float> epsilons;
	vector<float> sigmas;
	vector<float> radii;
	vector<unsigned int> species;
	vector<AreasDT*> areas;
	vector<vector<ForcesDT*> >forces;
	vector<float*> nonpolarFreeEnergies;
	vector<Source> source;
	vector<string> names;
	bool hasNeighbourList;
	vector<int> closestNeighbours;
	
	ForceField forcefield;
	
	
	//this will only be needed for usage with the tool. Amber and Adun will provide their own structures
	vector<AreasDT> realArea;
	vector<VectorCoordinates> realAtoms;
	vector<ForcesDT> realForceX;
	vector<ForcesDT> realForceY;
	vector<ForcesDT> realForceZ;
	

	NeighbourList *neighbourList;
	
	unsigned int identifySpecies(float eps, float sig);
	
	
};

#endif //MOLECULE_H_
