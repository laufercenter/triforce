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


using namespace std;
using namespace arma;

typedef vec Vector;

typedef struct
{
	double* x;
	double* y;
	double* z;
}
AtomicPointers; 


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



class Molecule{
	
public:
	Molecule();
	Molecule(Topology topology);
	void addAtom(double* x, double* y, double* z, double* area, double* forceX, double* forceY, double* forceZ, string type, int i=-1);
	void addAtom(double* x, double* y, double* z, double* area, double* forceX, double* forceY, double* forceZ, double sigma, double epsilon, int i=-1);
	void update();
	void addInternallyStoredAtom(double x, double y, double z, string type, int i=-1);
	void addInternallyStoredAtom(double x, double y, double z, double sigma, double epsilon, int i=-1);
	Vector getInternallyStoredAtomCoordinates(int i);
	void setInternallyStoredAtomCoordinates(int i, Vector &v);
	void perturbInternallyStoredAtomCoordinates(int i, Vector p);
	
	void print();
	vector<Vector> &fetchCoordinates();
	vector<double> &fetchRadii();
	vector<double*> &fetchAreaPointers();
	vector<vector<double*> > &fetchForcePointers();
	
private:
	
	string string2UpperCase(string s);
	void constructAtoms(int end);
	

	
	Topology topology;
	vector<AtomicPointers> atomicPointers;
	vector<Vector> atoms;
	vector<double> sigmas;
	vector<double> radii;
	vector<double> epsilons;
	vector<double*> areas;
	vector<vector<double*> >forces;
	vector<Source> source;
	
	
	ForceField forcefield;
	
	
	//this will only be needed for usage with the tool. Amber and Adun will provide their own structures
	vector<double> realArea;
	vector<double> realForceX;
	vector<double> realForceY;
	vector<double> realForceZ;
	

	
	void refreshInternalPointers();
	
	
};

#endif //MOLECULE_H_
