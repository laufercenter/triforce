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

#ifndef INTEGRATORNUMERICAL_H_
#define INTEGRATORNUMERICAL_H_

#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "tessellation.h"
#include "molecule.h"
#include "integrator.h"

#include <armadillo>


using namespace std;
using namespace arma;

class IntegratorNumerical: public Integrator{
	
public:
	IntegratorNumerical();
	IntegratorNumerical(int trials, float fd=0);
	float integrate(Molecule *molecule, Tessellation *tessellation);
	float integrate(Molecule *molecule, int index=-1,  void (*feedback)(float)=NULL);
	float integrateMolecule(Molecule *molecule, int index=-1);
	
	
private:
	Tessellation* tessellation;
	Molecule* molecule;
	int trials;
	FILE* file;
	vector<AreasDT*> areas;
	
	vector<Vector> atoms;
	vector<float> radii;
	
	float fd;
	
	
	float angle(Vector &a, Vector &b);
	Vector spherical2cartesian(Vector s);
	bool occludes(Vector v, vector<float> &lambdas, vector<Vector> &mus, vector<CircularInterfaceForm> &forms, float &conflict);
	Vector sphericalVector(float phi, float theta);
	
};

#endif //INTEGRATORNUMERICAL_H_
