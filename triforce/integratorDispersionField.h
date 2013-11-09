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

#ifndef INTEGRATORDF_H_
#define INTEGRATORDF_H_

#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "integrator.h"
#include "tessellation.h"
#include "interpolation.h"
#include "molecule.h"
#include "benchmark.h"

#include <armadillo>


using namespace std;
using namespace arma;

#define DIM_DISPERSIONFIELD 7


enum ParameterTypes{
	EPSILON,
	SIGMA
};



typedef struct
{
	Vector eps;
	Vector sig;
}
DispersionElement;


typedef struct
{
	Integral eps;
	Integral sig;
}
DispersionIntegral;

typedef struct
{
	float eps;
	float sig;
	float area;
}
Dispersion;

class IntegratorDispersionField: public IntegratorTriforce{
	
public:
	
	float integrate(Molecule *molecule, Tessellation *tessellation);
	
	
	
private:
	//(eps,sig) (6) (front,back) (data+6*grad)
	Interpolation* dataSolvationFreeEnergy;
	Interpolation* dataDispersion[2][N_SPECIES][2][DIM_DISPERSIONFIELD];
	vector<unsigned int> species;
	vector<float*> nonpolar;
	vector<float> epsilons;
	vector<float> sigmas;
	
	
	
	Integral integrateTriangle(int l, SASASegment &x, Vector integrationOrigin, float &phi);
	float integrateSASA(int l, SASASegmentList &s, float radius);
	DispersionElement lookUp(float PHI, float psi, float lambda, float &phi, CircularInterfaceForm form);

	
};

#endif //INTEGRATORDF_H_
