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

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "tessellation.h"
#include "interpolation.h"
#include "molecule.h"
#include "benchmark.h"

#include <armadillo>


using namespace std;
using namespace arma;

#define THRESHOLD_NEGATIVE  0.5
#define LOGISTIC_LIMIT 0.15
#define LOGISTIC_SMOOTHER_PARAMETER 8

typedef struct
{
	float integral;
	Vector force_i;
	Vector force_j;
	Vector force_k;
	Vector force_l;
}
Integral;

typedef struct
{
	int i;
	Vector force;
}
ForceElement;





class Integrator{
public:
	virtual float integrate(Molecule *molecule, Tessellation *tessellation)=0;
};




class IntegratorTriforce: public Integrator{
	
public:
	IntegratorTriforce();
	IntegratorTriforce(vector<Interpolation*> data);
	
	virtual float integrate(Molecule *molecule, Tessellation *tessellation);
	void outputIntegrationData(string filename, Vector &integrationOrigin, list<IntersectionNode*> &frontHemisphere, list<IntersectionNode*> &backHemisphere);
	Benchmark getBenchmark();
	
	
	
private:
	vector<Interpolation*> data;
	Tessellation* tessellation;
	Molecule* molecule;
	float tradius;
	vector<float*> areas;
	vector<vector<float*> > forces;
	vector<ForceElement> forcesDelayed;
	vector<Vector> atoms;
	vector<float> radii;
	Benchmark benchmark;
	
	
	float csc(float a);
	int sgn(float d);
	float complAngle(Vector &a, Vector &b);
	float complLongAngle(Vector &n, Vector &o, Vector &a);
	float angle(Vector &a, Vector &b);
	Integral integrateTriangle(int l, SASASegment &x, Vector integrationOrigin);
	float integrateSASA(int l, SASASegmentList &s, float radius);
	Vector lookUp(float PHI, float psi, float lambda, CircularInterfaceForm form);
	Vector lookUp2(float PHI, float psi, float lambda);
	Vector lookUp3(float PHI, float psi, float lambda);
	void addForce(int i, Vector force);
	void clearForces();
	float PHI2phi(float PHI, float psi, float lambda);
	float PHI2phi2(Vector integrationOrigin, float PHI, float psi, float lambda);
	bool isInPositiveEpsilonRange(float v, float eps);
	float V2phi(Vector &integrationOrigin, Vector cv, Vector &v);
	fmat33 rotz(float theta);

	Vector recoverCircularInterface(Vector p, float psi_b, float lambda_b, float PHI_b0);
	Vector recoverCircularInterface(float psi_a, float lambda_a, float PHI_a1, float psi_b, float lambda_b, float PHI_b0);
	Vector intersectionPoint(Vector c, float psi, float lambda, float PHI);
	bool isWithinNumericalLimits(float x, float l);
	float logisticSmoother(float lambda);
	float dlogisticSmoother(float lambda);
	Vector areaSmoother(Vector &x, float area, float radius);
	float sech(float x);
	void pushForces();
	void purgeForces();
	
	
};

#endif //INTEGRATOR_H_
