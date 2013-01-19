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

#include <armadillo>


using namespace std;
using namespace arma;

typedef struct
{
	double area;
	Vector force_i;
	Vector force_j;
	Vector force_k;
	Vector force_l;
}
Area;




class Integrator{
public:
	virtual double integrate(Molecule *molecule, Tessellation *tessellation)=0;
};




class IntegratorTriforce: public Integrator{
	
public:
	IntegratorTriforce();
	IntegratorTriforce(Interpolation *dataConcave, Interpolation *forcesConcave1, Interpolation *forcesConcave2, Interpolation *forcesConcave3);
	
	double integrate(Molecule *molecule, Tessellation *tessellation);
	void outputIntegrationData(string filename, Vector &integrationOrigin, list<IntersectionNode*> &frontHemisphere, list<IntersectionNode*> &backHemisphere);
	
	
	
private:
	vector<Interpolation*> data;
	Tessellation* tessellation;
	Molecule* molecule;
	
	vector<double*> areas;
	vector<vector<double*> > forces;	
	vector<Vector> atoms;
	vector<double> radii;
	
	
	double csc(double a);
	int sgn(double d);
	double complAngle(Vector &a, Vector &b);
	double complLongAngle(Vector &n, Vector &o, Vector &a);
	double angle(Vector &a, Vector &b);
	Area integrateTriangle(int l, SASANode &x0, SASANode &x1, Vector integrationOrigin);
	double integrateAtomicSASA(int l, SASAsForAtom sasasForAtom, double radius);
	double integrateSASA(int l, SASA &s, double radius);
	Vector lookUp(double PHI, double psi, double lambda, double q, CircularInterfaceForm form);
	Vector lookUp2(double PHI, double psi, double lambda);
	Vector lookUp3(double PHI, double psi, double lambda);
	void addForce(int i, Vector force);
	void clearForces();
	double PHI2phi(double PHI, double psi, double lambda);
	double PHI2phi2(Vector integrationOrigin, double PHI, double psi, double lambda);
	bool isInPositiveEpsilonRange(double v, double eps);
	double V2phi(Vector &integrationOrigin, Vector cv, Vector &v);
	mat33 rotz(double theta);
	
};

#endif //INTEGRATOR_H_
