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



class Integrator{
public:
	virtual double integrate(Molecule *molecule, Tessellation *tessellation)=0;
};




class IntegratorTriforce: public Integrator{
	
public:
	IntegratorTriforce();
	IntegratorTriforce(Interpolation *dataConvex, Interpolation *dataConcave);
	
	double integrate(Molecule *molecule, Tessellation *tessellation);
	void outputIntegrationData(string filename, Vector &integrationOrigin, list<IntersectionNode*> &frontHemisphere, list<IntersectionNode*> &backHemisphere);
	

	
private:
	Interpolation *dataConvex;
	vector<Interpolation*> forcesConvex;
	Interpolation *dataConcave;
	vector<Interpolation*> forcesConcave;
	
	
	Tessellation* tessellation;
	Molecule* molecule;
	vector<double*> areas;
	vector<vector<double*> > forces;
	
	//void splitSASA(list<IntersectionNode*> &sasa, vector<CircularRegion> &circles, int c, Vector &integrationOrigin, double radius, list<IntersectionNode*>** frontHemisphere, list<IntersectionNode*>** backHemisphere,  IntersectionGraph &intersectionGraph);
	//Vector halfSphereIntersectionPoint(Vector &integrationOrigin, CircularRegion &c, double radius, int sign);
	double csc(double a);
	int sgn(double d);
	double complAngle(Vector &a, Vector &b);
	double complLongAngle(Vector &n, Vector &o, Vector &a);
	double angle(Vector &a, Vector &b);
	double integrateTriangle(SASANode &x0, SASANode &x1, Vector integrationOrigin, double &totalAngle);
	//double integrateHemisphere(list<IntersectionNode*> &sasa, Vector &integrationOrigin, vector<CircularRegion> &circles, int ci);
	double integrateAtomicSASA(SASAsForAtom sasasForAtom);
	double integrateSASA(SASA &s);
	double PHI2phi(double PHI, double psi, double lambda);
	double V2phi(Vector &integrationOrigin, Vector cv, Vector &v);
	bool isInPositiveEpsilonRange(double v, double eps);
	mat33 rotz(double theta);
	double PHI2phi2(Vector integrationOrigin, double PHI, double psi, double lambda);
	
	
};

#endif //INTEGRATOR_H_
