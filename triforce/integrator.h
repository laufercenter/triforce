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

#include <armadillo>


using namespace std;
using namespace arma;

class Integrator{
	
public:
	Integrator();
	Integrator(Tessellation *tessellation);
	Integrator(Tessellation *tessellation, Interpolation *interpolation);
	
	double Integrate();
	
	
	
private:
	Interpolation *interpolation;
	Tessellation* tessellation;
	
	void splitSASA(list<IntersectionPoint*> &sasa, vector<CircularRegion> &circles, int c, Vector &integrationOrigin, list<IntersectionPoint*>** frontHemisphere, list<IntersectionPoint*>** backHemisphere );
	Vector halfSphereIntersectionPoint(Vector &integrationOrigin, CircularRegion &c);
	double csc(double a);
	double complLongAngle(Vector &n, Vector &o, Vector &a);
	int sgn(double d);
	double complAngle(Vector &a, Vector &b);
	double angle(Vector &a, Vector &b);
	double integrateTriangle(IntersectionPoint &x0, IntersectionPoint &x1, Vector integrationOrigin, vector<CircularRegion> &circles);
	double integrateHemisphere(list<IntersectionPoint*> &sasa, Vector &integrationOrigin, vector<CircularRegion> &circles);
	double integrateSASA(list<IntersectionPoint*> &sasa, vector<CircularRegion> &circles, Vector &integrationOrigin);
	Vector optimalIntegrationOrigin(list<IntersectionPoint*>* sasa);
	
	
};

#endif //INTEGRATOR_H_
