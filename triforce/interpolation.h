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

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include <string>
#include <vector>
#include <map>
#include "surface3d.h"

#include <armadillo>


using namespace std;
using namespace arma;


enum TaylorTermination{
	TAYLOR_LINEAR,
	TAYLOR_QUADRATIC,
	TAYLOR_CUBIC
};

class Interpolator{
	virtual float interpolate(Vector &x)=0;
};



class Interpolation: public Interpolator{
	
public:
	Interpolation(Data3D<float> *data, TaylorTermination degree);
	float interpolate(Vector &x);

	
private:
	Data3D<float> *data;
	TaylorTermination degree;

	float taylorExtension(VectorInt &r, Vector &x);
	vector<float> weights(vector<VectorInt> &sp, Vector &x, Vector &length);
	float multiPointTaylor(Vector &x);
	
	
};

#endif //INTERPOLATION_H_
