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

#ifndef INTERPOLATIONHCP_H_
#define INTERPOLATIONHCP_H_

#include <string>
#include <vector>
#include <map>
#include "surface3d.h"

#include <armadillo>


using namespace std;
using namespace arma;




class InterpolationHypercubePolytopical: public Interpolator{
	
public:
	InterpolationHypercubePolytopical(Data6D *data);
	float interpolate(Vector &x);

	
private:
	Data6D *data;

	vector<float> weights(vector<VectorInt> &sp, Vector &x, Vector &length);
	
	
};

#endif //INTERPOLATIONHCP_H_
