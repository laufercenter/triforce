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

#ifndef INTERPOLATION_POLYTOPICAL_H_
#define INTERPOLATION_POLYTOPICAL_H_

#include <string>
#include <vector>
#include <map>
#include "interpolation.h"

#include <armadillo>


using namespace std;
using namespace arma;





class InterpolationPolytopical: public Interpolator{
	
public:
	InterpolationPolytopical(Data<float> *data, Data<Vector> *weights);
	InterpolationPolytopical(Data<float> *data, Data<Vector> *weights, Interpolator *daisyChain);
	float interpolate(Vector &x);
	

	
private:
	Data<float> *data;
	Data<Vector> *weights;
	unsigned int dim;
	VectorInt dimensions;
	
	
	Vector &fetchAuxiliaryFloat(unsigned int i);
	VectorInt &fetchAuxiliaryInt(unsigned int i);

	vector<VectorInt> getSupportNodes(Vector &d);
	Interpolation *intp;
	VectorInt pGrid,pWeights;
	Vector distsGrid;
	
};

#endif //INTERPOLATION_POLYTOPICAL_H_
