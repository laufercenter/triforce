/* Copyright 2013, Nils J. D. Drechsel & Jordi Villà-Freixa
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

#ifndef MULTI_LAYERED_DEPTH_BUFFER_H_
#define MULTI_LAYERED_DEPTH_BUFFER_H_

#include <string>
#include <vector>
#include <set>
#include "boost/multi_array.hpp"

#include <armadillo>

#include "data3d.h"


using namespace std;
using namespace arma;




enum LineType{
	FRONT,
	BACK,
};



typedef multimap<double, LineType> DepthBufferLine;
typedef vector<DepthBufferLine>DepthBuffer;


class MultiLayeredDepthBuffer{
	
public:



	MultiLayeredDepthBuffer(int detail);
	Vector sphericalVector(double phi, double theta);
	Vector spherical2cartesian(Vector s);
	Vector cartesian2spherical(Vector v);
	bool circlef(double pos, Vector &c, double lambda, double &front, double &back);
	void insertIntoLineBuffer(DepthBufferLine &line, double front, double back);
	bool wouldChangeLineBuffer(DepthBufferLine &line, double front, double back);
	void addSphere(Vector &v, double lambda);
	bool passesBuffer(Vector &v, double lambda);

	

private:
	int len;
	DepthBuffer dbuffer;
	
	
};

#endif //MULTI_LAYERED_DEPTH_BUFFER_H_
