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

#include "tessellation.h"
#include "depth3d.h"


using namespace std;
using namespace arma;




enum LineType{
	FRONT,
	BACK,
};






typedef multimap<double, LineType> DepthBufferLine;
typedef vector<DepthBufferLine>DepthBuffer;
typedef vector<ScanlineMode>DepthBufferMode;


class MultiLayeredDepthBuffer{
	
public:



	MultiLayeredDepthBuffer(Depth3D &data, int m);
	void addSphere(CircularInterface &circle);
	bool passesBuffer(CircularInterface &circle);
	void print();

	

private:
	int len;
	Depth3D data;
	DepthBuffer dbuffer;
	DepthBufferMode dmode;
	
	Vector normalise(Vector x);
	int sgn(double d);
	Vector tessellationPlaneNormal,  tessellationAxisAuxiliary,  tessellationAxis;
	
	
	ScanlineMode insertIntoLineBuffer(DepthBufferLine &line, double front, double back);
	bool wouldChangeLineBuffer(DepthBufferLine &line, double front, double back);
	DepthBufferLine::iterator increaseLineInterator(DepthBufferLine::iterator it, DepthBufferLine &line);
	DepthBufferLine::iterator decreaseLineInterator(DepthBufferLine::iterator it, DepthBufferLine &line);
	

};

#endif //MULTI_LAYERED_DEPTH_BUFFER_H_
