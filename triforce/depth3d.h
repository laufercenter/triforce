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

#ifndef DEPTH3D_H_
#define DEPTH3D_H_

#include "data3d.h"




enum ScanlineMode{
	SCANLINE_FULL,
	SCANLINE_EMPTY,
	SCANLINE_PARTIAL,
	SCANLINE_EXTENDED
};


typedef struct
{
	vector<float> scanline0;
	vector<float> scanline1;
	vector<ScanlineMode> mode;
}
DepthInformation;



class Depth3D: public Data3D<float>{
	
public:
	Depth3D();
	Depth3D(Data3D<float>* d, int slack);
	
	
	DepthInformation &getFloorScanlines(float kappa, float psi, float lambda, bool invert);
	DepthInformation &getCeilScanlines(float kappa, float psi, float lambda, bool invert);
	float getInterpolatedDepth(float g, float kappa, float psi, float lambda, bool flip, int &p0);
	int closestGridPoint(float x);
	void closestGridPoint(Vector &x, VectorInt &p);
	void closestGridPoint(Vector &x, VectorInt &p, Vector &d);
	
	
private:
	int slack;
	Vector lengths;
	DepthInformation depthInfoBuffer;
	DepthInformation &getScanlines(float kappa, float psi, float lambda, bool invert, VectorInt &p);
	
	
};

#endif //DEPTH3D_H_
