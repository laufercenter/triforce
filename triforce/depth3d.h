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
	vector<double> scanline0;
	vector<double> scanline1;
	vector<ScanlineMode> mode;
}
DepthInformation;



class Depth3D: public Data3D{
	
public:
	Depth3D();
	Depth3D(Data3D* d, int slack);
	
	
	DepthInformation &getFloorScanlines(double kappa, double psi, double lambda, bool invert);
	DepthInformation &getCeilScanlines(double kappa, double psi, double lambda, bool invert);
	double getInterpolatedDepth(double g, double kappa, double psi, double lambda, bool flip, int &p0);
	int closestGridPoint(double x);
	void closestGridPoint(Vector &x, VectorInt &p);
	void closestGridPoint(Vector &x, VectorInt &p, Vector &l);
	
	
private:
	int slack;
	DepthInformation depthInfoBuffer;
	DepthInformation &getScanlines(double kappa, double psi, double lambda, bool invert, VectorInt &p);
	
	
};

#endif //DEPTH3D_H_
