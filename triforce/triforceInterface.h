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

#ifndef TRIFORCE_INTERFACE_H_
#define TRIFORCE_INTERFACE_H_

#include "datafile.h"
#include "data3d.h"
#include "depth3d.h"
#include "surface3d.h"
#include "interpolation.h"
#include "integrator.h"
#include "molecule.h"

using namespace std;
using namespace arma;



class TriforceInterface{
public:
	TriforceInterface(string path);
	
	double calculateSurfaceArea(Molecule &mol);
	
protected:
	DataFile *df0;
	DataFile *df1;
	DataFile *df2;
	DataFile *df3;
	DataFile *df4;
	
	Data3D *dat0;
	Data3D *dat1;
	Data3D *dat2;
	Data3D *dat3;
	Data3D *dat4;

	Surface3D* surf0;
	Surface3D* surf1;
	Surface3D* surf2;
	Surface3D* surf3;
	Depth3D* depth0;

	Interpolation *interpolator0;
	Interpolation *interpolator1;
	Interpolation *interpolator2;
	Interpolation *interpolator3;
	
	IntegratorTriforce *integrator;	
};




#endif //TRIFORCE_INTERFACE_H_
