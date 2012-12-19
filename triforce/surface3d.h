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

#ifndef SURFACE3D_H_
#define SURFACE3D_H_

#include "data3d.h"


class Surface3D: public Data3D{
	
public:
	Surface3D(Data3D* d);
	void surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths);
	
	
private:
	
};

#endif //SURFACE3D_H_
