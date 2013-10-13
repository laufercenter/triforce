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

#ifndef SHELL_NEIGHBOUR_LIST_H_
#define SHELL_NEIGHBOUR_LIST_H_

#include <string>
#include <vector>
#include <set>
#include "boost/multi_array.hpp"

#include <armadillo>

#include "data3d.h"


using namespace std;
using namespace arma;

typedef struct
{
	int x;
	int y;
	int z;
}
Coordinate;

typedef boost::multi_array<set<int>,3> ShellGrid;
typedef boost::multi_array<vector<Coordinate>,3> CubicalGrid;
typedef boost::multi_array<float,3> CubicalTemporaryGrid;


class ShellNeighbourList{
	
public:
	ShellNeighbourList(Vector center, Vector dim, float searchRadius, int sphericalDetail, int cubicalDetail);
	
	void addSphere(Vector &x, int id);
	set<int> getNeighbors(Vector &x);
	void deleteSphere(Vector &x, int id);
	

private:
	Vector dim;
	Vector center;
	VectorInt numSpheres;
	VectorInt numCubes;
	Vector origin;
	float shellRadius;
	float cubicalLength;
	float sphericalDistance;
	float sphericalDetail;
	float cubicalDetail;
	
	ShellGrid *shellGrid;
	CubicalGrid *cubicalGrid;

	
	
};

#endif //SHELL_NEIGHBOUR_LIST_H_
