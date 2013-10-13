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

#ifndef NEIGHBOUR_LIST_H_
#define NEIGHBOUR_LIST_H_

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

typedef boost::multi_array<set<int>,3> Grid;


class NeighbourList{
	
public:
	NeighbourList(Vector center, Vector dim, float searchRadius);
	
	void addSphere(Vector &x, int id);
	vector<int> getNeighbors(Vector &x);
	void deleteSphere(Vector &x, int id);
	

private:
	Vector dim;
	Vector center;
	VectorInt numCubes;
	Vector origin;
	float searchRadius;
	
	Grid *cubicalGrid;

	
	
};

#endif //NEIGHBOUR_LIST_H_