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


struct CoordinateComparator: public std::binary_function<Coordinate, Coordinate, bool>
{
	bool operator()(const Coordinate& lhs, const Coordinate& rhs) const
	{
		if(lhs.x == rhs.x){
			if(lhs.y==rhs.y){
				return lhs.z < rhs.z;
					
			}
			else{
				return lhs.y < rhs.y;
			}
		}
		else return lhs.x < rhs.x;
	}
};

typedef struct
{
	float eps;
	float sig;
}
LJPair;

typedef boost::multi_array<set<int>,3> Grid;
typedef boost::multi_array<LJPair,3> LJGrid;
typedef boost::multi_array<bool,3> DirtyGrid;


class NeighbourList{
	
public:
	NeighbourList(Vector center, Vector dim, float searchRadius, vector<float> *epsilons, vector<float> *sigmas);
	void addSphere(Vector &x, unsigned int id);
	vector<int> getNeighbors(Vector &x);
	void deleteSphere(unsigned int id);
	void update(vector<Vector> &atoms);
	bool isDirty(unsigned int i);
	map<Coordinate, Coordinate, CoordinateComparator>& getCellList();
	float lennardJonesImpact(Vector &x, Coordinate &distantCell);
	void print(FILE* outputfile0,FILE* outputfile1);

private:
	Vector dim;
	Vector center;
	VectorInt numCubes;
	Vector origin;
	float searchRadius;
	vector<Coordinate> spheres;
	vector<bool> dirty;
	
	Grid *cubicalGrid;
	LJGrid *lennardJonesGrid;
	vector<float> *epsilons;
	vector<float> *sigmas;
	map<Coordinate, Coordinate, CoordinateComparator> cellList;
	Vector cellPosition(Coordinate &x);
	float lennardJonesPotential(LJPair &ljPair, float x);
	
	
};

#endif //NEIGHBOUR_LIST_H_
