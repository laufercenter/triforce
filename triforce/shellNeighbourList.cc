#include "shellNeighbourList.h"


ShellNeighbourList::ShellNeighbourList(Vector &center, Vector &dim, double searchRadius, int sphericalDetail, int cubicalDetail){
	this->dim = dim;
	this->shellRadius = shellRadius;
	this->center = center;
	this->sphericalDetail = sphericalDetail;
	this->cubicalDetail = cubicalDetail;
	
	
	Coordinate enc;
	Vector c(3);
	Vector s(3);
	Vector v(3);
	double d;
	vector<Coordinate>::iterator it;
	CubicalTemporaryGrid *cubicalTemporaryGrid;
	
	
	numSpheres = VectorInt(3);
	numCubes = VectorInt(3);
	

	//prepare some quantities

	sphericalDistance = searchRadius / sphericalDetail;
	
	
	for(int i=0; i<3; i++)
		numSpheres(i) = dim(i) / sphericalDistance;
	
	shellRadius = searchRadius + (searchRadius / sphericalDetail);
	
	cubicalLength = shellRadius / cubicalDetail;
	
	for(int i=0; i<3; i++)
		numCubes(i) = dim(i) / cubicalLength;
	
	
	//generate a shell-grid
	shellGrid = new ShellGrid(boost::extents[numSpheres(0)][numSpheres(1)][numSpheres(2)]);
	cubicalGrid = new CubicalGrid(boost::extents[numCubes(0)][numCubes(1)][numCubes(2)]);
	cubicalTemporaryGrid = new CubicalTemporaryGrid(boost::extents[numCubes(0)][numCubes(1)][numCubes(2)]);
	
	
	//create associations (slow)
	
	for(int x=0; x<numCubes(0); ++x)
		for(int y=0; y<numCubes(1); ++y)
			for(int z=0; z<numCubes(2); ++z){
				(*cubicalTemporaryGrid)[x][y][z]=numeric_limits<double>::max();
			}
	
	
	for(int x=0; x<numCubes(0); ++x)
		for(int y=0; y<numCubes(1); ++y)
			for(int z=0; z<numCubes(2); ++z){
				c(0) = x*cubicalLength;
				c(1) = y*cubicalLength;
				c(2) = z*cubicalLength;
				
				//this can surely be improved
				for(int x1=0; x1<numSpheres(0); ++x1)
					for(int y1=0; y1<numSpheres(1); ++y1)
						for(int z1=0; z1<numSpheres(2); ++z1){
							s(0) = x*sphericalDistance;
							s(1) = y*sphericalDistance;
							s(2) = z*sphericalDistance;
							
							v = s-c;
							d = norm(v,2);
							if(d<shellRadius){
								enc.x = x1;
								enc.y = y1;
								enc.z = z1;
								
								if(d<(*cubicalTemporaryGrid)[x][y][z]){
									(*cubicalTemporaryGrid)[x][y][z]=d;
									it = ((*cubicalGrid)[x][y][z]).begin();
									++it;
									(*cubicalGrid)[x][y][z].insert(it,enc);
								}
								else (*cubicalGrid)[x][y][z].push_back(enc);
							}
						}
				
			}
				
	origin = center - dim*0.5;

	
}



void ShellNeighbourList::addSphere(Vector &v, int id){
	Vector c;
	int x,y,z;
	Coordinate s;
	c = (v-origin) / cubicalLength;
	x=c(0);
	y=c(1);
	z=c(2);
	
	for(int i=0; i<(*cubicalGrid)[x][y][z].size(); ++i){
		s = (*cubicalGrid)[x][y][z][i];
		(*shellGrid)[s.x][s.y][s.z].insert(id);
	}
	
}

set<int> ShellNeighbourList::getNeighbors(Vector &v){
	Vector c;
	int x,y,z;
	Coordinate s;
	c = (v-origin) / cubicalLength;
	x=c(0);
	y=c(1);
	z=c(2);
	
	s = (*cubicalGrid)[x][y][z][0];
	return (*shellGrid)[s.x][s.y][s.z];
	
}


void ShellNeighbourList::deleteSphere(Vector &v, int id){
	Vector c;
	int x,y,z;
	Coordinate s;
	c = (v-origin) / cubicalLength;
	x=c(0);
	y=c(1);
	z=c(2);
	
	for(int i=0; i<(*cubicalGrid)[x][y][z].size(); ++i){
		s = (*cubicalGrid)[x][y][z][i];
		(*shellGrid)[s.x][s.y][s.z].erase(id);
	}
	
}



