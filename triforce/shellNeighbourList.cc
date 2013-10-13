#include "shellNeighbourList.h"


ShellNeighbourList::ShellNeighbourList(Vector center, Vector dim, float searchRadius, int sphericalDetail, int cubicalDetail){
	this->dim = dim;
	this->shellRadius = shellRadius;
	this->center = center;
	this->sphericalDetail = sphericalDetail;
	this->cubicalDetail = cubicalDetail;
	
	printf("creating neighbourlist with dim (%f, %f, %f), center (%f, %f, %f) and searchRadius (%f)\n",dim(0),dim(1),dim(2),center(0),center(1),center(2),searchRadius);
	
	
	Coordinate enc;
	Vector c(3);
	Vector s(3);
	Vector v(3);
	float d0, d1;
	vector<Coordinate>::iterator it;
	CubicalTemporaryGrid *cubicalTemporaryGrid;
	
	
	numSpheres = VectorInt(3);
	numCubes = VectorInt(3);
	

	//prepare some quantities

	sphericalDistance = searchRadius / sphericalDetail;
	
	
	for(int i=0; i<3; i++)
		numSpheres(i) = dim(i) / sphericalDistance;
	
	shellRadius = searchRadius + (searchRadius / sphericalDetail);
	
	cubicalLength = searchRadius / cubicalDetail;
	
	for(int i=0; i<3; i++)
		numCubes(i) = dim(i) / cubicalLength;
	
	
	printf("auxiliary data: numSpheres (%d, %d, %d), numCubes (%d, %d, %d), shellRadius (%f), sphericalDistance (%f), cubicalLength (%f)\n",numSpheres(0),numSpheres(1),numSpheres(2), numCubes(0), numCubes(1), numCubes(2), shellRadius, sphericalDistance, cubicalLength);
	
	
	//generate a shell-grid
	shellGrid = new ShellGrid(boost::extents[numSpheres(0)][numSpheres(1)][numSpheres(2)]);
	cubicalGrid = new CubicalGrid(boost::extents[numCubes(0)][numCubes(1)][numCubes(2)]);
	cubicalTemporaryGrid = new CubicalTemporaryGrid(boost::extents[numCubes(0)][numCubes(1)][numCubes(2)]);
	
	
	//create associations (slow)
	
	for(int x=0; x<numCubes(0); ++x)
		for(int y=0; y<numCubes(1); ++y)
			for(int z=0; z<numCubes(2); ++z){
				(*cubicalTemporaryGrid)[x][y][z]=numeric_limits<float>::max();
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
							s(0) = x1*sphericalDistance;
							s(1) = y1*sphericalDistance;
							s(2) = z1*sphericalDistance;
							
							v = s-c;
							d0 = norm(v,2);
							d1 = max(norm(v,2)-cubicalLength*0.5, 0.0);
							if(d1<shellRadius){
								enc.x = x1;
								enc.y = y1;
								enc.z = z1;
								
								//printf("X: %d %d %d\n",x,y,z);
								
								if(d0<(*cubicalTemporaryGrid)[x][y][z]){
									(*cubicalTemporaryGrid)[x][y][z]=d0;
									if((*cubicalGrid)[x][y][z].size()==0){
										(*cubicalGrid)[x][y][z].push_back(enc);
									}
									else{
										it = ((*cubicalGrid)[x][y][z]).begin();
										++it;
										(*cubicalGrid)[x][y][z].insert(it,enc);
									}
								}
								else (*cubicalGrid)[x][y][z].push_back(enc);
							}
						}
				
			}
				
	origin = center - dim*0.5;
	delete(cubicalTemporaryGrid);
	
	printf("done\n");
	
/*	
	for(int x=0; x<numCubes(0); ++x)
		for(int y=0; y<numCubes(1); ++y)
			for(int z=0; z<numCubes(2); ++z){
				printf("%d %d %d [%d]:",x,y,z,(*cubicalGrid)[x][y][z].size());
				
				for(int i=0; i<(*cubicalGrid)[x][y][z].size(); ++i)
					printf(" (%d %d %d)",(*cubicalGrid)[x][y][z][i].x,(*cubicalGrid)[x][y][z][i].y,(*cubicalGrid)[x][y][z][i].z);
				printf("\n");
			}

*/	
}



void ShellNeighbourList::addSphere(Vector &v, int id){
	Vector c;
	int x,y,z;
	Coordinate s;
	c = (v-origin) / cubicalLength;
	x=c(0);
	y=c(1);
	z=c(2);
	
	//printf("adding sphere at (%d %d %d)\n",x,y,z);
	
	for(unsigned int i=0; i<(*cubicalGrid)[x][y][z].size(); ++i){
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

	//printf("retrieving neighbours of cell %d %d %d (%d)\n",x,y,z,(*cubicalGrid)[x][y][z].size());
	
	s = (*cubicalGrid)[x][y][z][0];
	//printf("retrieving neighbours in sphere %d %d %d\n",s.x,s.y,s.z);
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
	
	for(unsigned int i=0; i<(*cubicalGrid)[x][y][z].size(); ++i){
		s = (*cubicalGrid)[x][y][z][i];
		(*shellGrid)[s.x][s.y][s.z].erase(id);
	}
	
}



