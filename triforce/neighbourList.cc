#include "neighbourList.h"


NeighbourList::NeighbourList(Vector center, Vector dim, double searchRadius){
	this->dim = dim;
	searchRadius*=2;
	this->searchRadius = searchRadius;
	this->center = center;
	
	
	
	numCubes = VectorInt(3);
	

	//prepare some quantities
	
	for(int i=0; i<3; i++)
		numCubes(i) = max(1.0,ceil(dim(i) / searchRadius));
	
	//printf("creating neighbourlist with dim (%f, %f, %f), center (%f, %f, %f) searchRadius (%f) gridcells (%d %d %d)\n",dim(0),dim(1),dim(2),center(0),center(1),center(2),searchRadius,numCubes(0),numCubes(1),numCubes(2));
	
	
	
	cubicalGrid = new Grid(boost::extents[numCubes(0)][numCubes(1)][numCubes(2)]);
				
	origin = center - dim*0.5;
}



void NeighbourList::addSphere(Vector &v, int id){
	Vector c;
	int x,y,z;
	Coordinate s;
	c = (v-origin) / searchRadius;
	x=c(0);
	y=c(1);
	z=c(2);
	
	//printf("adding sphere at (%d %d %d) (%f %f %f)\n",x,y,z, v(0), v(1), v(2));
	
	
	(*cubicalGrid)[x][y][z].insert(id);
	
	
	
}

vector<int> NeighbourList::getNeighbors(Vector &v){
	Vector c;
	int x,y,z;
	int x1,y1,z1;
	vector<int> res;
	
	c = (v-origin) / searchRadius;
	x=c(0);
	y=c(1);
	z=c(2);
	

	
	
	for(int i=-1; i<=1; i++)
		for(int j=-1; j<=1; j++)
			for(int k=-1; k<=1; k++){
				x1 = x+i;
				y1 = y+j;
				z1 = z+k;
			
				//printf("trying to retrieve neighbours of cell %d %d %d\n",x1,y1,z1);
				
				if(	x1>=0 && x1<numCubes(0) &&
					y1>=0 && y1<numCubes(1) &&
					z1>=0 && z1<numCubes(2)){
						//printf("retrieving neighbours of cell %d %d %d (%d)\n",x1,y1,z1,(*cubicalGrid)[x1][y1][z1].size());
					
						res.insert(res.end(),(*cubicalGrid)[x1][y1][z1].begin(),(*cubicalGrid)[x1][y1][z1].end());
				}
			}
			
	
	return res;
	
}


void NeighbourList::deleteSphere(Vector &v, int id){
	Vector c;
	int x,y,z;
	c = (v-origin) / searchRadius;
	x=c(0);
	y=c(1);
	z=c(2);
	(*cubicalGrid)[x][y][z].erase(id);
	
	
	
}



