#include "neighbourList.h"


NeighbourList::NeighbourList(Vector center, Vector dim, float searchRadius){
	this->dim = dim;
	searchRadius*=2;
	this->searchRadius = searchRadius;
	this->center = center;
	
	
	
	numCubes = VectorInt(3);
	

	//prepare some quantities
	
	for(int i=0; i<3; i++)
		numCubes(i) = max(1.0f,ceil(dim(i) / searchRadius));
	
	//printf("creating neighbourlist with dim (%f, %f, %f), center (%f, %f, %f) searchRadius (%f) gridcells (%d %d %d)\n",dim(0),dim(1),dim(2),center(0),center(1),center(2),searchRadius,numCubes(0),numCubes(1),numCubes(2));
	
	
	
	cubicalGrid = new Grid(boost::extents[numCubes(0)][numCubes(1)][numCubes(2)]);
				
	origin = center - dim*0.5;
}



void NeighbourList::addSphere(Vector &v, int id){
	Vector c;
	int x,y,z;
	Coordinate C;
	c = (v-origin) / searchRadius;
	x=c(0);
	y=c(1);
	z=c(2);
	
	//printf("adding sphere at (%d %d %d) (%f %f %f)\n",x,y,z, v(0), v(1), v(2));
	
	
	(*cubicalGrid)[x][y][z].insert(id);
	
	C.x=0;
	C.y=0;
	C.z=0;
	if(spheres.size()<=id){
		spheres.resize(id+1,C);
		dirty.resize(id+1,false);
	}
		
	C.x=x;
	C.y=y;
	C.z=z;
	
	spheres[id]=C;
	
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
				
				if(	x1>=0 && x1<(int)numCubes(0) &&
					y1>=0 && y1<(int)numCubes(1) &&
					z1>=0 && z1<(int)numCubes(2)){
						//printf("retrieving neighbours of cell %d %d %d (%d)\n",x1,y1,z1,(*cubicalGrid)[x1][y1][z1].size());
					
						res.insert(res.end(),(*cubicalGrid)[x1][y1][z1].begin(),(*cubicalGrid)[x1][y1][z1].end());
				}
			}
			
	
	return res;
	
}


void NeighbourList::deleteSphere(int id){
	int x,y,z;
	Coordinate C;
	C=spheres[id];
	x=C.x;
	y=C.y;
	z=C.z;
	(*cubicalGrid)[x][y][z].erase(id);
	
	
	
}



void NeighbourList::update(vector<Vector> &atoms){
	Coordinate C;
	Vector c,v;
	int x,y,z;
	dirty.clear();
	dirty.resize(atoms.size(),false);
	for(unsigned int i=0; i<atoms.size(); ++i){
		v=atoms[i];
		c = (v-origin) / searchRadius;
		x=c(0);
		y=c(1);
		z=c(2);
		
		C=spheres[i];
		if(C.x!=x || C.y!=y || C.z!=z){
			deleteSphere(i);
			addSphere(v,i);
			dirty[i]=true;
		}
		
	}
	
}



bool NeighbourList::isDirty(unsigned int i){
	return dirty[i];
}


























