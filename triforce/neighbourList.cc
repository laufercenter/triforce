#include "neighbourList.h"


NeighbourList::NeighbourList(Vector center, Vector dim, float searchRadius, vector<float> *epsilons, vector<float> *sigmas){
	LJPair ljPair;
	this->dim = dim;
	searchRadius*=2; //get the diameter
	this->searchRadius = searchRadius;
	this->center = center;
	
	
	
	numCubes = VectorInt(3);
	

	//prepare some quantities
	
	for(int i=0; i<3; i++)
		numCubes(i) = max(1.0f,ceil(dim(i) / searchRadius));
	
	printf("creating neighbourlist with dim (%f, %f, %f), center (%f, %f, %f) searchRadius (%f) gridcells (%d %d %d)\n",dim(0),dim(1),dim(2),center(0),center(1),center(2),searchRadius,numCubes(0),numCubes(1),numCubes(2));
	
	
	
	cubicalGrid = new Grid(boost::extents[numCubes(0)][numCubes(1)][numCubes(2)]);
	lennardJonesGrid = new LJGrid(boost::extents[numCubes(0)][numCubes(1)][numCubes(2)]);
	ljPair.eps=0;
	ljPair.sig=0;
	std::fill(lennardJonesGrid->data(), lennardJonesGrid->data()+lennardJonesGrid->num_elements(),ljPair);
				
	origin = center - dim*0.5;
}




void NeighbourList::addSphere(Vector &v, unsigned int id){
	Vector c;
	int x,y,z;
	Coordinate C;
	c = (v-origin) / searchRadius;
	x=c(0);
	y=c(1);
	z=c(2);
	
	//printf("adding sphere at (%d %d %d) (%f %f %f)\n",x,y,z, v(0), v(1), v(2));
	
	
	(*cubicalGrid)[x][y][z].insert(id);
	(*lennardJonesGrid)[x][y][z].eps+=(*epsilons)[id];
	(*lennardJonesGrid)[x][y][z].sig+=(*sigmas)[id];
	
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
	
	if((*cubicalGrid)[x][y][z].size()==1) cellList[C]=C;
	
	
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



void NeighbourList::deleteSphere(unsigned int id){
	int x,y,z;
	Coordinate C;
	C=spheres[id];
	x=C.x;
	y=C.y;
	z=C.z;
	(*cubicalGrid)[x][y][z].erase(id);
	(*lennardJonesGrid)[x][y][z].eps-=(*epsilons)[id];
	(*lennardJonesGrid)[x][y][z].sig-=(*sigmas)[id];
	
	if((*cubicalGrid)[x][y][z].size()==0){
		cellList.erase(C);
		(*lennardJonesGrid)[x][y][z].eps=0;
		(*lennardJonesGrid)[x][y][z].sig=0;
	}
	
	
	
	
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




map<Coordinate, Coordinate, CoordinateComparator>& NeighbourList::getCellList(){
	return cellList;
}

Vector NeighbourList::cellPosition(Coordinate &x){
	Vector c,r;
	c(0)=x.x;
	c(1)=x.y;
	c(2)=x.z;
	c = c*searchRadius+origin;
	r(0)=searchRadius*0.5;
	r(1)=searchRadius*0.5;
	r(2)=searchRadius*0.5;
	c+=r; //this will move it into the center of the cell
	return c;
}

float NeighbourList::lennardJonesPotential(LJPair &ljPair, float x){
	float sig6;
	float sig12;
	float x6;
	float x12;
	
	sig6=pow(ljPair.sig,6);
	sig12=sig6*sig6;
	
	x6=pow(x,6);
	x12=x*x;
	
	
	return 4*ljPair.eps*(sig12/x12 - sig6/x6);
}


float NeighbourList::lennardJonesImpact(Vector &x, Coordinate &distantCell){
	Vector dc;
	float dist;
	LJPair ljPair;
	float n;
	
	ljPair = (*lennardJonesGrid)[distantCell.x][distantCell.y][distantCell.z];
	n = (*cubicalGrid)[distantCell.x][distantCell.y][distantCell.z].size();
	ljPair.eps/=n;
	ljPair.sig/=n;
	
	dc=cellPosition(distantCell);
	dist=norm(x-dc,2);
	return lennardJonesPotential(ljPair,dist)*n;
}

void NeighbourList::print(FILE* outputfile0,FILE* outputfile1){
	Vector x;
	LJPair ljPair;
	unsigned int n;
	fprintf(outputfile0,"X\tY\tZ\tx\ty\tz\tn\teps\tsig\n");
	map<Coordinate, Coordinate, CoordinateComparator>::iterator it;
	for(it=cellList.begin(); it!=cellList.end(); ++it){
		x=cellPosition(it->second);
		ljPair=(*lennardJonesGrid)[it->first.x][it->first.y][it->first.z];
		n = (*cubicalGrid)[it->first.x][it->first.y][it->first.z].size();
		fprintf(outputfile0,"%d\t%d\t%d\t%f\t%f\t%f\t%d\t%f\t%f\n",it->first.x, it->first.y, it->first.y, x(0), x(1), x(2), n, ljPair.eps, ljPair.sig);
	}
	
	fprintf(outputfile1,"index\tX\tY\tZ\n");
	for(unsigned int i=0; i<spheres.size(); ++i){
		fprintf(outputfile1,"%d\t%d\t%d\t%d\n",i, spheres[i].x, spheres[i].y, spheres[i].z);
	}
	
}

















