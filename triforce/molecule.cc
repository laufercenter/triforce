#include "molecule.h"

#include <algorithm>
#include <string>





Molecule::Molecule(){
	this->topology = new Topology();
	
}


Molecule::Molecule(Topology topology){
	this->topology=topology;

	
}




void Molecule::addRealAtom(double x, double y, double z, string type, int i){
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	addRealAtom(x,y,z,p.sigma,p.epsilon,i);
	
}
	
void Molecule::addRealAtom(double x, double y, double z, double sigma, double epsilon, int i){
	if(i<0) i = atomicPointers.size();
	
	if(i>=atomicPointers.size()){
		constructAtoms(i);
	}
	
	atoms[i](0)=x;
	atoms[i](1)=y;
	atoms[i](2)=z;
	
	realArea.push_back(0);
	realForceX.push_back(0);
	realForceY.push_back(0);
	realForceZ.push_back(0);
	
	addAtom(&atoms[i](0),&atoms[i](1),&atoms[i](2),&realArea[realArea.size()-1],&realForceX[realForceX.size()-1],&realForceY[realForceY.size()-1],&realForceZ[realForceZ.size()-1], sigma, epsilon,i);
}


void Molecule::constructAtoms(int end){
	CoordinatesPointers c;
	if(end>=atomicPointers.size()){
		atomicPointers.resize(end+1,c);
		atoms.resize(end+1,Vector(3));
		sigmas.resize(end+1,0);
		epsilons.resize(end+1,0);
		radii.resize(end+1,0);
		areas.resize(end+1,0);
		forces.resize(end+1,0);
	}
}


void Molecule::addAtom(double* x, double* y, double* z, double* area, double* forceX, double* forceY, double* forceZ, string type, int i){
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	addAtom(x,y,z,p.sigma,p.epsilon,i);	
	
	
}


void Molecule::addAtom(double* x, double* y, double* z, double* area, double* forceX, double* forceY, double* forceZ, double sigma, double epsilon, int i){
	AtomicPointers c;
	Parameters p;
	vector<double*> force;
	c.x=x;
	c.y=y;
	c.z=z;
	
	force.push_back(forceX);
	force.push_back(forceY);
	force.push_back(forceZ);
	
	
	if(i<0) i = atomicPointers.size();
	
	
	if(i>=atomicPointers.size()){
		constructAtoms(i);
	}
		
	atomicPointers[i]=c;
	sigmas[i]=sigma;
	epsilons[i]=epsilon;
	radii[i]=0.5*sigma+1.4;
	
	areas[i] = area;
	forces[i] = force;
		
}



void Molecule::update(){
	int i;
	AtomicPointers c;
	for(i=0; i<atomicPointers.size(); ++i){
		c=atomicPointers[i];
		atoms[i](0)=*c.x;
		atoms[i](1)=*c.y;
		atoms[i](2)=*c.z;
	}
}

vector<Vector> &Molecule::fetchCoordinates(){
	return atoms;
}
vector<double> &Molecule::fetchRadii(){
	return &radii;
}

vector<vector<double*> > &Molecule::fetchForcePointers(){
	return &forces;
}

vector<double*> &Molecule::fetchAreaPointers(){
	return &areas;
}


string Molecule::string2UpperCase(string s){
	string str=s;
	transform(str.begin(), str.end(),str.begin(), ::toupper);
	return str;
}


void Molecule::print(){
	printf("Molecule:\n");
	for(int i=0;i<atomicPointers.size();++i){
		printf("[%d]: (%f, %f, %f), eps: %f, sig: %f, pointers: (%d,%d,%d)\n",i,atoms[i](0),atoms[i](1),atoms[i](2),epsilons[i],sigmas[i],
		       atomicPointers[i].x,atomicPointers[i].y,atomicPointers[i].z);
		
	}
}
