#include "molecule.h"

#include <algorithm>
#include <string>





Molecule::Molecule(){

	
}


Molecule::Molecule(Topology topology){
	this->topology=topology;

	
}



void Molecule::addRealAtom(double x, double y, double z, string type, int i){
	if(i<0) i = coordinatesPointers.size();
	
	if(i>=coordinatesPointers.size()){
		constructAtoms(i);
	}
	
	atoms[i](0)=x;
	atoms[i](1)=y;
	atoms[i](2)=z;
	
	addAtom(&atoms[i](0),&atoms[i](1),&atoms[i](2),type,i);
}


void Molecule::constructAtoms(int end){
	CoordinatesPointers c;
	if(end>=coordinatesPointers.size()){
		coordinatesPointers.resize(end+1,c);
		atoms.resize(end+1,Vector(3));
		sigmas.resize(end+1,0);
		epsilons.resize(end+1,0);
		radii.resize(end+1,0);
	}
}


void Molecule::addAtom(double* x, double* y, double* z, string type, int i){
	CoordinatesPointers c;
	Parameters p;
	c.x=x;
	c.y=y;
	c.z=z;
	
	if(i<0) i = coordinatesPointers.size();
	
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	if(i>=coordinatesPointers.size()){
		constructAtoms(i);
	}
		
	coordinatesPointers[i]=c;
	sigmas[i]=p.sigma;
	epsilons[i]=p.epsilon;
	radii[i]=p.sigma+1.4;
		
		
}



void Molecule::update(){
	int i;
	CoordinatesPointers c;
	for(i=0; i<coordinatesPointers.size(); ++i){
		c=coordinatesPointers[i];
		atoms[i](0)=*c.x;
		atoms[i](1)=*c.y;
		atoms[i](2)=*c.z;
	}
}

vector<Vector> &Molecule::coordinates(){
	return atoms;
}
vector<double>* Molecule::fetchRadii(){
	return &radii;
}



string Molecule::string2UpperCase(string s){
	string str=s;
	transform(str.begin(), str.end(),str.begin(), ::toupper);
	return str;
}


void Molecule::print(){
	printf("Molecule:\n");
	for(int i=0;i<coordinatesPointers.size();++i){
		printf("[%d]: (%f, %f, %f), eps: %f, sig: %f, pointers: (%d,%d,%d)\n",i,atoms[i](0),atoms[i](1),atoms[i](2),epsilons[i],sigmas[i],
		       coordinatesPointers[i].x,coordinatesPointers[i].y,coordinatesPointers[i].z);
		
	}
}
