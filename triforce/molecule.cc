#include "molecule.h"

#include <algorithm>
#include <string>





Molecule::Molecule(){
	this->topology = *(new Topology());
	
}


Molecule::Molecule(Topology topology){
	this->topology=topology;

	
}


Vector Molecule::getInternallyStoredAtomCoordinates(int i){
	return atoms[i];
}

void Molecule::setInternallyStoredAtomCoordinates(int i, Vector &v){
	atoms[i] = v;
}

void Molecule::perturbInternallyStoredAtomCoordinates(int i, Vector p){
	atoms[i] = atoms[i] + p;
}

void Molecule::addInternallyStoredAtom(double x, double y, double z, string type, int i){
	Parameters p;
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	addInternallyStoredAtom(x,y,z,p.sigma,p.epsilon,i);
	
}
	
void Molecule::addInternallyStoredAtom(double x, double y, double z, double sigma, double epsilon, int i){
	if(i<0) i = atomicPointers.size();
	
	if(i>=atomicPointers.size()){
		constructAtoms(i);
	}
	
	atoms[i](0)=x;
	atoms[i](1)=y;
	atoms[i](2)=z;
	
	
	addAtom(&atoms[i](0),&atoms[i](1),&atoms[i](2),&realArea[realArea.size()-1],&realForceX[realForceX.size()-1],&realForceY[realForceY.size()-1],&realForceZ[realForceZ.size()-1], sigma, epsilon,i);
	
	source[i] = INTERNAL_SOURCE;
	
	//pointers might have been invalidated. We need to restore them
	refreshInternalPointers();
}

void Molecule::refreshInternalPointers(){
	int i;
	for(i=0; i < source.size(); ++i){
		if(source[i]==INTERNAL_SOURCE){
			atomicPointers[i].x = &(atoms[i](0));
			atomicPointers[i].y = &(atoms[i](1));
			atomicPointers[i].z = &(atoms[i](2));
			areas[i] = &(realArea[i]);
			forces[i][0] = &(realForceX[i]);
			forces[i][1] = &(realForceY[i]);
			forces[i][2] = &(realForceZ[i]);
		}
	}
}


void Molecule::constructAtoms(int end){
	AtomicPointers c;
	if(end>=atomicPointers.size()){
		atomicPointers.resize(end+1,c);
		atoms.resize(end+1,Vector(3));
		sigmas.resize(end+1,0);
		epsilons.resize(end+1,0);
		radii.resize(end+1,0);
		areas.resize(end+1,0);
		forces.resize(end+1,vector<double*>());
		realArea.resize(end+1,0);
		realForceX.resize(end+1,0);
		realForceY.resize(end+1,0);
		realForceZ.resize(end+1,0);
		source.resize(end+1,EXTERNAL_SOURCE);
	}
}


void Molecule::addAtom(double* x, double* y, double* z, double* area, double* forceX, double* forceY, double* forceZ, string type, int i){
	Parameters p;
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	addAtom(x,y,z,area, forceX, forceY, forceZ, p.sigma,p.epsilon,i);
	
	
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
	source[i] = EXTERNAL_SOURCE;
		
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
	return radii;
}

vector<vector<double*> > &Molecule::fetchForcePointers(){
	return forces;
}

vector<double*> &Molecule::fetchAreaPointers(){
	return areas;
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
	printf("Areas and forces:\n");
	for(int i=0;i<areas.size();++i){
		printf("[%d]: %f (%f, %f, %f) pointers: (%d, %d, %d)\n",i,*(areas[i]), *(forces[i][0]), *(forces[i][1]), *(forces[i][2]), forces[i][0], forces[i][1], forces[i][2]);
		//printf("[%d]: pointers: (%d, %d, %d, %d)\n",i,areas[i], forces[i][0], forces[i][1], forces[i][2]);
	}
	
}
