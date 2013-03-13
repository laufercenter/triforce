#include "molecule.h"

#include <algorithm>
#include <string>





Molecule::Molecule(){
	Molecule(*(new Topology()));
	
}


Molecule::Molecule(Topology topology){
	this->topology=topology;
	hasNeighbourList=false;
	
}

void Molecule::generateNeighbourList(){
	Vector maxPoint(3);
	Vector minPoint(3);
	Vector dim(3);
	Vector center(3);
	double maxRadius;
	Vector v;
	
	center.zeros();
	maxRadius=0;
	
	//calculate grid dimensions and center
	update();
	minPoint=maxPoint=atoms[0];
	for(unsigned int i=0; i<atoms.size(); ++i){
		v = atoms[i];
		for(int j=0; j<3; ++j){
			if(v(j) < minPoint(j)) minPoint(j) = v(j);
			if(v(j) > maxPoint(j)) maxPoint(j) = v(j);
		}
	}
	
	
	dim = (maxPoint - minPoint);
	for(int i=0; i<3; ++i)
		dim(i)=abs(dim(i));
	
	center=minPoint+(dim/2.0);
	

	dim*=1.5;
	
	for(unsigned int i=0; i<radii.size(); ++i){
		maxRadius=max(maxRadius,radii[i]);
	}
	
	
	neighbourList = new NeighbourList(center, dim, maxRadius);
	
	for(unsigned int i=0; i<atoms.size(); ++i){
		neighbourList->addSphere(atoms[i],i);
	}
	
	hasNeighbourList=true;
	
}

vector<int> Molecule::getNeighborListFor(int i){
	if(hasNeighbourList)
		return neighbourList->getNeighbors(atoms[i]);
	else{
		vector<int> s;
		for(unsigned int i=0; i<atoms.size(); ++i)
			s.push_back(i);
		
		return s;
	}
	
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

void Molecule::addInternallyStoredAtom(double x, double y, double z, string name, string type, int i){
	Parameters p;
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	addInternallyStoredAtom(x,y,z,p.sigma,p.epsilon,name,i);
	
}
	
void Molecule::addInternallyStoredAtom(double x, double y, double z, double sigma, double epsilon, string name, int i){
	if(i<0) i = atomicPointers.size();
	
	if(i>=(signed int) atomicPointers.size()){
		constructAtoms(i);
	}
	
	atoms[i](0)=x;
	atoms[i](1)=y;
	atoms[i](2)=z;
	
	
	addAtom(&atoms[i](0),&atoms[i](1),&atoms[i](2),&realArea[realArea.size()-1],&realForceX[realForceX.size()-1],&realForceY[realForceY.size()-1],&realForceZ[realForceZ.size()-1], sigma, epsilon,name,i);
	
	source[i] = INTERNAL_SOURCE;
	
	//pointers might have been invalidated. We need to restore them
	refreshInternalPointers();
}

void Molecule::refreshInternalPointers(){
	for(unsigned int i=0; i < source.size(); ++i){
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


void Molecule::constructAtoms(unsigned int end){
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
		names.resize(end+1,string(""));
	}
}


void Molecule::addAtom(double* x, double* y, double* z, double* area, double* forceX, double* forceY, double* forceZ, string name, string type, int i){
	Parameters p;
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	addAtom(x,y,z,area, forceX, forceY, forceZ, p.sigma,p.epsilon,name,i);
	
	
}


void Molecule::addAtom(double* x, double* y, double* z, double* area, double* forceX, double* forceY, double* forceZ, double sigma, double epsilon, string name, int i){
	AtomicPointers c;
	vector<double*> force;
	c.x=x;
	c.y=y;
	c.z=z;
	
	force.push_back(forceX);
	force.push_back(forceY);
	force.push_back(forceZ);
	
	
	if(i<0) i = atomicPointers.size();
	
	
	if(i>=(signed int) atomicPointers.size()){
		constructAtoms(i);
	}
		
	atomicPointers[i]=c;
	sigmas[i]=sigma;
	epsilons[i]=epsilon;
	radii[i]=0.5*sigma+1.4;
	
	areas[i] = area;
	forces[i] = force;
	source[i] = EXTERNAL_SOURCE;
	names[i]=name;	
}



void Molecule::update(){
	AtomicPointers c;
	for(unsigned i=0; i<atomicPointers.size(); ++i){
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
/*	fprintf(stderr,"Molecule:\n");
	for(int i=0;i<atomicPointers.size();++i){
		fprintf(stderr,"[%d]: (%f, %f, %f), eps: %f, sig: %f, pointers: (%d,%d,%d)\n",i,atoms[i](0),atoms[i](1),atoms[i](2),epsilons[i],sigmas[i],
		       atomicPointers[i].x,atomicPointers[i].y,atomicPointers[i].z);
	}
*/
	//fprintf(stderr,"Areas and forces:\n");
	fprintf(stdout,"index\tname\tarea\tgradx\tgrady\tgradz\tradius\tx\ty\tz\n");
	for(unsigned int i=0;i<areas.size();++i){
		fprintf(stdout,"%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i, names[i].c_str(),*(areas[i]), *(forces[i][0]), *(forces[i][1]), *(forces[i][2]), radii[i], atoms[i](0), atoms[i](1), atoms[i](2));
		//printf("[%d]: pointers: (%d, %d, %d, %d)\n",i,areas[i], forces[i][0], forces[i][1], forces[i][2]);
	}
	
}

void Molecule::printxyz(){
	for(unsigned int i=0;i<areas.size();++i){
		fprintf(stdout,"%.3f %.3f %.3f %.3f\n", atoms[i](0), atoms[i](1), atoms[i](2), radii[i]-1.4);
		//printf("[%d]: pointers: (%d, %d, %d, %d)\n",i,areas[i], forces[i][0], forces[i][1], forces[i][2]);
	}
	
}



void Molecule::printDifference(Molecule *mol){
	//we assume here that the two molecules have same number of atoms etc..
	vector<double*> areas2 = mol->fetchAreaPointers();
	vector<vector<double*> > forces2 = mol->fetchForcePointers();
	
	fprintf(stdout,"index\tname\tarea\tgradx\tgrady\tgradz\tradius\n");
	for(unsigned int i=0;i<areas.size();++i){
		fprintf(stdout,"%d\t%s\t%f\t%f\t%f\t%f\t%f\n",i, names[i].c_str(),*(areas[i])-*(areas2[i]), *(forces[i][0])-*(forces2[i][0]), *(forces[i][1])-*(forces2[i][1]), *(forces[i][2])-*(forces2[i][2]), radii[i]);
	}
	
}

