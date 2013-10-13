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
	float maxRadius;
	Vector v;
	
	if(atoms.size()==0) return;
	
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


unsigned int Molecule::identifySpecies(float eps, float sig){
	float realisticEps[6] = {3.39967, 2.64953,  3.25000, 2.95992, 3.74177, 3.56359};
	float realisticSig[6] = {0.35982400, 0.06568880, 0.71128000, 0.87864000, 0.83680000, 1.04600000};
	float bestScore = std::numeric_limits<float>::max();
	unsigned int bestSpecies=0;
	float score;
	//go through all the species, select the one that fits best
	for(unsigned int i=0; i<N_SPECIES; ++i){
		score = (eps-realisticEps[i]) * (eps-realisticEps[i]) + (sig-realisticSig[i]) * (eps-realisticSig[i]);
		if(score<bestScore){
			bestSpecies=i;
			bestScore=score;
		}
	}
	return bestSpecies;

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

void Molecule::addInternallyStoredAtom(float x, float y, float z, string name, string type, int i){
	Parameters p;
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	addInternallyStoredAtom(x,y,z, 0.5*p.sigma+1.4, name,i);
	
}

void Molecule::addInternallyStoredAtom(float x, float y, float z, float radius, string name, int i){
	addInternallyStoredAtom(x,y,z, radius, 0, 0, name,i);
}

void Molecule::addInternallyStoredAtom(float x, float y, float z, float eps, float sig, string name, int i){
	addInternallyStoredAtom(x,y,z, sig, eps, sig, name,i);
}
	
	
void Molecule::addInternallyStoredAtom(double x, double y, double z, double radius, string name, int i){
	if(i<0) i = atomicPointers.size();
	
	if(i>=(signed int) atomicPointers.size()){
		constructAtoms(i);
	}
	
	atoms[i](0)=x;
	atoms[i](1)=y;
	atoms[i](2)=z;
	
	
	addAtom(&atoms[i](0),&atoms[i](1),&atoms[i](2),&realArea[realArea.size()-1],&realForceX[realForceX.size()-1],&realForceY[realForceY.size()-1],&realForceZ[realForceZ.size()-1], radius, name,i);
	
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
	c.x=0;
	c.y=0;
	c.z=0;
	if(end>=atomicPointers.size()){
		atomicPointers.resize(end+1,c);
		atoms.resize(end+1,Vector(3));
		epsilons.resize(end+1,0);
		sigmas.resize(end+1,0);
		radii.resize(end+1,0);
		species.resize(end+1,0);
		areas.resize(end+1,0);
		forces.resize(end+1,vector<float*>());
		realArea.resize(end+1,0);
		realForceX.resize(end+1,0);
		realForceY.resize(end+1,0);
		realForceZ.resize(end+1,0);
		source.resize(end+1,EXTERNAL_SOURCE);
		names.resize(end+1,string(""));
	}
}


void Molecule::addAtom(float* x, float* y, float* z, float* area, float* forceX, float* forceY, float* forceZ, string name, string type, int i){
	Parameters p;
	p=topology.getAssociatedCell(string2UpperCase(type));
	
	addAtom(x,y,z,area, forceX, forceY, forceZ, 0.5*p.sigma+1.4, name,i);
	
	
}


	
	
void Molecule::addAtom(double* x, double* y, double* z, double* area, double* forceX, double* forceY, double* forceZ, double radius, string name, int i){
	AtomicPointers c;
	vector<float*> force;
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
	radii[i]=radius;
	
	//sea water modifications
	epsilons[i]=sqrt(eps);
	sigmas[i]=(sig+0.82)*0.5;
	
	areas[i] = area;
	forces[i] = force;
	source[i] = EXTERNAL_SOURCE;
	names[i]=name;
	
	//determine species
	species[i] = identifySpecies(eps,sig);
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
vector<float> &Molecule::fetchRadii(){
	return radii;
}

vector<float> &Molecule::fetchEpsilons(){
	return epsilons;
}

vector<float> &Molecule::fetchSigmas(){
	return sigmas;
}

vector<unsigned int> &Molecule::fetchSpecies(){
	return species;
}

vector<vector<float*> > &Molecule::fetchForcePointers(){
	return forces;
}

vector<float*> &Molecule::fetchAreaPointers(){
	return areas;
}

vector<float*> &Molecule::fetchNonpolarPointers(){
	return nonpolarFreeEnergies;
}

string Molecule::string2UpperCase(string s){
	string str=s;
	transform(str.begin(), str.end(),str.begin(), ::toupper);
	return str;
}


void Molecule::print(FILE* outputfile){
/*	fprintf(stderr,"Molecule:\n");
	for(int i=0;i<atomicPointers.size();++i){
		fprintf(stderr,"[%d]: (%f, %f, %f), eps: %f, sig: %f, pointers: (%d,%d,%d)\n",i,atoms[i](0),atoms[i](1),atoms[i](2),epsilons[i],sigmas[i],
		       atomicPointers[i].x,atomicPointers[i].y,atomicPointers[i].z);
	}
*/
	//fprintf(stderr,"Areas and forces:\n");
	fprintf(outputfile,"index\tname\tarea\tgradx\tgrady\tgradz\tradius\tx\ty\tz\n");
	for(unsigned int i=0;i<areas.size();++i){
		fprintf(outputfile,"%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i, names[i].c_str(),*(areas[i]), *(forces[i][0]), *(forces[i][1]), *(forces[i][2]), radii[i], atoms[i](0), atoms[i](1), atoms[i](2));
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
	vector<float*> areas2 = mol->fetchAreaPointers();
	vector<vector<float*> > forces2 = mol->fetchForcePointers();
	
	fprintf(stdout,"index\tname\tarea\tgradx\tgrady\tgradz\tradius\n");
	for(unsigned int i=0;i<areas.size();++i){
		fprintf(stdout,"%d\t%s\t%f\t%f\t%f\t%f\t%f\n",i, names[i].c_str(),*(areas[i])-*(areas2[i]), *(forces[i][0])-*(forces2[i][0]), *(forces[i][1])-*(forces2[i][1]), *(forces[i][2])-*(forces2[i][2]), radii[i]);
	}
	
}

