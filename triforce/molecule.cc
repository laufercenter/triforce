#include "molecule.h"

#include <algorithm>
#include <string>
#include <sstream>


using namespace std;
using namespace arma;





Molecule::Molecule(ForceField forcefield){
	ifstream *ifs;
	string line;
	vector<string> *content;
	Parameters p;
	string atom;
	
	Molecule::forcefield=forcefield;
	
	
	//load parameter database
	switch(forcefield){
		case Amber99SBildn: d = DataFileDigest("Amber99SBildb.csv",Parameters); break;
	}
	
	dict = d.digest();
	
	while(ifs){
		std::getline(*ifs,line);
		content=split(line,' ');
		
		atom=string2UpperCase((*content)[0]);
		
		p.mass=string2double((*content)[1]);
		p.epsilon=string2double((*content)[2]);
		p.sigma=string2double((*content)[3]);
		
		dict[atom]=p;
				
	}
	
	
}








void Molecule::addAtom(int i, double* x, double* y, double* z, string type){
	CoordinatesPointers c;
	Parameters p;
	c.x=x;
	c.y=y;
	c.z=z;
	
	p=dict[string2UpperCase(type)];
	
	
	if(i != coordinatesPointers.size()){
		if(i>coordinatesPointers.size()){
			coordinatesPointers.resize(i+1,c);
			atoms.resize(i+1,vec(3));
			sigmas.resize(i+1,0);
			epsilons.resize(i+1,0);
		}
		
		coordinatesPointers[i]=c;
		sigmas[i]=p.sigma;
		epsilons[i]=p.epsilon;
		
	}
	else{
		coordinatesPointers.push_back(c);
		atoms.push_back(vec(3));
		sigmas.push_back(p.sigma);
		epsilons.push_back(p.epsilon);
	}
		
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

vector<vec> &Molecule::coordinates(){
	return atoms;
}

