/**
	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
	Email: nils.drechsel@gmail.com
	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
*/


#include "topology.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;






Topology::Topology(TopologyMode mode){
	this->mode=mode;
	
}

void Topology::setCell(string ident, Parameters p){
		data[ident]=p;
}

Parameters Topology::getCell(string ident){
	if(data.find(ident)==data.end())
		throw AssociationException();

	return data[ident];
}

Parameters Topology::getAssociatedCell(string ident){
	string a = associations[ident];
	return getCell(a);
}


bool Topology::contains(string ident){
	if(data.find(ident)!=data.end()) return true;
	else return false;
}


void Topology::setMassValue(string ident, float v){
		data[ident].mass=v;
}

void Topology::setEpsilonValue(string ident, float v){
		data[ident].epsilon=v;
}

void Topology::setSigmaValue(string ident, float v){
		data[ident].sigma=v;
}

void Topology::setAssociation(string ident0, string ident1){
		associations[ident0]=ident1;
}

string Topology::getAssociation(string ident){
	return associations[ident];
}

void Topology::print(){
	MapVector::iterator it;
	MapString::iterator it2;
	Parameters p;
	printf("DATA:\n");
	for(it=data.begin();it!=data.end();++it){
		p=it->second;
		printf("[%s]: %f %f %f\n",it->first.c_str(),p.mass,p.epsilon,p.sigma);
	}
	printf("ASSOCIATIONS:\n");
	for(it2=associations.begin();it2!=associations.end();++it2){
		printf("[%s]: %s\n",it2->first.c_str(),it2->second.c_str());
	}
	
}

