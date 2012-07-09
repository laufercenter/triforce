#include "topology.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;




Topology::Topology(){
	
}

void Topology::setCell(string ident, Parameters p){
		data[ident]=p;
}

Parameters Topology::getCell(string ident){
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


void Topology::setMassValue(string ident, double v){
		data[ident].mass=v;
}

void Topology::setEpsilonValue(string ident, double v){
		data[ident].epsilon=v;
}

void Topology::setSigmaValue(string ident, double v){
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

