#include "topology.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;




Topology::Topology(){
	
}

void Topology::setCell(string ident, vector<double> v){
		data[ident]=v;
}

vector<double> Topology::getCell(string ident){
	return data[ident];
}


bool Topology::contains(string ident){
	if(data.find(ident)!=data.end()) return true;
	else return false;
}


void Topology::setCellValue(string ident, int i, double v){
		data[ident][i]=v;
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
	printf("DATA:\n");
	for(it=data.begin();it!=data.end();++it){
		printf("[%s]:",it->first.c_str());
		for(int i=0; i<it->second.size();++i)
			printf(" %f",it->second[i]);
		printf("\n");
	}
	printf("ASSOCIATIONS:\n");
	for(it2=associations.begin();it2!=associations.end();++it2){
		printf("[%s]: %s\n",it2->first.c_str(),it2->second.c_str());
	}
	
}

