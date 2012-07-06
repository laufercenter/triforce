#include "datamapcsv.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;




DataMapCSV::DataMapCSV(){
	
}

void DataMapCSV::setCell(string ident, vector<double> v){
		data[ident]=v;
}



void DataMapCSV::print(){
	MapVector::iterator it;
	
	for(it=data.begin();it!=data.end();++it){
		printf("[%s]:",it->first.c_str());
		for(int i=0; i<it->second.size();++i)
			printf(" %f",it->second[i]);
		printf("\n");
	}
}

