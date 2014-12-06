/**
	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
	Email: nils.drechsel@gmail.com
	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
*/




#include "benchmark.h"


Benchmark::Benchmark(string section){
	this->section=section;
}


Benchmark::Benchmark(){
	this->section=string("---");
}


clock_t Benchmark::ms(){
	return clock() / (CLOCKS_PER_SEC / 1000);
}

void Benchmark::start(string phase){
	this->phase=phase;
	clock_start=ms();		
}
void Benchmark::stop(){
	clock_t clock_end;
	float t;
	pair<string,float> p;
	clock_end=ms();
	
	t=(float)(clock_end-clock_start);
	if(times.find(phase)==times.end()){
		p.first=string(phase);
		p.second=t;
		times.insert(p);
	}
	else{
		times[phase]=times[phase]+t;
	}
	
}



void Benchmark::print(FILE* outputfile){
	EntityList::iterator it;
	fprintf(outputfile,"%s\n",section.c_str());
	for(it=times.begin(); it!=times.end(); ++it){
		fprintf(outputfile,"\t%s\t\t%.3f ms\n",it->first.c_str(), it->second);
	}
	for(it=stats.begin(); it!=stats.end(); ++it){
		fprintf(outputfile,"\t%s\t\t%.3f\n",it->first.c_str(), it->second);
	}
}


void Benchmark::addQuantity(string quantity, float x){
	pair<string,float> p;
	if(stats.find(quantity)==stats.end()){
		p.first=string(quantity);
		p.second=x;
		stats.insert(p);
	}
	else{
		stats[quantity]=stats[quantity]+x;
	}
	
}


