#include "integratorStatistical.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;


IntegratorStatistical::IntegratorStatistical(){
	IntegratorStatistical(100);
}



IntegratorStatistical::IntegratorStatistical(int trials){
	this->trials=trials;
	
}


float IntegratorStatistical::angle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}



float IntegratorStatistical::integrate(Molecule *molecule, Tessellation *tessellation){
	this->molecule = molecule;
	this->tessellation = tessellation;
	
	SASAsForMolecule sasas;
	float area;
	
	sasas = tessellation->sasas();
	
	area = 0;
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		area += integrateAtomicSASA(sasas[i]);
	}
	
	return area;
	
	
}


float IntegratorStatistical::integrateAtomicSASA(SASAsForAtom sasasForAtom){
	float radius;
	float area;
	Vector v(3);
	int form;
	int sasaCount;
	int sesaCount;
	bool occluded;
	Vector n(3);
	float l;
	
	
	radius = sasasForAtom.radius;
	area = 0;
	sasaCount=0;
	sesaCount=0;
	srand(120);
	
	for(int t= 0; t<trials; ++t){
		v = randu<fvec>(3);
		v(0) = v(0)-0.5;
		v(1) = v(1)-0.5;
		v(2) = v(2)-0.5;
		occluded = false;

		for(int i=0; !occluded && i<sasasForAtom.sasas.size();++i){
			for(int j=0; !occluded && j<sasasForAtom.sasas[i].sasa.size(); ++j){
				n = sasasForAtom.sasas[i].sasa[j].normalForCircularInterface;
				l = sasasForAtom.sasas[i].sasa[j].lambda.rotation;
				form = sasasForAtom.sasas[i].sasa[j].form;
				
				
				//printf("C[%d,%d] (%f,%f,%f) {%f}\n",sasasForAtom.sasas[i].sasa[j].id0,sasasForAtom.sasas[i].sasa[j].id1,sasasForAtom.sasas[i].sasa[j].normalForCircularInterface(0),sasasForAtom.sasas[i].sasa[j].normalForCircularInterface(1),sasasForAtom.sasas[i].sasa[j].normalForCircularInterface(2),sasasForAtom.sasas[i].sasa[j].lambda);
				
				if(form != SPLITTER){
					
					if(form==CONVEX){
						if(angle(v,n) <= l) occluded=true; 
						
						//printf("ANGLE: %f, LAMBDA: %f\n",angle(v,n),l);
					}
					else{
						n = -n;
						if(angle(v,n) >= l) occluded=true; 
					}
				}
			}
			
		}
		
		if(!occluded) sasaCount++;
		else sesaCount++;
	}
	
	printf("SASA: %d, SESA: %d\n",sasaCount,sesaCount);
	
	area = 4*M_PI*(float)sasaCount / ((float)(sasaCount+sesaCount));
	
	return area;
	
}
