#include "depth3d.h"

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;


Depth3D::Depth3D(Data3D* dforward, Data3D* dbackward){
	
	parameter0Dim=d->parameter0Dim;
	parameter1Dim=d->parameter1Dim;
	parameter2Dim=d->parameter2Dim;
	minParameter0=d->minParameter0;
	minParameter1=d->minParameter1;
	minParameter2=d->minParameter2;
	cellLengthParameter0=d->cellLengthParameter0;
	cellLengthParameter1=d->cellLengthParameter1;
	cellLengthParameter2=d->cellLengthParameter2;
	headerParameter0=d->headerParameter0;
	headerParameter1=d->headerParameter1;
	headerParameter2=d->headerParameter2;
	data=d->data;


}




void Data3D::closestGridPoint(Vector &x, VectorInt &p){
	int i_parameter1;
	int i_parameter2;
	
	
	i_parameter1 = static_cast<int>(floor((x(0)-minParameter1)/cellLengthParameter1));
	
	i_parameter2 = static_cast<int>(floor((x(1)-minParameter2)/cellLengthParameter2));
	
	
	p(0) = i_parameter1;
	p(1) = i_parameter2;
}





DepthInformation Depth3D::getScanlines(double kappa, double psi, double lambda, CircularInterfaceForm form){
	VectorInt p(3);
	Vector x(2);
	int offset;
	int j;
	double k;
	DepthInformation dat;
	
	x(0)=psi;
	x(1)=lambda;
	dat.scanline0.resize(parameter0Dim,-1);
	dat.scanline1.resize(parameter0Dim,-1);
	dat.mode.resize(parameter0Dim,-1);
	
	closestGridPoint(x,p);
	
	
	//scanline0 represents the beginning of the occluded region, scanline1 the end
	if(form==CONVEX){
		for(int i=0; i<parameter0Dim; ++i){
			k = data[i][p(0)][p(1)];
			if(k<0) dat.mode[i] = OCCLUDED;
			else if(k>2*M_PI) dat.mode[i] = EMPTY; 
			else{
				dat.scanline1 = k+kappa;
				if(dat.scanline1>2*M_PI) dat.scanline1-=s*M_PI;
				
				dat.scanline0 = (2*M_PI-k)+kappa;
				if(dat.scanline0>2*M_PI) dat.scanline0-=s*M_PI;
			}
		}
	}
	else{
		for(int i=0; i<parameter0Dim; ++i){
			k = data[i][p(0)][p(1)];
			if(k<0) dat.mode[i] = EMPTY;
			else if(k>2*M_PI) dat.mode[i] = OCCLUDED; 
			else{
				dat.scanline0 = k+kappa;
				if(dat.scanline0>2*M_PI) dat.scanline0-=s*M_PI;
				
				dat.scanline1 = (2*M_PI-k)+kappa;
				if(dat.scanline1>2*M_PI) dat.scanline1-=s*M_PI;
			}
		}
	}
	
	return dat;
}



