#include "depth3d.h"

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;


#define doublepi 6.283185307179586231996
#define halfpi 1.570796326794896557999


Depth3D::Depth3D(){
}

Depth3D::Depth3D(Data3D* d, int slack){
	
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
	
	this->slack=slack;

	depthInfoBuffer.scanline0.resize(parameter0Dim,-1);
	depthInfoBuffer.scanline1.resize(parameter0Dim,-1);
	depthInfoBuffer.mode.resize(parameter0Dim,SCANLINE_EMPTY);

}



int Depth3D::closestGridPoint(double x){
	unsigned int i_parameter0;
	
	i_parameter0 = static_cast<unsigned int>(floor((x-minParameter0)/cellLengthParameter0));

		return i_parameter0;
	
}


void Depth3D::closestGridPoint(Vector &x, VectorInt &p){
	unsigned int i_parameter1;
	unsigned int i_parameter2;
	
	
	i_parameter1 = static_cast<unsigned int>(floor((x(0)-minParameter1)/cellLengthParameter1));
	
	i_parameter2 = static_cast<unsigned int>(floor((x(1)-minParameter2)/cellLengthParameter2));
	
	
	p(0) = i_parameter1;
	p(1) = i_parameter2;
}



void Depth3D::closestGridPoint(Vector &x, VectorInt &p, Vector &l){
	unsigned int i_parameter0;
	unsigned int i_parameter1;
	unsigned int i_parameter2;
	
	l(0) = cellLengthParameter0;
	l(1) = cellLengthParameter1;
	l(2) = cellLengthParameter2;
	
	
	i_parameter0 = static_cast<unsigned int>(floor((x(0)-minParameter0)/cellLengthParameter0));
	
	i_parameter1 = static_cast<unsigned int>(floor((x(1)-minParameter1)/cellLengthParameter1));
	
	i_parameter2 = static_cast<unsigned int>(floor((x(2)-minParameter2)/cellLengthParameter2));
	
	
	p(0) = i_parameter0;
	p(1) = i_parameter1;
	p(2) = i_parameter2;
}




double Depth3D::getInterpolatedDepth(double g, double kappa, double psi, double lambda, bool flip, int &p0){
	VectorInt p(3);
	Vector x(3);
	Vector l(3);
	double a,b,c;
	double base;
	double offset;
	double sigma;
	
	x(0)=g;
	x(1)=psi;
	x(2)=lambda;
	closestGridPoint(x,p,l);
	p0=p(0);
	
	printf("P: %d %d %d\n",p(0),p(1),p(2));
	
	if(p(0)==parameter0Dim-1){
		a = (*data)[p(0)][p(1)][p(2)];
		a = min(max(a,0.0), M_PI);
		
		if(flip) a = (doublepi)-a;
		return a + kappa;
	}
	else{
		//perform linear interpolation
		base=l(0)*p(0);
		offset=g-base;
		sigma=offset/l(0);
		a = (*data)[p(0)][p(1)][p(2)];
		a = min(max(a,0.0), M_PI);

		b = (*data)[p(0)+1][p(1)][p(2)];
		b = min(max(b,0.0), M_PI);
		
		c = a*(1.0-sigma) + b*sigma;
		if(flip) c = (doublepi)-c;

		return  c + kappa;
		
	}
		
	
	
	
}


DepthInformation &Depth3D::getFloorScanlines(double kappa, double psi, double lambda, bool invert){
	VectorInt p(3);
	Vector x(2);
	x(0)=psi;
	x(1)=lambda;
	closestGridPoint(x,p);

	if(!invert)
		p(1)=max(0, (int)p(1)-slack);
		//p(1)=max(0, (int)p(1));
	
	if(invert)
		p(1)=min((int)parameter2Dim-1, (int)p(1)+slack);
		//p(1)=min((int)parameter2Dim-1, (int)p(1));
	
	return getScanlines(kappa, psi, lambda, invert, p);
}



DepthInformation &Depth3D::getCeilScanlines(double kappa, double psi, double lambda, bool invert){
	VectorInt p(3);
	Vector x(2);
	x(0)=psi;
	x(1)=lambda;
	
	closestGridPoint(x,p);
	
	return getScanlines(kappa, psi, lambda, invert, p);
}
	
	
	
DepthInformation &Depth3D::getScanlines(double kappa, double psi, double lambda, bool invert, VectorInt &p){
	Vector x(2);
	double k;
	Vector l(3); 
	
	x(0)=psi;
	x(1)=lambda;
	
	
	//scanline0 represents the beginning of the occluded region, scanline1 the end
	if(!invert){
		for(unsigned int i=0; i<parameter0Dim; ++i){
			k = (*data)[i][p(0)][p(1)];
			if(k<0) depthInfoBuffer.mode[i] = SCANLINE_FULL;
			else if(k>doublepi) depthInfoBuffer.mode[i] = SCANLINE_EMPTY; 
			else{
				depthInfoBuffer.mode[i] = SCANLINE_PARTIAL;
				depthInfoBuffer.scanline1[i] = k+kappa;
				if(depthInfoBuffer.scanline1[i]>doublepi) depthInfoBuffer.scanline1[i]-=doublepi;
				if(depthInfoBuffer.scanline1[i]>doublepi) depthInfoBuffer.scanline1[i]-=doublepi;
				
				if(depthInfoBuffer.scanline1[i]<0) depthInfoBuffer.scanline1[i]+=doublepi;
				if(depthInfoBuffer.scanline1[i]<0) depthInfoBuffer.scanline1[i]+=doublepi;
				
				
				depthInfoBuffer.scanline0[i] = (doublepi-k)+kappa;
				if(depthInfoBuffer.scanline0[i]>doublepi) depthInfoBuffer.scanline0[i]-=doublepi;
				if(depthInfoBuffer.scanline0[i]>doublepi) depthInfoBuffer.scanline0[i]-=doublepi;
				
				if(depthInfoBuffer.scanline0[i]<0) depthInfoBuffer.scanline0[i]+=doublepi;
				if(depthInfoBuffer.scanline0[i]<0) depthInfoBuffer.scanline0[i]+=doublepi;
				
			}
		}
	}
	else{
		for(unsigned int i=0; i<parameter0Dim; ++i){
			k = (*data)[i][p(0)][p(1)];
			if(k<0) depthInfoBuffer.mode[i] = SCANLINE_EMPTY;
			else if(k>doublepi) depthInfoBuffer.mode[i] = SCANLINE_FULL; 
			else{
				depthInfoBuffer.mode[i] = SCANLINE_PARTIAL;
				depthInfoBuffer.scanline0[i] = k+kappa;
				if(depthInfoBuffer.scanline0[i]>doublepi) depthInfoBuffer.scanline0[i]-=doublepi;
				if(depthInfoBuffer.scanline0[i]>doublepi) depthInfoBuffer.scanline0[i]-=doublepi;
				
				depthInfoBuffer.scanline1[i] = (doublepi-k)+kappa;
				if(depthInfoBuffer.scanline1[i]>doublepi) depthInfoBuffer.scanline1[i]-=doublepi;
				if(depthInfoBuffer.scanline1[i]>doublepi) depthInfoBuffer.scanline1[i]-=doublepi;
			}
		}
	}
	
	return depthInfoBuffer;
}



