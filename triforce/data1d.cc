#include "data1d.h"

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;

#define THRESHOLD_NUMERICAL 0.00001

Data1D::Data1D(){
}

Data1D::Data1D(unsigned int parameter0Dim, bool containsAuxiliaryData){
	this->parameter0Dim = parameter0Dim;
	
		
	this->containsAuxiliaryData = containsAuxiliaryData;
		
	headerParameter0 = new Table1dDouble(boost::extents[parameter0Dim]);
	data = new Table1dDouble(boost::extents[parameter0Dim]);
	if(containsAuxiliaryData){
		auxiliary = new Table1dDouble(boost::extents[parameter0Dim]);
	}
	
}


void Data1D::init(){
	cellLengthParameter0 = abs((*headerParameter0)[1]-(*headerParameter0)[0]);
	
	minParameter0 = (*headerParameter0)[0];
}


void Data1D::setHeaderParameter0Cell(unsigned int x, float value){
	(*headerParameter0)[x] = value;
}

void Data1D::setAuxiliaryCell(unsigned int x, float value){
	(*auxiliary)[x] = value;
}

float Data1D::getHeaderParameter0Cell(unsigned int x){
	return (*headerParameter0)[x];
}


float Data1D::getAuxiliaryCell(unsigned int x){
	return (*auxiliary)[x];
}

void Data1D::setDataCell(unsigned int x, float value){
	(*data)[x] = value;
}


Vector Data1D::getHeaderVector(unsigned int parameter0){
	Vector v=Vector(1);
	v(0) = (*headerParameter0)[parameter0];
	
	return v;
}




float Data1D::getDataCell(unsigned int x){
	return (*data)[x];
}


bool Data1D::isWithinNumericalLimits(float x, float t){
	if(abs(x-t) <= THRESHOLD_NUMERICAL) return true;
	else return false;
}





void Data1D::closestGridPoint(float x, unsigned int &p, float &l){
	unsigned int i_parameter0;
	
	l = cellLengthParameter0;
	
	i_parameter0 = static_cast<unsigned int>(floor((x-minParameter0)/cellLengthParameter0));
	
	
	p = i_parameter0;
}


void Data1D::print(){
	printf("no.:");
	for(unsigned int i=0; i<parameter0Dim; ++i) printf("\t%u",i);
	printf("\nheader:");
	for(unsigned int i=0; i<parameter0Dim; ++i) printf("\t%f",(*headerParameter0)[i]);
	printf("\ndata:");
	for(unsigned int i=0; i<parameter0Dim; ++i) printf("\t%f",(*data)[i]);
	printf("\n");
}












