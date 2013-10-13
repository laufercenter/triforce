#include "data6d.h"

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;

#define THRESHOLD_NUMERICAL 0.00001

Data6D::Data6D(){
}

Data6D::Data6D(unsigned int parameter0Dim, unsigned int parameter1Dim, unsigned int parameter2Dim, unsigned int parameter3Dim, unsigned int parameter4Dim, unsigned int parameter5Dim){
	this->parameter0Dim = parameter0Dim;
	this->parameter1Dim = parameter1Dim;
	this->parameter2Dim = parameter2Dim;
	this->parameter3Dim = parameter3Dim;
	this->parameter4Dim = parameter4Dim;
	this->parameter5Dim = parameter5Dim;
	
	headerParameter0 = new Table1dDouble(boost::extents[parameter0Dim]);
	headerParameter1 = new Table1dDouble(boost::extents[parameter1Dim]);
	headerParameter2 = new Table1dDouble(boost::extents[parameter2Dim]);
	headerParameter2 = new Table1dDouble(boost::extents[parameter3Dim]);
	headerParameter2 = new Table1dDouble(boost::extents[parameter4Dim]);
	headerParameter2 = new Table1dDouble(boost::extents[parameter5Dim]);
	data = new Table6dDouble(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim][parameter3Dim][parameter4Dim][parameter5Dim]);
}


void Data6D::init(){
	cellLengthParameter0 = abs((*headerParameter0)[1]-(*headerParameter0)[0]);
	cellLengthParameter1 = abs((*headerParameter1)[1]-(*headerParameter1)[0]);
	cellLengthParameter2 = abs((*headerParameter2)[1]-(*headerParameter2)[0]);
	cellLengthParameter3 = abs((*headerParameter3)[1]-(*headerParameter3)[0]);
	cellLengthParameter4 = abs((*headerParameter4)[1]-(*headerParameter4)[0]);
	cellLengthParameter5 = abs((*headerParameter5)[1]-(*headerParameter5)[0]);
	
	minParameter0 = (*headerParameter0)[0];
	minParameter1 = (*headerParameter1)[0];
	minParameter2 = (*headerParameter2)[0];
	minParameter3 = (*headerParameter3)[0];
	minParameter4 = (*headerParameter4)[0];
	minParameter5 = (*headerParameter5)[0];
}


void Data6D::setHeaderParameter0Cell(unsigned int x, float value){
	(*headerParameter0)[x] = value;
}

void Data6D::setHeaderParameter1Cell(unsigned int x, float value){
	(*headerParameter1)[x] = value;
}

void Data6D::setHeaderParameter2Cell(unsigned int x, float value){
	(*headerParameter2)[x] = value;
}

void Data6D::setHeaderParameter3Cell(unsigned int x, float value){
	(*headerParameter3)[x] = value;
}

void Data6D::setHeaderParameter4Cell(unsigned int x, float value){
	(*headerParameter4)[x] = value;
}

void Data6D::setHeaderParameter5Cell(unsigned int x, float value){
	(*headerParameter5)[x] = value;
}

float Data6D::getHeaderParameter0Cell(unsigned int x){
	return (*headerParameter0)[x];
}

float Data6D::getHeaderParameter1Cell(unsigned int x){
	return (*headerParameter1)[x];
}

float Data6D::getHeaderParameter2Cell(unsigned int x){
	return (*headerParameter2)[x];
}

float Data6D::getHeaderParameter3Cell(unsigned int x){
	return (*headerParameter3)[x];
}

float Data6D::getHeaderParameter4Cell(unsigned int x){
	return (*headerParameter4)[x];
}

float Data6D::getHeaderParameter5Cell(unsigned int x){
	return (*headerParameter5)[x];
}

void Data6D::setDataCell(unsigned int x, unsigned int y, unsigned int z, unsigned int u, unsigned int v, unsigned int w, float value){
	(*data)[x][y][z][u][v][w] = value;
}



Vector Data6D::getHeaderVector(unsigned int parameter0, unsigned int parameter1, unsigned int parameter2, unsigned int parameter3, unsigned int parameter4, unsigned int parameter5){
	//printf("HEADER VECTOR: %d %d %d\n",parameter0,parameter1,parameter2);
	Vector v=Vector(6);
	v(0) = (*headerParameter0)[parameter0];
	v(1) = (*headerParameter1)[parameter1];
	v(2) = (*headerParameter2)[parameter2];
	v(3) = (*headerParameter2)[parameter3];
	v(4) = (*headerParameter2)[parameter4];
	v(5) = (*headerParameter2)[parameter5];
	
	return v;
}





float Data6D::getDataCell(unsigned int x, unsigned int y, unsigned int z, unsigned int u, unsigned int v, unsigned int w){
	return (*data)[x][y][z][u][v][w];
}

bool Data6D::isWithinNumericalLimits(float x, float t){
	if(abs(x-t) <= THRESHOLD_NUMERICAL) return true;
	else return false;
}





void Data6D::closestGridPoint(Vector &x, VectorInt &p, Vector &l){
	unsigned int i_parameter0;
	unsigned int i_parameter1;
	unsigned int i_parameter2;
	unsigned int i_parameter3;
	unsigned int i_parameter4;
	unsigned int i_parameter5;
	
	l(0) = cellLengthParameter0;
	l(1) = cellLengthParameter1;
	l(2) = cellLengthParameter2;
	l(3) = cellLengthParameter3;
	l(4) = cellLengthParameter4;
	l(5) = cellLengthParameter5;
	
	//printf("CELL LENGTHS: %f %f %f\n",l(0),l(1),l(2));
	
	//printf("X: %f %f %f, l: %f %f %f, min: %f %f %f, celll: %f %f %f\n",x(0),x(1),x(2),l(0),l(1),l(2),minParameter0, minParameter1, minParameter2, cellLengthParameter0, cellLengthParameter1, cellLengthParameter2);
	
	i_parameter0 = static_cast<unsigned int>(floor((x(0)-minParameter0)/cellLengthParameter0));
	
	i_parameter1 = static_cast<unsigned int>(floor((x(1)-minParameter1)/cellLengthParameter1));
	
	i_parameter2 = static_cast<unsigned int>(floor((x(2)-minParameter2)/cellLengthParameter2));

	i_parameter3 = static_cast<unsigned int>(floor((x(3)-minParameter3)/cellLengthParameter3));

	i_parameter4 = static_cast<unsigned int>(floor((x(4)-minParameter4)/cellLengthParameter4));

	i_parameter5 = static_cast<unsigned int>(floor((x(5)-minParameter5)/cellLengthParameter5));
	
	//printf("CLOSESTGRRRID: %f %f %f\n",x(2),x(2)-minParameter2,(x(2)-minParameter2)/cellLengthParameter2);
	
	p(0) = i_parameter0;
	p(1) = i_parameter1;
	p(2) = i_parameter2;
	p(3) = i_parameter3;
	p(4) = i_parameter4;
	p(5) = i_parameter5;
}



void Data6D::surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths){
	VectorInt v2(6);
	Vector v3(64);
	float d;
	VectorInt v(6);
	unsigned int discontinuity;
	float vpsi,vlambda;
	Vector excess(6);
	closestGridPoint(x, v, lengths);
	Vector p(6);
	
	VectorInt r0,r1(6);
	
	//printf("CLOSEST GRIDPOINT: %d %d %d\n",v(0),v(1),v(2));
	
	if(v(0)>=parameter0Dim || v(1)>=parameter1Dim || v(2)>=parameter2Dim || v(3)>=parameter3Dim || v(4)>=parameter4Dim || v(5)>=parameter5Dim) return;
	
	
	for(unsigned int i=0; i<6; ++i)
		p(i) = x(i) / lengths(i);
	
	excess = p-v;
	
	for(unsigned int i=0; i<6; ++i)
		if(excess(i)<=0.5) r0(i)=0;
			else r0(i)=1;
			
	r.push_back(r0);
	for(unsigned int i=0; i<6; ++i){
		if(r0(i)==0){
			r0(i)=1;
			push_back(v+r0);
			r0(i)=0;
		}
		else{
			r0(i)=0;
			push_back(v+r0);
			r0(i)=1;
		}
	}
	
	return;
	
	/*
	
	
	//decide on which side of the discontinuities x lies.
	
	if(x(1) < x(2)){
		discontinuity = 0;
	}
	else if(x(1)+x(2) < M_PI){
		discontinuity = 1;
	}
	else{
		discontinuity = 2;
	}
	
	
	
	
	//printf("DISCONTINUITY: %d\n",discontinuity);
	
	
	r.clear();
	for(unsigned int i=0;i<2;++i)
		for(unsigned int j=0;j<2;++j)
			for(unsigned int k=0;k<2;++k)
				for(unsigned int l=0;l<2;++l)
					for(unsigned int m=0;m<2;++m)
						for(unsigned int n=0;n<2;++n){
							v2(0)=v(0)+i;
							v2(1)=v(1)+j;
							v2(2)=v(2)+k;
							v2(3)=v(2)+l;
							v2(4)=v(2)+m;
							v2(5)=v(2)+n;
							//printf("v2: %d %d %d, par: %d %d %d, bool: %d\n",v2(0),v2(1),v2(2),parameter0Dim,parameter1Dim,parameter2Dim, v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim);
							if(v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim && v2(3)<parameter3Dim && v2(4)<parameter4Dim && v2(5)<parameter5Dim){
								//if(!isnan((*data)[v2(0)][v2(1)][v2(2)][v2(3)][v2(4)][v2(5)])){
									vpsi = (*headerParameter1)[v2(1)];
									vlambda = (*headerParameter2)[v2(2)];
									
									if(discontinuity==0){
										if(vpsi < vlambda) r.push_back(v2);
									}
									else if(discontinuity==1){
										if(vpsi >= vlambda && vpsi + vlambda < M_PI && !isWithinNumericalLimits(vpsi + vlambda,M_PI)) r.push_back(v2);
									}
									else{
										if(vpsi >= vlambda && vpsi + vlambda >= M_PI) r.push_back(v2);
									}
										//printf("ACCEPTED %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));
									
									//else	printf("REJECTED 0 %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));

								//}
								//else	printf("REJECTED 1 %d %d %d\n",v2(0),v2(1),v2(2));

							}
							//else	printf("PRE-REJECTED 2 %d %d %d\n",v2(0),v2(1),v2(2));
						}
					
				
			
			
			
			
	if(r.size()==0){
		for(unsigned int i=0;i<2;++i)
			for(unsigned int j=0;j<2;++j)
				for(unsigned int k=0;k<2;++k)
					for(unsigned int l=0;l<2;++l)
						for(unsigned int m=0;m<2;++m)
							for(unsigned int n=0;n<2;++n){
					
								v2(0)=v(0)+i;
								v2(1)=v(1)+j;
								v2(2)=v(2)+k;
								v2(3)=v(2)+l;
								v2(4)=v(2)+m;
								v2(5)=v(2)+n;
								
								if(v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim && v2(3)<parameter3Dim && v2(4)<parameter4Dim && v2(5)<parameter5Dim){
									//if(!isnan((*data)[v2(0)][v2(1)][v2(2)][v2(3)][v2(4)][v2(5)])){
										//printf("RE-ACCEPTED %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));
											r.push_back(v2);
									}

								}
							}
		
	}

	*/		
			
}

void Data6D::print(){
	printf("header dim0");
	for(unsigned int i=0; i<parameter0Dim; ++i) printf("\t%f",(*headerParameter0)[i]);
	printf("\nheader dim1");
	for(unsigned int j=0; j<parameter1Dim; ++j) printf("\t%f",(*headerParameter1)[j]);
	printf("\nheader dim2");
	for(unsigned int k=0; k<parameter2Dim; ++k) printf("\t%f",(*headerParameter2)[k]);
	printf("\nheader dim3");
	for(unsigned int k=0; k<parameter3Dim; ++k) printf("\t%f",(*headerParameter3)[k]);
	printf("\nheader dim4");
	for(unsigned int k=0; k<parameter4Dim; ++k) printf("\t%f",(*headerParameter4)[k]);
	printf("\nheader dim5");
	for(unsigned int k=0; k<parameter5Dim; ++k) printf("\t%f",(*headerParameter5)[k]);
}


void Data6D::printDataCell(unsigned int i, unsigned int j, unsigned int k, unsigned int l, unsigned int m, unsigned int n){
	printf("Cell[%d,%d,%d,%d,%d,%d]: %f\n",i,j,k,l,m,n,(*data)[i][j][k][l][m][n]);
}

