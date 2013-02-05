#include "data3d.h"

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;

#define THRESHOLD_NUMERICAL 0.00001

Data3D::Data3D(){
}

Data3D::Data3D(int parameter0Dim, int parameter1Dim, int parameter2Dim, int derivativeLevel, bool containsAuxiliaryData){
	int maxdim=0;
	this->parameter0Dim = parameter0Dim;
	this->parameter1Dim = parameter1Dim;
	this->parameter2Dim = parameter2Dim;
	
	this->derivativeLevel = derivativeLevel;
	this->containsAuxiliaryData = containsAuxiliaryData;
		
		
	headerParameter0 = new Table1dDouble(boost::extents[parameter0Dim]);
	headerParameter1 = new Table1dDouble(boost::extents[parameter1Dim]);
	headerParameter2 = new Table1dDouble(boost::extents[parameter2Dim]);
	data = new Table3dDouble(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim]);
	auxiliary = new Table3dDouble(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim]);
	gradient = new Table3dVector(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim]);
	hessian = new Table3dMatrix(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim]);
	for(int x=0; x<parameter0Dim; x++)
		for(int y=0; y<parameter1Dim; y++)
			for(int z=0; z<parameter2Dim; z++){
				boost::array<Table3dVector::index,3> idx = {{x,y,z}};
				(*gradient)(idx) = Vector(3);
				(*hessian)(idx) = Matrix(3,3);
			}
}


void Data3D::init(){
	cellLengthParameter0 = abs((*headerParameter0)[1]-(*headerParameter0)[0]);
	cellLengthParameter1 = abs((*headerParameter1)[1]-(*headerParameter1)[0]);
	cellLengthParameter2 = abs((*headerParameter2)[1]-(*headerParameter2)[0]);
	
	minParameter0 = (*headerParameter0)[0];
	minParameter1 = (*headerParameter1)[0];
	minParameter2 = (*headerParameter2)[0];
}


void Data3D::setHeaderParameter0Cell(int x, double value){
	(*headerParameter0)[x] = value;
}

void Data3D::setHeaderParameter1Cell(int x, double value){
	(*headerParameter1)[x] = value;
}

void Data3D::setHeaderParameter2Cell(int x, double value){
	(*headerParameter2)[x] = value;
}


void Data3D::setDataCell(int x, int y, int z, double value){
	(*data)[x][y][z] = value;
}

void Data3D::setAuxiliaryCell(int x, int y, int z, double value){
	(*auxiliary)[x][y][z] = value;
}

void Data3D::setGradientCell(int x, int y, int z, int i, double value){
	(*gradient)[x][y][z](i) = value;
}
void Data3D::setHessianCell(int x, int y, int z, int i, int j, double value){
	(*hessian)[x][y][z](i,j) = value;
}

Vector Data3D::getHeaderVector(int parameter0, int parameter1, int parameter2){
	//printf("HEADER VECTOR: %d %d %d\n",parameter0,parameter1,parameter2);
	Vector v=Vector(3);
	v(0) = (*headerParameter0)[parameter0];
	v(1) = (*headerParameter1)[parameter1];
	v(2) = (*headerParameter2)[parameter2];
	
	return v;
}



double Data3D::parameter2GridLength(){
	return abs((*headerParameter2)[parameter2Dim-1]-(*headerParameter2)[0]);
}

double Data3D::parameter1GridLength(int parameter2){
	return abs((*headerParameter1)[parameter1Dim-1]-(*headerParameter1)[0]);
}

double Data3D::parameter0GridLength(int parameter1, int parameter2){
	return abs((*headerParameter0)[parameter0Dim-1]-(*headerParameter0)[0]);
}



double Data3D::getDataCell(int x, int y, int z){
	return (*data)[x][y][z];
}

double Data3D::getAuxiliaryCell(int x, int y, int z){
	return (*auxiliary)[x][y][z];
}

Vector &Data3D::getGradient(int x, int y, int z){
	return (*gradient)[x][y][z];
}

Matrix &Data3D::getHessian(int x, int y, int z){
	return (*hessian)[x][y][z];
}

bool Data3D::isWithinNumericalLimits(double x, double t){
	if(abs(x-t) <= THRESHOLD_NUMERICAL) return true;
	else return false;
}





void Data3D::closestGridPoint(Vector &x, VectorInt &p, Vector &l){
	double lengthparameter0;
	int i_parameter0;
	double lengthparameter1;
	int i_parameter1;
	double lengthparameter2;
	int i_parameter2;
	
	l(0) = cellLengthParameter0;
	l(1) = cellLengthParameter1;
	l(2) = cellLengthParameter2;
	
	//printf("CELL LENGTHS: %f %f %f\n",l(0),l(1),l(2));
	
	//printf("X: %f %f %f, l: %f %f %f, min: %f %f %f, celll: %f %f %f\n",x(0),x(1),x(2),l(0),l(1),l(2),minParameter0, minParameter1, minParameter2, cellLengthParameter0, cellLengthParameter1, cellLengthParameter2);
	
	i_parameter0 = static_cast<int>(floor((x(0)-minParameter0)/cellLengthParameter0));
	
	i_parameter1 = static_cast<int>(floor((x(1)-minParameter1)/cellLengthParameter1));
	
	i_parameter2 = static_cast<int>(floor((x(2)-minParameter2)/cellLengthParameter2));
	
	//printf("CLOSESTGRRRID: %f %f %f\n",x(2),x(2)-minParameter2,(x(2)-minParameter2)/cellLengthParameter2);
	
	p(0) = i_parameter0;
	p(1) = i_parameter1;
	p(2) = i_parameter2;
}



void Data3D::surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths){
	VectorInt v2(3);
	int i,j,k;
	VectorInt v(3);
	bool neg;
	
	closestGridPoint(x, v, lengths);
	
	
	//printf("CLOSEST GRIDPOINT: %d %d %d (%f, %f, %f)\n",v(0),v(1),v(2),(*headerParameter0)[v(0)],(*headerParameter1)[v(1)],(*headerParameter2)[v(2)]);
	//printf("CLOSEST GRIDPOINT: %d %d %d\n",v(0),v(1),v(2));
	
	
	r.clear();
	for(i=0;i<2;++i)
		for(j=0;j<2;++j)
			for(k=0;k<2;++k){
				v2(0)=v(0)+i;
				v2(1)=v(1)+j;
				v2(2)=v(2)+k;
				if(v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim){
					if(!isnan((*data)[v2(0)][v2(1)][v2(2)])){
							r.push_back(v2);
							printf("ACCEPTED\n");
						
					}
					else printf("REJECTED NAN\n");
				}
				else printf("REJECTED OUT OF LIMIT\n");
			}
			
}


void Data3D::printDataCell(int i, int j, int k){
	printf("Cell[%d,%d,%d]: %f\n",i,j,k,(*data)[i][j][k]);
}

void Data3D::printGradientCell(int i, int j, int k){
	printf("Gradient[%d,%d,%d]: (%f, %f, %f)\n",i,j,k,(*gradient)[i][j][k](0),(*gradient)[i][j][k](1),(*gradient)[i][j][k](2));
}

void Data3D::printHessianCell(int i, int j, int k){
	printf("Hessian[%d,%d,%d]: \n");
	printf("|%f, %f, %f|\n",i,j,k,(*hessian)[i][j][k](0,0),(*hessian)[i][j][k](0,1),(*hessian)[i][j][k](0,2));
	printf("|%f, %f, %f|\n",i,j,k,(*hessian)[i][j][k](1,0),(*hessian)[i][j][k](1,1),(*hessian)[i][j][k](1,2));
	printf("|%f, %f, %f|\n",i,j,k,(*hessian)[i][j][k](2,0),(*hessian)[i][j][k](2,1),(*hessian)[i][j][k](2,2));
}
