/* Copyright 2012, Nils J. D. Drechsel & Jordi Villà-Freixa
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef DATA3D_H_
#define DATA3D_H_

#include "data.h"
#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <armadillo>
#include <algorithm>
#include <limits>
#include <math.h>


using namespace boost;


using namespace std;
using namespace arma;

typedef fvec Vector;
typedef Col<unsigned int> VectorInt;

typedef fmat Matrix;




template<class T>
class Data3D: public Data<T>{
	
	
typedef boost::multi_array<float,1> Table1dDouble;
typedef boost::multi_array<float,2> Table2dDouble;
typedef boost::multi_array<T,3> Table3dDouble;
typedef boost::multi_array<Vector,3> Table3dVector;
typedef boost::multi_array<Matrix,3> Table3dMatrix;
	
public:
	Data3D();
	Data3D(VectorInt &parameterDim, unsigned int derivativeLevel, bool containsAuxiliaryData);
	void setHeaderCell(unsigned int p,unsigned int x, float value);
	float getHeaderCell(unsigned int p,unsigned int x);
	void setDataCell(VectorInt &x, T value);
	void setAuxiliaryCell(VectorInt &x, T value);
	void setGradientCell(VectorInt &x, unsigned int i, float value);
	void setHessianCell(VectorInt &x, unsigned int i, unsigned int j, float value);
	void print();
	T getDataCell(VectorInt &x);
	T getDataCell(VectorInt &x, bool checkLimits);
	T getAuxiliaryCell(VectorInt &x);
	Vector &getGradient(VectorInt &x);
	Matrix &getHessian(VectorInt &x);
	
	Vector getHeaderVector(VectorInt &p);
	
/*	void printDataCell(unsigned int i, unsigned int j, unsigned int k);
	void printGradientCell(unsigned int i, unsigned int j, unsigned int k);
	void printHessianCell(unsigned int i, unsigned int j, unsigned int k);
*/	

	void closestGridPoint(Vector &x, VectorInt &p, Vector &d);
	void init();
	VectorInt getDimensions();	
	Vector getReciprocalCellLengths();	

	
	unsigned int parameter0Dim;
	unsigned int parameter1Dim;
	unsigned int parameter2Dim;
	float minParameter0;
	float minParameter1;
	float minParameter2;
	float cellLengthParameter0;
	float cellLengthParameter1;
	float cellLengthParameter2;
	float reciprocalCellLengthParameter0;
	float reciprocalCellLengthParameter1;
	float reciprocalCellLengthParameter2;
	Table1dDouble *headerParameter0;
	Table1dDouble *headerParameter1;
	Table1dDouble *headerParameter2;
	Table3dDouble *data;
	Table3dVector *gradient;
	Table3dMatrix *hessian;
	Table3dDouble *auxiliary;
	unsigned int derivativeLevel;
	bool containsAuxiliaryData;
	
	
	
};




template <class T>
Data3D<T>::Data3D(){
}

template <class T>
Data3D<T>::Data3D(VectorInt &parameterDim, unsigned int derivativeLevel, bool containsAuxiliaryData){
	
	this->parameter0Dim = parameterDim(0);
	this->parameter1Dim = parameterDim(1);
	this->parameter2Dim = parameterDim(2);
	
	this->derivativeLevel = derivativeLevel;
	this->containsAuxiliaryData = containsAuxiliaryData;
		
		
	headerParameter0 = new Table1dDouble(boost::extents[parameter0Dim]);
	headerParameter1 = new Table1dDouble(boost::extents[parameter1Dim]);
	headerParameter2 = new Table1dDouble(boost::extents[parameter2Dim]);
	data = new Table3dDouble(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim]);
	if(containsAuxiliaryData){
		auxiliary = new Table3dDouble(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim]);
	}
	if(derivativeLevel>=2){
		gradient = new Table3dVector(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim]);
		if(derivativeLevel>=3){
			hessian = new Table3dMatrix(boost::extents[parameter0Dim][parameter1Dim][parameter2Dim]);
		}
		
		for(unsigned int x=0; x<parameter0Dim; x++)
			for(unsigned int y=0; y<parameter1Dim; y++)
				for(unsigned int z=0; z<parameter2Dim; z++){
					boost::array<Table3dVector::index,3> idx = {{x,y,z}};
					(*gradient)(idx) = Vector(3);
					if(derivativeLevel>=3){
						(*hessian)(idx) = Matrix(3,3);
					}
				}
		
		
	}
}


template <class T>
VectorInt Data3D<T>::getDimensions(){
	VectorInt v(3);
	v(0)=parameter0Dim;
	v(1)=parameter1Dim;
	v(2)=parameter2Dim;
	return v;
}

template <class T>
Vector Data3D<T>::getReciprocalCellLengths(){
	Vector v(3);
	v(0)=reciprocalCellLengthParameter0;
	v(1)=reciprocalCellLengthParameter1;
	v(2)=reciprocalCellLengthParameter2;
	return v;
}

template <class T>
void Data3D<T>::init(){
	cellLengthParameter0 = abs((*headerParameter0)[1]-(*headerParameter0)[0]);
	cellLengthParameter1 = abs((*headerParameter1)[1]-(*headerParameter1)[0]);
	cellLengthParameter2 = abs((*headerParameter2)[1]-(*headerParameter2)[0]);

	reciprocalCellLengthParameter0 = 1.0/abs((*headerParameter0)[1]-(*headerParameter0)[0]);
	reciprocalCellLengthParameter1 = 1.0/abs((*headerParameter1)[1]-(*headerParameter1)[0]);
	reciprocalCellLengthParameter2 = 1.0/abs((*headerParameter2)[1]-(*headerParameter2)[0]);

	
	minParameter0 = (*headerParameter0)[0];
	minParameter1 = (*headerParameter1)[0];
	minParameter2 = (*headerParameter2)[0];
}

template <class T>
void Data3D<T>::setHeaderCell(unsigned int p, unsigned int x, float value){
	switch(p){
		case 0:	(*headerParameter0)[x] = value; break;
		case 1:	(*headerParameter1)[x] = value; break;
		case 2:	(*headerParameter2)[x] = value; break;
	}
}


template <class T>
float Data3D<T>::getHeaderCell(unsigned int p, unsigned int x){
	switch(p){
		case 0:	return (*headerParameter0)[x]; break;
		case 1:	return (*headerParameter1)[x]; break;
		case 2: return (*headerParameter2)[x]; break;
	}	
	return 0;
}


template <class T>
void Data3D<T>::setDataCell(VectorInt &x, T value){
	(*data)[x(0)][x(1)][x(2)] = value;
}


template <class T>
void Data3D<T>::setAuxiliaryCell(VectorInt &x, T value){
	(*auxiliary)[x(0)][x(1)][x(2)] = value;
}


template <class T>
void Data3D<T>::setGradientCell(VectorInt &x, unsigned int i, float value){
	(*gradient)[x(0)][x(1)][x(2)](i) = value;
}

template <class T>
void Data3D<T>::setHessianCell(VectorInt &x, unsigned int i, unsigned int j, float value){
	(*hessian)[x(0)][x(1)][x(2)](i,j) = value;
}

template <class T>
Vector Data3D<T>::getHeaderVector(VectorInt &p){
	//printf("HEADER VECTOR: %d %d %d\n",parameter0,parameter1,parameter2);
	Vector v=Vector(3);
	v(0) = (*headerParameter0)[p(0)];
	v(1) = (*headerParameter1)[p(1)];
	v(2) = (*headerParameter2)[p(2)];
	
	return v;
}




template <class T>
T Data3D<T>::getDataCell(VectorInt &x){
	return (*data)[x(0)][x(1)][x(2)];
}

template <class T>
T Data3D<T>::getDataCell(VectorInt &x, bool checkLimits){
	
	return (*data)[min(parameter0Dim-1,max((unsigned int)0,x(0)))][min(parameter1Dim-1,max((unsigned int)0,x(1)))][min(parameter2Dim-1,max((unsigned int)0,x(2)))];
}


template <class T>
T Data3D<T>::getAuxiliaryCell(VectorInt &x){
	return (*auxiliary)[x(0)][x(1)][x(2)];
}


template <class T>
Vector &Data3D<T>::getGradient(VectorInt &x){
	return (*gradient)[x(0)][x(1)][x(2)];
}


template <class T>
Matrix &Data3D<T>::getHessian(VectorInt &x){
	return (*hessian)[x(0)][x(1)][x(2)];
}








template <class T>
void Data3D<T>::closestGridPoint(Vector &x, VectorInt &p, Vector &d){
	unsigned int i_parameter0;
	unsigned int i_parameter1;
	unsigned int i_parameter2;
	
	//printf("CELL LENGTHS: %f %f %f\n",l(0),l(1),l(2));
	
	//printf("X: %f %f %f, l: %f %f %f, min: %f %f %f, celll: %f %f %f\n",x(0),x(1),x(2),l(0),l(1),l(2),minParameter0, minParameter1, minParameter2, cellLengthParameter0, cellLengthParameter1, cellLengthParameter2);
	
	d=Vector(3);
	d(0) = (x(0)-minParameter0)*reciprocalCellLengthParameter0;
	d(1) = (x(1)-minParameter1)*reciprocalCellLengthParameter1;
	d(2) = (x(2)-minParameter2)*reciprocalCellLengthParameter2;
	
	i_parameter0 = static_cast<unsigned int>(floor(d(0)));
	
	i_parameter1 = static_cast<unsigned int>(floor(d(1)));
	
	i_parameter2 = static_cast<unsigned int>(floor(d(2)));
	
	d(0) = d(0)-i_parameter0;
	d(1) = d(1)-i_parameter1;
	d(2) = d(2)-i_parameter2;
	
	//printf("CLOSESTGRRRID: %f %f %f\n",x(2),x(2)-minParameter2,(x(2)-minParameter2)/cellLengthParameter2);
	
	p=VectorInt(3);
	p(0) = min(parameter0Dim-1, max((unsigned int) 0, i_parameter0));
	p(1) = min(parameter1Dim-1, max((unsigned int) 0, i_parameter1));
	p(2) = min(parameter2Dim-1, max((unsigned int) 0, i_parameter2));
}





template <class T>
void Data3D<T>::print(){
	printf("header dim0");
	for(unsigned int i=0; i<parameter0Dim; ++i) printf("\t%f",(*headerParameter0)[i]);
	printf("\nheader dim1");
	for(unsigned int j=0; j<parameter1Dim; ++j) printf("\t%f",(*headerParameter1)[j]);
	printf("\nheader dim2");
	for(unsigned int k=0; k<parameter2Dim; ++k) printf("\t%f",(*headerParameter2)[k]);
}


/*
template <class T>
void Data3D<T>::printDataCell(unsigned int i, unsigned int j, unsigned int k){
	printf("Cell[%d,%d,%d]: %f\n",i,j,k,(*data)[i][j][k]);
}


template <class T>
void Data3D<T>::printGradientCell(unsigned int i, unsigned int j, unsigned int k){
	printf("Gradient[%d,%d,%d]: (%f, %f, %f)\n",i,j,k,(*gradient)[i][j][k](0),(*gradient)[i][j][k](1),(*gradient)[i][j][k](2));
}


template <class T>
void Data3D<T>::printHessianCell(unsigned int i, unsigned int j, unsigned int k){
	printf("Hessian[%d,%d,%d]: \n",i,j,k);
	printf("|%f, %f, %f|\n",(*hessian)[i][j][k](0,0),(*hessian)[i][j][k](0,1),(*hessian)[i][j][k](0,2));
	printf("|%f, %f, %f|\n",(*hessian)[i][j][k](1,0),(*hessian)[i][j][k](1,1),(*hessian)[i][j][k](1,2));
	printf("|%f, %f, %f|\n",(*hessian)[i][j][k](2,0),(*hessian)[i][j][k](2,1),(*hessian)[i][j][k](2,2));
}
*/

#endif //DATA3D_H_
