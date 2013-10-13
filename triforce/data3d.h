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

#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>

#include <armadillo>


using namespace std;
using namespace arma;

typedef vec Vector;
typedef Col<unsigned int> VectorInt;

typedef fmat Matrix;

//#define THRESHOLD_NUMERICAL 0.00001

//typedef boost::numeric::ublas::vector<float> Vector;
//typedef boost::numeric::ublas::matrix<float> Matrix;



typedef boost::multi_array<float,1> Table1dDouble;
typedef boost::multi_array<float,2> Table2dDouble;
typedef boost::multi_array<float,3> Table3dDouble;
typedef boost::multi_array<Vector,3> Table3dVector;
typedef boost::multi_array<Matrix,3> Table3dMatrix;


class Data3D{
	
public:
	Data3D();
	Data3D(unsigned int parameter0Dim, unsigned int parameter1Dim, unsigned int parameter2Dim, unsigned int derivativeLevel, bool containsAuxiliaryData);
	void setHeaderParameter0Cell(unsigned int x, double value);
	void setHeaderParameter1Cell(unsigned int x, double value);
	void setHeaderParameter2Cell(unsigned int x, double value);
	double getHeaderParameter0Cell(unsigned int x);
	double getHeaderParameter1Cell(unsigned int x);
	double getHeaderParameter2Cell(unsigned int x);
	void setDataCell(unsigned int x, unsigned int y, unsigned int z, double value);
	void setAuxiliaryCell(unsigned int x, unsigned int y, unsigned int z, double value);
	void setGradientCell(unsigned int x, unsigned int y, unsigned int z, unsigned int i, double value);
	void setHessianCell(unsigned int x, unsigned int y, unsigned int z, unsigned int i, unsigned int j, double value);
	void print();
	double getDataCell(unsigned int x, unsigned int y, unsigned int z);
	double getAuxiliaryCell(unsigned int x, unsigned int y, unsigned int z);
	Vector &getGradient(unsigned int x, unsigned int y, unsigned int z);
	Matrix &getHessian(unsigned int x, unsigned int y, unsigned int z);
	
	Vector getHeaderVector(unsigned int parameter0, unsigned int parameter1, unsigned int parameter2);
	//Vector bisectFloor(Vector &x);
	double parameter2GridLength();
	double parameter1GridLength(unsigned int parameter2);
	double parameter0GridLength(unsigned int parameter1, unsigned int parameter2);
	
	void printDataCell(unsigned int i, unsigned int j, unsigned int k);
	void printGradientCell(unsigned int i, unsigned int j, unsigned int k);
	void printHessianCell(unsigned int i, unsigned int j, unsigned int k);
	
	virtual void closestGridPoint(Vector &x, VectorInt &p, Vector &l);
	virtual void surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths);
	bool isWithinNumericalLimits(float x, float t);
	void init();
	

	
	unsigned int parameter0Dim;
	unsigned int parameter1Dim;
	unsigned int parameter2Dim;
	double minParameter0;
	double minParameter1;
	double minParameter2;
	double cellLengthParameter0;
	double cellLengthParameter1;
	double cellLengthParameter2;
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

#endif //DATA3D_H_
