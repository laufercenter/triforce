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

#ifndef DATA6D_H_
#define DATA6D_H_

#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>

#include <armadillo>


using namespace std;
using namespace arma;

typedef fvec Vector;
typedef Col<unsigned int> VectorInt;

typedef fmat Matrix;

//#define THRESHOLD_NUMERICAL 0.00001

//typedef boost::numeric::ublas::vector<float> Vector;
//typedef boost::numeric::ublas::matrix<float> Matrix;



typedef boost::multi_array<float,1> Table1dDouble;
typedef boost::multi_array<float,2> Table2dDouble;
typedef boost::multi_array<float,6> Table6dDouble;
typedef boost::multi_array<Vector,6> Table6dVector;
typedef boost::multi_array<Matrix,6> Table6dMatrix;


class Data6D{
	
public:
	Data6D();
	Data6D(unsigned int parameter0Dim, unsigned int parameter1Dim, unsigned int parameter2Dim, unsigned int parameter3Dim, unsigned int parameter4Dim, unsigned int parameter5Dim);
	void setHeaderParameter0Cell(unsigned int x, float value);
	void setHeaderParameter1Cell(unsigned int x, float value);
	void setHeaderParameter2Cell(unsigned int x, float value);
	void setHeaderParameter3Cell(unsigned int x, float value);
	void setHeaderParameter4Cell(unsigned int x, float value);
	void setHeaderParameter5Cell(unsigned int x, float value);
	float getHeaderParameter0Cell(unsigned int x);
	float getHeaderParameter1Cell(unsigned int x);
	float getHeaderParameter2Cell(unsigned int x);
	float getHeaderParameter3Cell(unsigned int x);
	float getHeaderParameter4Cell(unsigned int x);
	float getHeaderParameter5Cell(unsigned int x);
	void setDataCell(unsigned int x, unsigned int y, unsigned int z, unsigned int u, unsigned int v, unsigned int w, float value);
	//void setAuxiliaryCell(unsigned int x, unsigned int y, unsigned int z, unsigned int u, unsigned int v, unsigned int w, float value);
	//void setGradientCell(unsigned int x, unsigned int y, unsigned int z, unsigned int u, unsigned int v, unsigned int w, unsigned int i, float value);
	//void setHessianCell(unsigned int x, unsigned int y, unsigned int z, unsigned int u, unsigned int v, unsigned int w, unsigned int i, unsigned int j, float value);
	void print();
	float getDataCell(unsigned int x, unsigned int y, unsigned int z, unsigned int u, unsigned int v, unsigned int w);
	//float getAuxiliaryCell(unsigned int x, unsigned int y, unsigned int z);
	//Vector &getGradient(unsigned int x, unsigned int y, unsigned int z);
	//Matrix &getHessian(unsigned int x, unsigned int y, unsigned int z);
	
	Vector getHeaderVector(unsigned int parameter0, unsigned int parameter1, unsigned int parameter2, unsigned int parameter3, unsigned int parameter4, unsigned int parameter5);
	//Vector bisectFloor(Vector &x);
	
	void printDataCell(unsigned int i, unsigned int j, unsigned int k);
	//void printGradientCell(unsigned int i, unsigned int j, unsigned int k);
	//void printHessianCell(unsigned int i, unsigned int j, unsigned int k);
	
	virtual void closestGridPoint(Vector &x, VectorInt &p, Vector &l);
	virtual void surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths);
	bool isWithinNumericalLimits(float x, float t);
	void init();
	

	
	unsigned int parameter0Dim;
	unsigned int parameter1Dim;
	unsigned int parameter2Dim;
	unsigned int parameter3Dim;
	unsigned int parameter4Dim;
	unsigned int parameter5Dim;
	float minParameter0;
	float minParameter1;
	float minParameter2;
	float minParameter3;
	float minParameter4;
	float minParameter5;
	float cellLengthParameter0;
	float cellLengthParameter1;
	float cellLengthParameter2;
	float cellLengthParameter3;
	float cellLengthParameter4;
	float cellLengthParameter5;
	Table1dDouble *headerParameter0;
	Table1dDouble *headerParameter1;
	Table1dDouble *headerParameter2;
	Table1dDouble *headerParameter3;
	Table1dDouble *headerParameter4;
	Table1dDouble *headerParameter5;
	Table6dDouble *data;
	Table6dVector *gradient;
	//Table6dMatrix *hessian;
	//Table6dDouble *auxiliary;
	//unsigned int derivativeLevel;
	//bool containsAuxiliaryData;
	
	
	
};

#endif //DATA6D_H_
