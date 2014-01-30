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
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02119, USA.
 */

#ifndef DATA1D_H_
#define DATA1D_H_

#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>

#include <armadillo>


using namespace std;
using namespace arma;

typedef fvec Vector;
//typedef Col<unsigned int> VectorInt;

typedef fmat Matrix;

//#define THRESHOLD_NUMERICAL 0.00001

//typedef boost::numeric::ublas::vector<float> Vector;
//typedef boost::numeric::ublas::matrix<float> Matrix;



typedef boost::multi_array<float,1> Table1dFloat;


class Data1D{
	
public:
	Data1D();
	Data1D(unsigned int parameter0Dim, bool containsAuxiliaryData);
	void setHeaderParameter0Cell(unsigned int x, float value);
	float getHeaderParameter0Cell(unsigned int x);
	void setDataCell(unsigned int x, float value);
	void setAuxiliaryCell(unsigned int x, float value);
	void print();
	float getDataCell(unsigned int x);
	float getAuxiliaryCell(unsigned int x);
	
	Vector getHeaderVector(unsigned int parameter0);
	
	
	virtual void closestGridPoint(float x, unsigned int &p, float &l);
	bool isWithinNumericalLimits(float x, float t);
	void init();
	

	
	unsigned int parameter0Dim;
	float minParameter0;
	float cellLengthParameter0;
	Table1dFloat *headerParameter0;
	Table1dFloat *data;
	Table1dFloat *auxiliary;
	
	bool containsAuxiliaryData;
	
	
	
};

#endif //DATA1D_H_
