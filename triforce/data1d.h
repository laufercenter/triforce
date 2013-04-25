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

typedef vec Vector;
//typedef Col<unsigned int> VectorInt;

typedef mat Matrix;

//#define THRESHOLD_NUMERICAL 0.00001

//typedef boost::numeric::ublas::vector<double> Vector;
//typedef boost::numeric::ublas::matrix<double> Matrix;



typedef boost::multi_array<double,1> Table1dDouble;


class Data1D{
	
public:
	Data1D();
	Data1D(unsigned int parameter0Dim, bool containsAuxiliaryData);
	void setHeaderParameter0Cell(unsigned int x, double value);
	double getHeaderParameter0Cell(unsigned int x);
	void setDataCell(unsigned int x, double value);
	void setAuxiliaryCell(unsigned int x, double value);
	void print();
	double getDataCell(unsigned int x);
	double getAuxiliaryCell(unsigned int x);
	
	Vector getHeaderVector(unsigned int parameter0);
	
	
	virtual void closestGridPoint(double x, unsigned int &p, double &l);
	bool isWithinNumericalLimits(double x, double t);
	void init();
	

	
	unsigned int parameter0Dim;
	double minParameter0;
	double cellLengthParameter0;
	Table1dDouble *headerParameter0;
	Table1dDouble *data;
	Table1dDouble *auxiliary;
	
	bool containsAuxiliaryData;
	
	
	
};

#endif //DATA1D_H_
