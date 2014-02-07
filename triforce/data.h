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

#ifndef DATA_H_
#define DATA_H_

#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>

#include <armadillo>
#define PURE =0

using namespace std;
using namespace arma;

typedef fvec Vector;
typedef Col<unsigned int> VectorInt;

typedef fmat Matrix;

//#define THRESHOLD_NUMERICAL 0.00001

//typedef boost::numeric::ublas::vector<float> Vector;
//typedef boost::numeric::ublas::matrix<float> Matrix;




template<class T>
class Data{
	
	
public:
	//Data(VectorInt &parameterDim, unsigned int derivativeLevel, bool containsAuxiliaryData) PURE;
	virtual void setHeaderCell(unsigned int p,unsigned int x, float value) PURE;
	virtual float getHeaderCell(unsigned int p,unsigned int x) PURE;
	virtual void setDataCell(VectorInt &x, T value) PURE;
	virtual void setAuxiliaryCell(VectorInt &x, T value) PURE;
	virtual void setGradientCell(VectorInt &x, unsigned int i, float value) PURE;
	virtual void setHessianCell(VectorInt &x, unsigned int i, unsigned int j, float value) PURE;
	virtual void print() PURE;
	virtual T getDataCell(VectorInt &x) PURE;
	virtual T getDataCell(VectorInt &x, bool checkLimits) PURE;
	virtual T getAuxiliaryCell(VectorInt &x) PURE; 
	virtual Vector &getGradient(VectorInt &x) PURE;
	virtual Matrix &getHessian(VectorInt &x) PURE;
	
	virtual Vector getHeaderVector(VectorInt &p) PURE;
	
/*	void printDataCell(unsigned int i, unsigned int j, unsigned int k);
	void printGradientCell(unsigned int i, unsigned int j, unsigned int k);
	void printHessianCell(unsigned int i, unsigned int j, unsigned int k);
*/	

	virtual void closestGridPoint(Vector &x, VectorInt &p, Vector &d) PURE;
	virtual void init() PURE;
	virtual VectorInt getDimensions() PURE;
	virtual Vector getReciprocalCellLengths() PURE;
	

	
	
	
	
};


#endif //DATA_H_
