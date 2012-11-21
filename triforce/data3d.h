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
typedef Col<int> VectorInt;

typedef mat Matrix;

//#define THRESHOLD_NUMERICAL 0.00001

//typedef boost::numeric::ublas::vector<double> Vector;
//typedef boost::numeric::ublas::matrix<double> Matrix;



typedef boost::multi_array<double,1> Table1dDouble;
typedef boost::multi_array<double,2> Table2dDouble;
typedef boost::multi_array<double,3> Table3dDouble;
typedef boost::multi_array<Vector,3> Table3dVector;
typedef boost::multi_array<Matrix,3> Table3dMatrix;


class Data3D{
	
public:
	Data3D(int PHIDim, int psiDim, int lambdaDim);
	void setHeaderPHICell(int x, double value);
	void setHeaderPsiCell(int x, double value);
	void setHeaderLambdaCell(int x, double value);
	void setDataCell(int x, int y, int z, double value);
	void setGradientCell(int x, int y, int z, int i, double value);
	void setHessianCell(int x, int y, int z, int i, int j, double value);
	void print();
	double getDataCell(int x, int y, int z);
	Vector &getGradient(int x, int y, int z);
	Matrix &getHessian(int x, int y, int z);
	
	Vector getHeaderVector(int PHI, int psi, int lambda);
	//Vector bisectFloor(Vector &x);
	Vector cellLength();
	double lambdaGridLength();
	double psiGridLength(int lambda);
	double PHIGridLength(int psi, int lambda);
	
	void printDataCell(int i, int j, int k);
	void printGradientCell(int i, int j, int k);
	void printHessianCell(int i, int j, int k);
	
	void closestGridPoint(Vector &x, VectorInt &p, Vector &l);
	void surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths);
	bool isWithinNumericalLimits(double x, double t);
	

	
private:
	
	int PHIDim;
	int psiDim;
	int lambdaDim;
	Table1dDouble *headerPHI;
	Table1dDouble *headerPsi;
	Table1dDouble *headerLambda;
	Table3dDouble *data;
	Table3dVector *gradient;
	Table3dMatrix *hessian;
	
	
	
};

#endif //DATA3D_H_
