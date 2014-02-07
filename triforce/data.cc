

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;

#define THRESHOLD_NUMERICAL 0.00001

template <class T>
Data3D<T>::Data3D(){
}

template <class T>
Data3D<T>::Data3D(unsigned int parameter0Dim, unsigned int parameter1Dim, unsigned int parameter2Dim, unsigned int derivativeLevel, bool containsAuxiliaryData){
	this->parameter0Dim = parameter0Dim;
	this->parameter1Dim = parameter1Dim;
	this->parameter2Dim = parameter2Dim;
	
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
void Data3D<T>::init(){
	cellLengthParameter0 = abs((*headerParameter0)[1]-(*headerParameter0)[0]);
	cellLengthParameter1 = abs((*headerParameter1)[1]-(*headerParameter1)[0]);
	cellLengthParameter2 = abs((*headerParameter2)[1]-(*headerParameter2)[0]);
	
	minParameter0 = (*headerParameter0)[0];
	minParameter1 = (*headerParameter1)[0];
	minParameter2 = (*headerParameter2)[0];
}

template <class T>
void Data3D<T>::setHeaderParameter0Cell(unsigned int x, float value){
	(*headerParameter0)[x] = value;
}

template <class T>
void Data3D<T>::setHeaderParameter1Cell(unsigned int x, float value){
	(*headerParameter1)[x] = value;
}

template <class T>
void Data3D<T>::setHeaderParameter2Cell(unsigned int x, float value){
	(*headerParameter2)[x] = value;
}

template <class T>
float Data3D<T>::getHeaderParameter0Cell(unsigned int x){
	return (*headerParameter0)[x];
}

template <class T>
float Data3D<T>::getHeaderParameter1Cell(unsigned int x){
	return (*headerParameter1)[x];
}

template <class T>
float Data3D<T>::getHeaderParameter2Cell(unsigned int x){
	return (*headerParameter2)[x];
}

template <class T>
void Data3D<T>::setDataCell(unsigned int x, unsigned int y, unsigned int z, T value){
	(*data)[x][y][z] = value;
}


template <class T>
void Data3D<T>::setAuxiliaryCell(unsigned int x, unsigned int y, unsigned int z, float value){
	(*auxiliary)[x][y][z] = value;
}


template <class T>
void Data3D<T>::setGradientCell(unsigned int x, unsigned int y, unsigned int z, unsigned int i, float value){
	(*gradient)[x][y][z](i) = value;
}

template <class T>
void Data3D<T>::setHessianCell(unsigned int x, unsigned int y, unsigned int z, unsigned int i, unsigned int j, float value){
	(*hessian)[x][y][z](i,j) = value;
}

template <class T>
Vector Data3D<T>::getHeaderVector(unsigned int parameter0, unsigned int parameter1, unsigned int parameter2){
	//printf("HEADER VECTOR: %d %d %d\n",parameter0,parameter1,parameter2);
	Vector v=Vector(3);
	v(0) = (*headerParameter0)[parameter0];
	v(1) = (*headerParameter1)[parameter1];
	v(2) = (*headerParameter2)[parameter2];
	
	return v;
}




template <class T>
float Data3D<T>::parameter2GridLength(){
	return abs((*headerParameter2)[parameter2Dim-1]-(*headerParameter2)[0]);
}


template <class T>
float Data3D<T>::parameter1GridLength(unsigned int parameter2){
	return abs((*headerParameter1)[parameter1Dim-1]-(*headerParameter1)[0]);
}


template <class T>
float Data3D<T>::parameter0GridLength(unsigned int parameter1, unsigned int parameter2){
	return abs((*headerParameter0)[parameter0Dim-1]-(*headerParameter0)[0]);
}


template <class T>
T Data3D<T>::getDataCell(unsigned int x, unsigned int y, unsigned int z){
	return (*data)[x][y][z];
}


template <class T>
float Data3D<T>::getAuxiliaryCell(unsigned int x, unsigned int y, unsigned int z){
	return (*auxiliary)[x][y][z];
}


template <class T>
Vector &Data3D<T>::getGradient(unsigned int x, unsigned int y, unsigned int z){
	return (*gradient)[x][y][z];
}


template <class T>
Matrix &Data3D<T>::getHessian(unsigned int x, unsigned int y, unsigned int z){
	return (*hessian)[x][y][z];
}


template <class T>
bool Data3D<T>::isWithinNumericalLimits(float x, float t){
	if(abs(x-t) <= THRESHOLD_NUMERICAL) return true;
	else return false;
}






template <class T>
void Data3D<T>::closestGridPoint(Vector &x, VectorInt &p, Vector &l){
	unsigned int i_parameter0;
	unsigned int i_parameter1;
	unsigned int i_parameter2;
	Vector d(3);
	
	l(0) = cellLengthParameter0;
	l(1) = cellLengthParameter1;
	l(2) = cellLengthParameter2;
	
	//printf("CELL LENGTHS: %f %f %f\n",l(0),l(1),l(2));
	
	//printf("X: %f %f %f, l: %f %f %f, min: %f %f %f, celll: %f %f %f\n",x(0),x(1),x(2),l(0),l(1),l(2),minParameter0, minParameter1, minParameter2, cellLengthParameter0, cellLengthParameter1, cellLengthParameter2);
	
	d=Vector(3);
	d(0) = floor((x(0)-minParameter0)/cellLengthParameter0);
	d(1) = floor((x(1)-minParameter1)/cellLengthParameter1);
	d(2) = floor((x(2)-minParameter2)/cellLengthParameter2);
	
	i_parameter0 = static_cast<unsigned int>(d(0));
	
	i_parameter1 = static_cast<unsigned int>(d(1));
	
	i_parameter2 = static_cast<unsigned int>(d(2));
	
	//printf("CLOSESTGRRRID: %f %f %f\n",x(2),x(2)-minParameter2,(x(2)-minParameter2)/cellLengthParameter2);
	
	p(0) = i_parameter0;
	p(1) = i_parameter1;
	p(2) = i_parameter2;
}




template <class T>
void Data3D<T>::surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths){
	VectorInt v2(3);
	unsigned int i,j,k;
	VectorInt v(3);
	
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


template <class T>
void Data3D<T>::print(){
	printf("header dim0");
	for(unsigned int i=0; i<parameter0Dim; ++i) printf("\t%f",(*headerParameter0)[i]);
	printf("\nheader dim1");
	for(unsigned int j=0; j<parameter1Dim; ++j) printf("\t%f",(*headerParameter1)[j]);
	printf("\nheader dim2");
	for(unsigned int k=0; k<parameter2Dim; ++k) printf("\t%f",(*headerParameter2)[k]);
}



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
