#include "data3d.h"

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;

#define THRESHOLD_NUMERICAL 0.00001



Data3D::Data3D(int PHIDim, int psiDim, int lambdaDim){
	int maxdim=0;
	this->PHIDim = PHIDim;
	this->psiDim = psiDim;
	this->lambdaDim = lambdaDim;
	
	
		
		
	headerPHI = new Table1dDouble(boost::extents[PHIDim]);
	headerPsi = new Table1dDouble(boost::extents[psiDim]);
	headerLambda = new Table1dDouble(boost::extents[lambdaDim]);
	data = new Table3dDouble(boost::extents[PHIDim][psiDim][lambdaDim]);
	gradient = new Table3dVector(boost::extents[PHIDim][psiDim][lambdaDim]);
	hessian = new Table3dMatrix(boost::extents[PHIDim][psiDim][lambdaDim]);
	for(int x=0; x<PHIDim; x++)
		for(int y=0; y<psiDim; y++)
			for(int z=0; z<lambdaDim; z++){
				boost::array<Table3dVector::index,3> idx = {{x,y,z}};
				(*gradient)(idx) = Vector(3);
				(*hessian)(idx) = Matrix(3,3);
			}
}


void Data3D::setHeaderPHICell(int x, double value){
	(*headerPHI)[x] = value;
}

void Data3D::setHeaderPsiCell(int x, double value){
	(*headerPsi)[x] = value;
}

void Data3D::setHeaderLambdaCell(int x, double value){
	(*headerLambda)[x] = value;
}


void Data3D::setDataCell(int x, int y, int z, double value){
	(*data)[x][y][z] = value;
}

void Data3D::setGradientCell(int x, int y, int z, int i, double value){
	(*gradient)[x][y][z](i) = value;
}
void Data3D::setHessianCell(int x, int y, int z, int i, int j, double value){
	(*hessian)[x][y][z](i,j) = value;
}

Vector Data3D::getHeaderVector(int PHI, int psi, int lambda){
	//printf("HEADER VECTOR: %d %d %d\n",PHI,psi,lambda);
	Vector v=Vector(3);
	v(0) = (*headerPHI)[PHI];
	v(1) = (*headerPsi)[psi];
	v(2) = (*headerLambda)[lambda];
	
	return v;
}

Vector Data3D::cellLength(){
	Vector v(3);
	
	v(0) = abs((*headerPHI)[1]-(*headerPHI)[0]);
	v(1) = abs((*headerPsi)[1]-(*headerPsi)[0]);
	v(2) = abs((*headerLambda)[1]-(*headerLambda)[0]);
	
	return v;
}


double Data3D::lambdaGridLength(){
	return abs((*headerLambda)[lambdaDim-1]-(*headerLambda)[0]);
}

double Data3D::psiGridLength(int lambda){
	return abs((*headerPsi)[psiDim-1]-(*headerPsi)[0]);
}

double Data3D::PHIGridLength(int psi, int lambda){
	return abs((*headerPHI)[PHIDim-1]-(*headerPHI)[0]);
}



double Data3D::getDataCell(int x, int y, int z){
	return (*data)[x][y][z];
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



/*
Vector Data3D::bisectFloor(Vector &x){
	int l,r,m;
	int s;
	double vl,vr;
	Vector v=Vector(3);
	for(int i=0;i<3;i++){
		l=0;
		vl=0;
		r=dimensions[i]-1;
		vr=numeric_limits<double>::max();
		//printf("start bisection(%d): %d %d (%f)\n",i,l,r,x(i));
		s=0;
		while(r-l > 1){
			m = l+((r-l)/2);
			if((*header)[i][m]<=x(i)){
				l=m;
				vl=(*header)[i][m];
			}
			if((*header)[i][m]>x(i)){
				r=m;
				vr=(*header)[i][m];
			}
			//printf("step: l%d m%d(%f) r%d \n",l,m,(*header)[i][m],r);
			++s;
			if(s>20) exit(-1);
		}
		v(i)=l;		
	}
	return v;
}
*/





void Data3D::closestGridPoint(Vector &x, VectorInt &p, Vector &l){
	double lengthPHI;
	int i_PHI;
	double lengthPsi;
	int i_psi;
	double lengthLambda;
	int i_lambda;
	
	l = cellLength();
	
	//printf("CELL LENGTHS: %f %f %f\n",l(0),l(1),l(2));
	
	lengthPHI = l(0);
	i_PHI = static_cast<int>(x(0)/lengthPHI);
	
	lengthPsi = l(1);
	i_psi = static_cast<int>(x(1)/lengthPsi);
	
	lengthLambda = l(2);
	i_lambda = static_cast<int>(x(2)/lengthLambda);
	
	p(0) = i_PHI;
	p(1) = i_psi;
	p(2) = i_lambda;
}




void Data3D::surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths){
	VectorInt v2(3);
	int i,j,k;
	VectorInt v(3);
	bool neg;
	
	closestGridPoint(x, v, lengths);
	
	if((*headerPsi)[v(1)]+(*headerLambda)[v(2)] < M_PI && !isWithinNumericalLimits((*headerPsi)[v(1)]+(*headerLambda)[v(2)],M_PI)) neg=true;
	else neg=false;
	
	
	//printf("CLOSEST GRIDPOINT: %d %d %d\n",v(0),v(1),v(2));
	
	
	r.clear();
	for(i=0;i<2;++i)
		for(j=0;j<2;++j)
			for(k=0;k<2;++k){
				v2(0)=v(0)+i;
				v2(1)=v(1)+j;
				v2(2)=v(2)+k;
				if(v2(0)<PHIDim && v2(1)<psiDim && v2(2)<lambdaDim){
					if(!isnan((*data)[v2(0)][v2(1)][v2(2)])){
						if(!neg || ((*headerPsi)[v2(1)]+(*headerLambda)[v2(2)]<M_PI && !isWithinNumericalLimits((*headerPsi)[v2(1)]+(*headerLambda)[v2(2)],M_PI))){
							r.push_back(v2);
							//printf("ACCEPTED %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerPsi)[v2(1)],(*headerLambda)[v2(2)],neg,isWithinNumericalLimits((*headerPsi)[v2(1)]+(*headerLambda)[v2(2)],M_PI));
						}
						//else	printf("REJECTED 0 %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerPsi)[v2(1)],(*headerLambda)[v2(2)],neg,isWithinNumericalLimits((*headerPsi)[v2(1)]+(*headerLambda)[v2(2)],M_PI));

					}
					//else	printf("REJECTED 1 %d %d %d\n",v2(0),v2(1),v2(2));

				}
				//else	printf("PRE-REJECTED 2 %d %d %d\n",v2(0),v2(1),v2(2));
			}

			
}

void Data3D::print(){
	/*
	printf("DATA...\n");
	for(int z=0; z<dimensions[2]; z++){
		for(int y=0; y<dimensions[1]; y++){
			for(int x=0; x<dimensions[0]; x++){
				double v=(*data)[x][y][z];
				printf("%f ",v);
			}
			printf("\n");
		}
		printf("\n\n");
	}
	
	printf("GRADIENTS...\n");
	for(int z=0; z<dimensions[2]; z++){
		for(int y=0; y<dimensions[1]; y++){
			for(int x=0; x<dimensions[0]; x++){
				Vector v=(*gradient)[x][y][z];
				printf("[%f %f %f] ",v(0),v(1),v(2));
			}
			printf("\n");
		}
		printf("\n\n");
	}
	*/
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
