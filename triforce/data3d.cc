#include "data3d.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;




Data3D::Data3D(int PHIDim0, int PHIDim1, int PHIDim2, int psiDim0, int psiDim1, int lambdaDim){
	int maxdim=0;
	this->PHIDim0 = PHIDim0;
	this->PHIDim1 = PHIDim1;
	this->PHIDim2 = PHIDim2;
	this->psiDim0 = psiDim0;
	this->psiDim1 = psiDim1;
	this->lambdaDim = lambdaDim;
	
	
		
		
	headerPHI = new Table3dDouble(boost::extents[PHIDim0][PHIDim1][PHIDim2]);
	headerPsi = new Table2dDouble(boost::extents[psiDim0][psiDim1]);
	headerLambda = new Table1dDouble(boost::extents[lambdaDim]);
	data = new Table3dDouble(boost::extents[PHIDim0][psiDim0][lambdaDim]);
	gradient = new Table3dVector(boost::extents[PHIDim0][psiDim0][lambdaDim]);
	hessian = new Table3dMatrix(boost::extents[PHIDim0][psiDim0][lambdaDim]);
	for(int x=0; x<PHIDim0; x++)
		for(int y=0; y<psiDim0; y++)
			for(int z=0; z<lambdaDim; z++){
				boost::array<Table3dVector::index,3> idx = {{x,y,z}};
				(*gradient)(idx) = Vector(3);
				(*hessian)(idx) = Matrix(3,3);
			}
}


void Data3D::setHeaderPHICell(int x, int y, int z, double value){
	(*headerPHI)[x][y][z] = value;
}

void Data3D::setHeaderPsiCell(int y, int z, double value){
	(*headerPsi)[y][z] = value;
}

void Data3D::setHeaderLambdaCell(int z, double value){
	(*headerLambda)[z] = value;
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
	Vector v=Vector(3);
	v(0) = (*headerPHI)[PHI][psi][lambda];
	v(1) = (*headerPsi)[psi][lambda];
	v(2) = (*headerLambda)[lambda];
	
	return v;
}

double Data3D::lambdaCellLength(){
	return abs((*headerLambda)[1]-(*headerLambda)[0]);
}

double Data3D::psiCellLength(int lambda){
	return abs((*headerPsi)[1][lambda]-(*headerPsi)[0][lambda]);
}

double Data3D::PHICellLength(int psi, int lambda){
	return abs((*headerPHI)[1][psi][lambda]-(*headerPHI)[0][psi][lambda]);
}

double Data3D::lambdaGridLength(){
	return abs((*headerLambda)[lambdaDim-1]-(*headerLambda)[0]);
}

double Data3D::psiGridLength(int lambda){
	return abs((*headerPsi)[psiDim0-1][lambda]-(*headerPsi)[0][lambda]);
}

double Data3D::PHIGridLength(int psi, int lambda){
	return abs((*headerPHI)[PHIDim0-1][psi][lambda]-(*headerPHI)[0][psi][lambda]);
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
void Data3D::lowestGridPointAndCellLengths(Vector &x, VectorInt &p, Vector &lengths){
	double lengthLambda, lengthPsi, lengthPHI;
	int i_lambda, i_psi, i_PHI;
	
	//first, lambda:
	lengthLambda = lambdaCellLength();
	i_lambda = floor(x(2)/(double)lengthLambda);
	
	lengthPsi = psiCellLength(i_lambda);
	i_psi = floor(x(1)/(double)lengthPsi);
	
	lengthPHI = PHICellLength(i_psi,i_lambda);
	i_PHI = floor(x(1)/(double)lengthPHI);
	
	p(0) = i_PHI;
	p(1) = i_psi;
	p(2) = i_lambda;
	
	lengths(0) = lengthPHI;
	lengths(1) = lengthPsi;
	lengths(2) = lengthLambda;
	
}


void Data3D::surroundingPointsandCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths){
	VectorInt v;
	VectorInt *v2;
	int i,j,k;
	
	
	lowestGridPointAndCellLengths(x,v,lengths);
	r.clear();
	for(i=0;i<2;++i)
		for(j=0;j<2;++j)
			for(k=0;k<2;++k){
				v2=new VectorInt(3);
				(*v2)(0)=v(0)+i;
				(*v2)(1)=v(1)+j;
				(*v2)(2)=v(2)+k;
				r.push_back(*v2);
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
