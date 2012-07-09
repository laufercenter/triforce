#include "data3d.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;




Data3D::Data3D(vector<int> &dimensions){
	int maxdim=0;
	this->dimensions.insert(this->dimensions.begin(),dimensions.begin(),dimensions.end());
	
	
	for(int i=0; i<dimensions.size(); i++)
		if(dimensions[i]>maxdim) maxdim=dimensions[i];
		
	header = Table2dDouble(boost::extents[3][maxdim]);
	data = Table3dDouble(boost::extents[dimensions[0]][dimensions[1]][dimensions[2]]);
	gradient = Table3dVector(boost::extents[dimensions[0]][dimensions[1]][dimensions[2]]);
	hessian = Table3dMatrix(boost::extents[dimensions[0]][dimensions[1]][dimensions[2]]);
	for(int x=0; x<dimensions[0]; x++)
		for(int y=0; y<dimensions[1]; y++)
			for(int z=0; z<dimensions[2]; z++){
				boost::array<Table3dVector::index,3> idx = {{x,y,z}};
				gradient(idx) = Vector(3);
				hessian(idx) = Matrix(3,3);
			}
}


void Data3D::setHeaderCell(int row, int col, double value){
	header[row][col] = value;
}

void Data3D::setDataCell(int x, int y, int z, double value){
	data[x][y][z] = value;
}

void Data3D::setGradientCell(int x, int y, int z, int i, double value){
	gradient[x][y][z](i) = value;
}
void Data3D::setHessianCell(int x, int y, int z, int i, int j, double value){
	hessian[x][y][z](i,j) = value;
}

Vector Data3D::getHeaderVector(int x, int y, int z){
	Vector v=Vector(3);
	v(0) = header[0][x];
	v(1) = header[1][y];
	v(2) = header[2][z];
	
	return v;
}

double Data3D::getDataCell(int x, int y, int z){
	return data[x][y][z];
}

Vector &Data3D::getGradient(int x, int y, int z){
	return gradient[x][y][z];
}

Matrix &Data3D::getHessian(int x, int y, int z){
	return hessian[x][y][z];
}

Vector Data3D::bisectFloor(Vector &x){
	int l,r,m;
	double vl,vr;
	Vector v=Vector(3);
	for(int i=0;i<3;i++){
		l=0;
		vl=0;
		r=dimensions[i]-1;
		vr=numeric_limits<double>::max();
		while(r-l > 1){
			m = (r-l)/2;
			if(header[i][m]<=x(i)){
				l=m;
				vl=header[i][m];
			}
			else if(header[i][m]>x(i)){
				r=m;
				vr=header[i][m];
			}
		}
		v(i)=l;		
	}
	return v;
}


vector<Vector> Data3D::surroundingPoints(Vector &x){
	Vector v;
	vector<Vector> r;
	Vector *v2;
	int i,j,k;
	
	
	v=bisectFloor(x);
	
	for(i=0;i<2;++i)
		for(j=0;j<2;++j)
			for(k=0;k<2;++k){
				v2=new Vector(3);
				(*v2)(0)=v(0)+i;
				(*v2)(1)=v(1)+j;
				(*v2)(2)=v(2)+k;
				r.push_back(*v2);
			}
			
	return r;
}

Vector Data3D::standardDistance(){
	Vector r = Vector(3);
	for(int i=0;i<3;i++)
		r(i) = abs(header[i][1]-header[i][0]);
	
	return r;
}

void Data3D::print(){
	for(int z=0; z<dimensions[0]; z++){
		for(int x=0; x<dimensions[1]; x++){
			for(int y=0; y<dimensions[2]; y++){
				double v=data[x][y][z];
				printf("%f ",v);
			}
			printf("\n");
		}
		printf("\n\n");
	}
}




