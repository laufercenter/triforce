#include "data3d.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;
using namespace boost;




Data3D::Data3D(vector<int> &dimensions){
	int maxdim=0;
	this->dimensions.insert(this->dimensions.begin(),dimensions.begin(),dimensions.end());
	
	
	for(int i=0; i<dimensions.size(); i++)
		if(dimensions[i]>maxdim) maxdim=dimensions[i];
		
	header = new Table2dDouble(boost::extents[3][maxdim]);
	data = new Table3dDouble(boost::extents[dimensions[0]][dimensions[1]][dimensions[2]]);
	gradient = new Table3dVector(boost::extents[dimensions[0]][dimensions[1]][dimensions[2]]);
	hessian = new Table3dMatrix(boost::extents[dimensions[0]][dimensions[1]][dimensions[2]]);
	for(int x=0; x<dimensions[0]; x++)
		for(int y=0; y<dimensions[1]; y++)
			for(int z=0; z<dimensions[2]; z++){
				boost::array<Table3dVector::index,3> idx = {{x,y,z}};
				(*gradient)(idx) = Vector(3);
				(*hessian)(idx) = Matrix(3,3);
			}
}


void Data3D::setHeaderCell(int row, int col, double value){
	(*header)[row][col] = value;
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

void Data3D::print(){
	for(int z=0; z<dimensions[0]; z++){
		for(int x=0; x<dimensions[1]; x++){
			for(int y=0; y<dimensions[2]; y++){
				double v=(*data)[x][y][z];
				printf("%f ",v);
			}
			printf("\n");
		}
		printf("\n\n");
	}
}




