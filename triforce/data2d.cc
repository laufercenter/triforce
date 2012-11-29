#include "data2d.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;




Data2D::Data2D(vector<int> &dimensions){
	int maxdim=0;
	this->dimensions.insert(this->dimensions.begin(),dimensions.begin(),dimensions.end());
	
	
	for(int i=0; i<dimensions.size(); i++)
		if(dimensions[i]>maxdim) maxdim=dimensions[i];
		
		
	header = new Table2dDouble(boost::extents[2][maxdim]);
	data = new Table2dDouble(boost::extents[dimensions[0]][dimensions[1]]);
}


void Data2D::setHeaderCell(int row, int col, double value){
	(*header)[row][col] = value;
}

void Data2D::setDataCell(int x, int y, double value){
	(*data)[x][y] = value;
}


Vector Data2D::getHeaderVector(int x, int y){
	Vector v=Vector(2);
	v(0) = (*header)[0][x];
	v(1) = (*header)[1][y];
	
	return v;
}

double Data2D::getDataCell(int x, int y){
	return (*data)[x][y];
}


Vector Data2D::bisectFloor(Vector &x){
	int l,r,m;
	int s;
	double vl,vr;
	Vector v=Vector(2);
	for(int i=0;i<2;i++){
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


vector<Vector> Data2D::surroundingPoints(Vector &x){
	Vector v;
	vector<Vector> r;
	Vector *v2;
	int i,j,k;
	
	
	v=bisectFloor(x);
	
	for(i=0;i<2;++i)
		for(j=0;j<2;++j){
			v2=new Vector(2);
			(*v2)(0)=v(0)+i;
			(*v2)(1)=v(1)+j;
			r.push_back(*v2);
		}
			
	return r;
}

Vector Data2D::standardDistance(){
	Vector r = Vector(2);
	for(int i=0;i<2;i++)
		r(i) = fabs((*header)[i][1]-(*header)[i][0]);
	
	return r;
}

void Data2D::print(){
	for(int y=0; y<dimensions[0]; y++){
		for(int x=0; x<dimensions[1]; x++){
			double v=(*data)[x][y];
			printf("%f ",v);
		}
		printf("\n");
	}
}

void Data2D::printDataCell(int i, int j){
	printf("Cell[%d,%d]: %f\n",i,j,(*data)[i][j]);
}

