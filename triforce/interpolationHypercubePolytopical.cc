#include "interpolation.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;



InterpolationHypercubePolytopical::InterpolationHypercubePolytopical(Data6D *data){
	this->data=data;
}


vector<float> InterpolationHypercubePolytopical::weights(vector<VectorInt> &sp, Vector &x, Vector &lengths){
	Vector d=Vector(6);
	vector<float> *r;
	Vector stddist;
	float w;
	float maxw=0;
	r = new vector<float>;
	Vector p;
	
	
	//printf("stddist: %f, %f, %f\n",stddist(0),stddist(1),stddist(2));
	for(unsigned int i=0;i<sp.size();++i){
		p = data->getHeaderVector(sp[i](0),sp[i](1),sp[i](2),sp[i](3),sp[i](4),sp[i](5));
		for(unsigned int j=0;j<6;++j){
			if(i==0 || sp[0](j)==sp[i](j))
				d(j)=(fabs(p(j)-x(j)) / lengths(j));
			else
				d(j)=(fabs(p(j)-x(j)) / (lengths(j)*0.5));
		}
		w = 1.0-max(d(5),max(d(4),max(d(3),max(d(2),max(d(1),d(0))))));
		maxw=maxw+w;
		r->push_back(w);
		
		
	}
	for(unsigned int i=0;i<sp.size();++i){
		r->at(i) = r->at(i)/maxw;
		
		//printf("weight: %f\n", r->at(i));
		
	}
	return *r;
	
}



float InterpolationHypercubePolytopical::interpolate(Vector &x, vector<VectorInt> &sp, vector<float> &w){
	float v;
	Vector lengths(6);
	float t;
	
	v=0;
	
	if(sp.size()==0){
	
		data->surroundingPointsAndCellLengths(x,sp,lengths);
	
	
		if(sp.size()==0){
			printf("NO SPs FOUND. IT's OUTRAGEOUS! (%f, %f, %f)\n",x(0),x(1),x(2));
			//exit(-1);
			return 0;
		}
			
		w = weights(sp,x,lengths);
	}
	
	for(unsigned int i=0;i<sp.size();i++){
		//data->printDataCell(sp[i](0),sp[i](1),sp[i](2));
		//data->printGradientCell(sp[i](0),sp[i](1),sp[i](2));
		//data->printHessianCell(sp[i](0),sp[i](1),sp[i](2));
		t = data0->getDataCell(sp[i](0),sp[i](1),sp[i](2),sp[i](3),sp[i](4),sp[i](5));
		
		//printf("taylor %d: %f [%f]\n",i,t,w[i]);
		v+=w[i]*t;
	}
	
	return v;
	
}


