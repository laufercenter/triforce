#include "interpolation.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;



Interpolation::Interpolation(Data3D *data, TaylorTermination degree){
	this->data=data;
	this->degree=degree;
}

float Interpolation::taylorExtension(VectorInt &r, Vector &x){
	return taylorExtension(r(0),r(1),r(2),x);
}

float Interpolation::taylorExtension(int i_PHI, int i_psi, int i_lambda, Vector &x){
	Vector p;
	float f;
	Vector g;
	Matrix h;
	float v;
	Vector d;
	Vector htmp;
	bool quadratic=false;
	bool cubic=false;
	float quadratic_term=0;
	float cubic_term=0;
	if(degree==TAYLOR_QUADRATIC || degree==TAYLOR_CUBIC) quadratic=true;
	if(degree==TAYLOR_CUBIC) cubic=true;
	
	p = data->getHeaderVector(i_PHI,i_psi, i_lambda);
	f = data->getDataCell(i_PHI,i_psi, i_lambda);
	
	if(quadratic)
		g = data->getGradient(i_PHI,i_psi, i_lambda);
	if(cubic)
		h = data->getHessian(i_PHI,i_psi, i_lambda);
	
	if(quadratic){
		d =(x-p);
		quadratic_term = dot(g,d);
	}
	
	if(cubic){
		htmp = h * d;
		cubic_term = 0.5*dot(d,htmp);
	}
	v = f + quadratic_term + cubic_term;
	
	//printf("taylor extension (%d %d %d): %f\n",i_PHI,i_psi,i_lambda, v);
	
	//printf("taylor extension (%d %d %d): d:%f grad:(%f %f %f) v:%f\n",i_PHI,i_psi,i_lambda,f,g(0),g(1),g(2), v);
	
	
	//v = dataConvex[i_PHI, i_psi, i_lambda] + c(t(gradientsConvex[,i_PHI, i_psi, i_lambda]) %*% (x-p)) + 0.5 *(c(t((x-p)) %*% hessiansConvex[,,i_PHI, i_psi, i_lambda] %*% (x-p)))
	return v;
}

vector<float> Interpolation::weights(vector<VectorInt> &sp, Vector &x, Vector &lengths){
	Vector d=Vector(3);
	vector<float> *r;
	Vector stddist;
	float w;
	float maxw=0;
	r = new vector<float>;
	Vector p;
	
	//printf("stddist: %f, %f, %f\n",stddist(0),stddist(1),stddist(2));
	for(unsigned int i=0;i<sp.size();++i){
		p = data->getHeaderVector(sp[i](0),sp[i](1),sp[i](2));
		for(unsigned int j=0;j<3;++j)
			d(j)=(fabs(p(j)-x(j)) / lengths(j));
		w = 1.0-max(d(0),max(d(1),d(2)));
		maxw=maxw+w;
		r->push_back(w);
		
		
	}
	for(unsigned int i=0;i<sp.size();++i){
		r->at(i) = r->at(i)/maxw;
		
		//printf("weight: %f\n", r->at(i));
		
	}
	return *r;
	
}



float Interpolation::interpolate(Vector &x){
	return multiPointTaylor(x);
}




float Interpolation::multiPointTaylor(Vector &x){
	vector<VectorInt> sp;
	vector<float> w;
	float v=0;
	Vector lengths(3);
	float closestSPWeight;
	int closestSP;
	
	
	data->surroundingPointsAndCellLengths(x,sp,lengths);
	
	
	if(sp.size()==0){
		printf("NO SPs FOUND. IT's OUTRAGEOUS! (%f, %f, %f)\n",x(0),x(1),x(2));
		//exit(-1);
		return 0;
	}
		
	w = weights(sp,x,lengths);
	closestSP = -1;
	closestSPWeight = 0;
	for(unsigned int i=0;i<sp.size();i++){
		//data->printDataCell(sp[i](0),sp[i](1),sp[i](2));
		//data->printGradientCell(sp[i](0),sp[i](1),sp[i](2));
		//data->printHessianCell(sp[i](0),sp[i](1),sp[i](2));
		
		float t = taylorExtension(sp[i],x);
		
		//printf("taylor %d: %f [%f]\n",i,t,w[i]);
		v+=w[i]*t;
		
		if(w[i] >= closestSPWeight)
			closestSP = i;
	}
	
	//phi = data->getAuxiliaryCell(sp[closestSP](0),sp[closestSP](1),sp[closestSP](2));
	
	
	/*
	for(int i=0;i<sp.size();i++){
		float t = taylorExtension(sp[i],x);
		
		if(abs(t-v) > 0.1){
			printf("THIS INTERPOLATION IS OUTRAGEOUS: %f\n",abs(t-v));
			//exit(-1);
		}
	}
	*/
	
	//printf("---- %f\n",v);
	
	return v;
}
