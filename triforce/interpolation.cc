/**
	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
	Email: nils.drechsel@gmail.com
	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
*/

#include "interpolation.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;



Interpolation::Interpolation(Data<float> *data, TaylorTermination degree){
	this->data=data;
	this->degree=degree;
	dimensions=data->getDimensions();
	dim=dimensions.size();
	lengths=data->getReciprocalCellLengths();
	createTemplates();
	daisyChained=false;
}


Interpolation::Interpolation(Data<float> *data, TaylorTermination degree, Interpolator *daisyChain){
	this->data=data;
	this->degree=degree;
	dimensions=data->getDimensions();
	dim=dimensions.size();
	lengths=data->getReciprocalCellLengths();
	createTemplates();
	
	this->daisyChain=daisyChain;
	daisyChained=true;
}


void Interpolation::createTemplates(){
	float k;
	VectorInt v(dim);
	
	templates.clear();
	k=2<<dim;
	for(unsigned int i=0; i<k; ++i){
		v.zeros();
		for(unsigned int j=0; j<dim; ++j){
			if((i & (1<<j)) > 0) v(j)=1;
		}
		templates.push_back(v);
	}
}


float Interpolation::taylorExtension(VectorInt &r, Vector &x){
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
	
	p = data->getHeaderVector(r);
	f = data->getDataCell(r);
	
	if(quadratic)
		g = data->getGradient(r);
	if(cubic)
		h = data->getHessian(r);
	
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

Vector Interpolation::weights(vector<VectorInt> &sp, Vector &x){
	Vector d(3);
	Vector r(sp.size());
	Vector stddist;
	float wght;
	float maxw=0;
	float maxw_r;
	Vector p;
	
	//printf("stddist: %f, %f, %f\n",stddist(0),stddist(1),stddist(2));
	for(unsigned int i=0;i<sp.size();++i){
		p = data->getHeaderVector(sp[i]);
		for(unsigned int j=0;j<3;++j)
			d(j)=(fabs(p(j)-x(j)) * lengths(j));
		wght = 1.0-max(d(0),max(d(1),d(2)));
		maxw=maxw+wght;
		r(i)=wght;
	}
	maxw_r=1.0/maxw;
	for(unsigned int i=0;i<sp.size();++i){
		r(i) = r(i)*maxw_r;
		
		//printf("weight: %f\n", r->at(i));
		
	}
	return r;
	
}



float Interpolation::interpolate(Vector &x){
	return multiPointTaylor(x);
}




float Interpolation::multiPointTaylor(Vector &x){
	VectorInt v2;
	VectorInt v;
	bool error;
	Vector d(dim);
	float e;
	
	//printf("INTERPOLATING (%f, %f, %f)\n",x(0),x(1),x(2));
	
	data->closestGridPoint(x, v, d);
	

	if(daisyChained){
		sp=fetchSupportNodes();
	}
	else{
		sp.clear();
		sp.reserve(templates.size());
		for(unsigned int i=0; i<templates.size(); ++i){
			v2=v+templates[i];
			error=false;
			for(unsigned int j=0; j<dim && !error; ++j){
				if(v2(j)>=dimensions(j)) error=true;
			}
			if(!error) sp.push_back(v+templates[i]);
		}
	}
	
	
	if(sp.size()==0){
		printf("NO SPs FOUND. IT's OUTRAGEOUS! (%f, %f, %f)\n",x(0),x(1),x(2));
		//exit(-1);
		return 0;
	}
	
	if(daisyChained){
		w = daisyChain->fetchWeights();
		
	}
	else{
		w = weights(sp,x);
	}
	
	
	e=0;
	for(unsigned int i=0;i<sp.size();i++){
		//data->printDataCell(sp[i](0),sp[i](1),sp[i](2));
		//data->printGradientCell(sp[i](0),sp[i](1),sp[i](2));
		//data->printHessianCell(sp[i](0),sp[i](1),sp[i](2));
		
		float t = taylorExtension(sp[i],x);
		
		//printf("taylor %d: %f [%f]\n",i,t,w[i]);
		e+=w[i]*t;
		
		//printf("E2[%d] (%d %d %d): %f\n",i, sp[i](0), sp[i](1), sp[i](2), t);
		
		
	}
	
	
	
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
	
	return e;
}




