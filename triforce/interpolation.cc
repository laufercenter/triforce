#include "interpolation.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;


Interpolation::Interpolation(Data3D *data){
	this->data=data;
}

double Interpolation::taylorExtension(Vector &r, Vector &x){
	return taylorExtension(r(0),r(1),r(2),x);
}

double Interpolation::taylorExtension(int i_PHI, int i_psi, int i_lambda, Vector &x){
	Vector p;
	double f;
	Vector g;
	Vector h;
	double v;
	Vector d;
	Vector htmp;
	p = data->getHeaderVector(i_PHI,i_psi, i_lambda);
	f = data->getDataCell(i_PHI,i_psi, i_lambda);
	g = data->getGradient(i_PHI,i_psi, i_lambda);
	h = data->getHessian(i_PHI,i_psi, i_lambda);
	
	d =(x-p);
	htmp = h * d;
	v = f + dot(g,d) + 0.5*dot(d,htmp);
	//v = dataConvex[i_PHI, i_psi, i_lambda] + c(t(gradientsConvex[,i_PHI, i_psi, i_lambda]) %*% (x-p)) + 0.5 *(c(t((x-p)) %*% hessiansConvex[,,i_PHI, i_psi, i_lambda] %*% (x-p)))
	return v;
}

vector<double> Interpolation::weights(vector<Vector> &sp, Vector &x){
	Vector d=Vector(3);
	vector<double> *r;
	Vector stddist;
	double w;
	double maxw=0;
	r = new vector<double>;
	Vector p;
	stddist =data->standardDistance();
	for(int i=0;i<sp.size();++i){
		p = data->getHeaderVector(sp[i](0),sp[i](1),sp[i](2));
		for(int j=0;j<3;++j)
			d(j)=(abs(p(j)-x(j)) / stddist(j));
		w = max(d(0),max(d(1),d(2)));
		maxw=maxw+w;
		r->push_back(w);
		
	}
	for(int i=0;i<sp.size();++i){
		r->at(i) = r->at(i)/maxw;
	}
	return *r;
	
}




double Interpolation::multiPointTaylor(Vector &x){
	vector<Vector> sp;
	vector<double> w;
	double v=0;
	
	sp = data->surroundingPoints(x);
	w = weights(sp,x);
	for(int i=0;i<sp.size();i++)
		v+=w[i]*taylorExtension(sp[i],x);
	
	return v;
}
