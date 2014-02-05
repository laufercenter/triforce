#include "interpolationPolytopical.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;



InterpolationPolytopical::InterpolationPolytopical(Data<float> *data, Data<Vector> *weights, unsigned int dim){
	this->data=data;
	this->weights=weights;
	this->dim=dim;
}


vector<VectorInt> InterpolationPolytopical::getSupportNodes(Vector &d){
	VectorInt p0(3),p(3);
	multimap<float,unsigned int> o;
	multimap<float,unsigned int>::reverse_iterator it;
	vector<VectorInt> sp;
	Vector d0;
	pair<float, int> a;
	unsigned int k;
	
	
	p0.zeros();
	//find closest node
	for(unsigned int i=0; i<dim; ++i){
		if(d(i)>0.5){
			p0(i)=1;
			a.first=1.0-d(i);
			a.second=i;
			o.insert(a);
		}
	}
	//push back closest node
	sp.push_back(p0);
	
	//get neighbour nodes
	for(it=o.rbegin(); it!=o.rend(); ++it){
		p=p0;
		k=it->second;
		if(p(k)==0) p(k)=1;
		else p(k)=1;
		sp.push_back(p);
	}

	return sp;
}

float InterpolationPolytopical::interpolate(Vector &x){
	VectorInt p(dim),p1(dim);
	Vector d(dim),d1(dim);
	Vector l(dim),l1(dim);
	VectorInt g;
	Vector w;
	float c;
	float c_weight;
	vector<VectorInt> sp;
	float w_sum,w_sum_r;
	float e;
	//this gives us the node to the "lower left" in the grid plus normalised distance to the point and the length of the cell
	data->closestGridPoint(x, p, d, l);
	weights->closestGridPoint(d, p1, d1, l1);
	
	//centroid weight;
	c_weight=1;
	for(unsigned int i=0; i<dim; ++i){
		c_weight=min(c_weight,static_cast<float>(1-fabs(d(i)-0.5)*2));
	}
	//centroid
	c=data->getAuxiliaryCell(p);
	
	sp=getSupportNodes(d);
	w=weights->getDataCell(p1);
	
	//calculate sum of weights
	w_sum=c_weight;
	for(unsigned int i=0; i<dim; ++i){
		w_sum+=w(i);
	}
	w_sum_r = 1/w_sum;
	

	e=(c_weight*w_sum_r) * c;
	for(unsigned int i=0; i<dim; ++i){
		g = p+sp[i];
		e+= w(i)*w_sum_r * data->getDataCell(g);
	}
	
	return e;
	
	
}

