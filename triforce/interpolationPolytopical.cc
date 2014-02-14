#include "interpolationPolytopical.h"

#include <algorithm>
#include <string>


using namespace std;
using namespace arma;



InterpolationPolytopical::InterpolationPolytopical(Data<float> *data, Data<Vector> *weights){
	this->data=data;
	this->weights=weights;
	dimensions=data->getDimensions();
	dim=dimensions.size();
//	intp=new Interpolation(data,TAYLOR_QUADRATIC);
	daisyChained=false;
	
}

InterpolationPolytopical::InterpolationPolytopical(Data<float> *data, Data<Vector> *weights, Interpolator *daisyChain){
	this->data=data;
	this->weights=weights;
	dimensions=data->getDimensions();
	dim=dimensions.size();
//	intp=new Interpolation(data,TAYLOR_QUADRATIC);
	this->daisyChain=daisyChain;
	daisyChained=true;
	
}

vector<VectorInt> InterpolationPolytopical::getSupportNodes(Vector &d){
	VectorInt p0(3),p(3);
	multimap<float,unsigned int> o;
	multimap<float,unsigned int>::reverse_iterator it;
	vector<VectorInt> sn;
	pair<float, int> a;
	unsigned int k;
	float d0;
	
	sn.clear();
	sn.reserve(dim+1);
	
	p0.zeros();
	//find closest node
	for(unsigned int i=0; i<dim; ++i){
		if(d(i)>0.5){
			p0(i)=1;
			d0 = 1.0-d(i);
		}
		else d0=d(i);
		a.first=d0;
		a.second=i;
		o.insert(a);
	}
	//push back closest node
	sn.push_back(p0);
	
	
	//get neighbour nodes
	for(it=o.rbegin(); it!=o.rend(); ++it){
		p=VectorInt(p0);
		k=it->second;
		if(p(k)==0) p(k)=1;
		else p(k)=0;
		sn.push_back(p);
	}

	return sn;
}

float InterpolationPolytopical::interpolate(Vector &x){
	Vector distsWeights(dim);
	VectorInt g;
	float c;
	float c_weight;
	float w_sum,w_sum_r;
	float e;
	//this gives us the node to the "lower left" in the grid plus normalised distance to the point and the length of the cell
	
	if(daisyChained){
		distsGrid=daisyChain->fetchAuxiliaryFloat(0);
		pGrid=daisyChain->fetchAuxiliaryInt(0);
		pWeights=daisyChain->fetchAuxiliaryInt(1);
	}
	else{
		pGrid=VectorInt(3);
		pWeights=VectorInt(3);
		distsGrid=Vector(3);
		data->closestGridPoint(x, pGrid, distsGrid);
		weights->closestGridPoint(distsGrid, pWeights, distsWeights);
	}

/*	
	printf("CLOSEST GRID POINTS p:(%d %d %d) p1: (%d %d %d)\n",p(0),p(1),p(2),p1(0),p1(1),p1(2));
	printf("d:(%f %f %f)\n",d(0),d(1),d(2));
*/	
	//centroid weight;
	c_weight=1;
	for(unsigned int i=0; i<dim; ++i){
		c_weight=min(c_weight,static_cast<float>(1-fabs(distsGrid(i)-0.5)*2));
	}
	//centroid
	c=data->getAuxiliaryCell(pGrid);
	
	if(daisyChained){
		sp=daisyChain->fetchSupportNodes();
		w=daisyChain->fetchWeights();
	}
	else{
		sp=getSupportNodes(distsGrid);
		w=weights->getDataCell(pWeights);
	}
	
	//printf("CHAIN: %d %d %d\n",daisyChained, this, sp.size());
	
	
	//calculate sum of weights
	w_sum=c_weight;
	for(unsigned int i=0; i<dim; ++i){
		w_sum+=w(i);
	}
	w_sum_r = 1/w_sum;
	

	e=(c_weight*w_sum_r) * c;
	printf("E[-] (- - -): %f\n", c);
	for(unsigned int i=0; i<dim; ++i){
		g = pGrid+sp[i];
		e+= w(i)*w_sum_r * data->getDataCell(g,true);
		printf("E[%d] (%d %d %d): %f\n",i, g(0), g(1), g(2), data->getDataCell(g,true));
	}
	
	
	return e;
	
	
}

Vector &InterpolationPolytopical::fetchAuxiliaryFloat(unsigned int i){
	return distsGrid;
}
VectorInt &InterpolationPolytopical::fetchAuxiliaryInt(unsigned int i){
	switch(i){
		return pGrid;break;
		return pWeights; break;
	}
	return pGrid;
}

