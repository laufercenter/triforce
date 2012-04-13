#include "spline.h"


typedef mat Matrix;

Spline::Spline(vector<Vector> geometry){
	int i,l;
	
	Vector v,g;
	

	//first define constant cutmull-rom base matrix
	//double c[16]={-1,2,-1,0,3,-5,0,2,-3,4,1,0,1,-1,0,0};
	double c[16]={0,2,0,0,-1,0,1,0,2,-5,4,-1,-1,3,-3,1};
	B = Matrix(c,4,4) * 0.5;

	l=geometry.size();
	
	//store x values in separate vector
	for(i=0; i<l; i++){
		X.push_back(geometry[i](0));
	}
	
	

	//padding of geometry vector
	v=Vector(2);
	v(0)=geometry[0][(0) - 1;
	v(1)=geometry[0][(1);
	geometry.insert(geometry.begin(),v);
	v(0)=geometry[geometry.size()-1][(0) + 1;
	v(1)=geometry[geometry.size()-1][(1);
	geometry.insert(geometry.end(),v);
	
	//creating geometry matrices
	for(i=0; i<l; i++){
		g=Vector(4);
		g(0) = geometry[i](1);
		g(1) = geometry[i+1](1);
		g(2) = geometry[i+2](1);
		g(3) = geometry[i+3](1);
		
		G.push_back(B * g);
		
	}
	
	

}


int Spline::logSearch(double x){
	int l,r,c;
	l=0;
	r=X.size()-1;
	c=(l-r)/2
	while(c!=l){
		if(X[c]>x) r=c;
		else l=c;
	}
	return l;
}

double Spline::f(double x){
	int c;
	double x0,x1,t;
	Vector T;
	
	//determine bin in which x resides
	c = logSearch(x);
	
	//determine where it is in the bin [0..1]
	x0 = X[c];
	x1 = X[c+1];
	t = (c-x0) / (x1-x0);
	
	//calculate
	t_square= 
	t_cubic = t_square * t;
	T = Vector(4);
	T(0)=1;
	T(1)=t;
	T(2)=t_square;
	T(3)=t_cubic;
	
	res = T * G[c];
	
	return res;
	
		
}


