#include "integrator.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;


Integrator::Integrator(){
	
}



Integrator::Integrator(Tessellation *tessellation, Interpolation *interpolation){
	this->tessellation = tessellation;
	this->interpolation = interpolation;
	
}

Vector Integrator::optimalIntegrationOrigin(list<IntersectionPoint*>* sasa){
	Vector v(3);
	v(0) = 1;
	v(1) = 0;
	v(2) = 0;
	
}

double Integrator::Integrate(){
	vector<vector<CircularRegion>* >* circularRegions;
	vector<vector<list<IntersectionPoint*>*>*>* intersections;
	vector<list<IntersectionPoint*>*>* sasas;
	list<IntersectionPoint*>* sasa;
	vector<CircularRegion>* circles;
	Vector integrationOrigin;
	double area=0;
	
	intersections = tessellation->intersectionPoints();
	circularRegions = tessellation->circularRegions();
	
	//iterate over all atoms
	for(int i=0;i<intersections->size();++i){
		//iterate over all sasas
		sasas = intersections->at(i);
		circles = circularRegions->at(i);
		for(int j=0;j<=sasas->size();++j){
			sasa = sasas->at(j);
			integrationOrigin = optimalIntegrationOrigin(sasa);
			area += integrateSASA(*sasa,*circles, integrationOrigin);
		}
		
	}
	
	return area;
	
	
}

double Integrator::integrateSASA(list<IntersectionPoint*> &sasa, vector<CircularRegion> &circles, Vector &integrationOrigin){
	CircularRegion c;
	list<IntersectionPoint*>* frontHemisphere;
	list<IntersectionPoint*>* backHemisphere;
	int ci;
	double area=0;
	
	c.vector = integrationOrigin;
	c.openingAngle = M_PI/2;
	c.form=CONVEX;
	circles.push_back(c);
	ci = circles.size()-1;
	
	splitSASA(sasa, circles, ci, integrationOrigin, &frontHemisphere, &backHemisphere);
	
	area += integrateHemisphere(*frontHemisphere, integrationOrigin, circles);	
	area += integrateHemisphere(*backHemisphere, integrationOrigin, circles);
	
	delete frontHemisphere;
	delete backHemisphere;
	
	return area;
	
}


double Integrator::integrateHemisphere(list<IntersectionPoint*> &sasa, Vector &integrationOrigin, vector<CircularRegion> &circles){
	list<IntersectionPoint*>::iterator it;
	IntersectionPoint *x0, *x1;
	double area=0;
	
	x0 = *(--sasa.end());
	for(it = sasa.begin(); it!=sasa.end(); ++it){
		x1 = *it;
		
		area += integrateTriangle(*x0, *x1, integrationOrigin, circles);
		
	}
	
	return area;

}

double Integrator::integrateTriangle(IntersectionPoint &x0, IntersectionPoint &x1, Vector integrationOrigin, vector<CircularRegion> &circles){
	double psi;
	double lambda;
	double PHI0;
	double PHI1;
	CircularRegion c;
	c = circles[x0.with];
	double v0,v1;
	double area = 0;
	Vector vec0,vec1;
	Vector n, o;
	double maxArea;
	
	//calculate psi
	psi = angle(c.vector, integrationOrigin);
	//calculate lambda
	lambda = c.openingAngle;
	
	n = cross(c.vector, integrationOrigin);
	o = cross(n,c.vector);
	
	//calculate PHI0
	vec0 = x0.vector - c.vector;
	PHI0 = complLongAngle(n, o, vec0);

	//calculate PHI1
	vec1 = x1.vector - c.vector;
	PHI1 = complLongAngle(n, o, vec1);
	
	
	maxArea = interpolation->interpolate(0, psi, lambda);
	
	if(PHI0 >= 0){
		if(PHI1 >= 0){
			if(PHI0 <= PHI1){
				area += interpolation->interpolate(PHI0, psi, lambda);
				area -= interpolation->interpolate(PHI1, psi, lambda);
			}
			else{
				area += interpolation->interpolate(PHI0, psi, lambda);
				area += maxArea;
				area += maxArea;
				area -= interpolation->interpolate(PHI1, psi, lambda);
			}
		}
		else{
			area += interpolation->interpolate(PHI0, psi, lambda);
			area += interpolation->interpolate(PHI1, psi, lambda);
		}
	}
	else{
		if(PHI1 >= 0){
			area += maxArea;
			area -= interpolation->interpolate(PHI0, psi, lambda);
			area += maxArea;
			area -= interpolation->interpolate(PHI1, psi, lambda);
			
		}
		else{
			if(PHI0 >= PHI1){
				area += interpolation->interpolate(PHI1, psi, lambda);
				area -= interpolation->interpolate(PHI0, psi, lambda);
			}
			else{
				area += maxArea;
				area -= interpolation->interpolate(PHI0, psi, lambda);
				area += maxArea;
				area += interpolation->interpolate(PHI1, psi, lambda);
				
			}
		}
		
	}
	
	
	return area;
	
	
	
	
}


double Integrator::angle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}

double Integrator::complAngle(Vector &a, Vector &b){
	return asin(norm_dot(a,b));
}

int Integrator::sgn(double d){
	if(d>=0) return 1;
	else return -1;
}


/**
 * n is a normal vector to the plane connecting integration origin, origin and center of circular region
 * o is a vector pointing in direction of integration origin but intersects with the interface of circular region
 * a is one of the PHI vectors from the center of circular region to its interface
 * 
 * function will return the angle between a and the plane between 0 and PI (instead of 0 and PI/2)
 * 
 */
double Integrator::complLongAngle(Vector &n, Vector &o, Vector &a){
	double v;
	int s;
	
	v = asin(norm_dot(n,a));
	s = sgn(dot(o,a));
	
	if(s<0 && v>0) v = M_PI/2.0 - v;
	else if(s<0 && v<0) v = -M_PI/2.0 - v;
	   
	return v;
		
}

/* R code
complAngle <-function(a,b){
	la = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
	lb = sqrt(b[1]^2 + b[2]^2 + b[3]^2)
	a = a/la
	b = b/lb
	
	asin(a[1]*b[1]+a[2]*b[2]+a[3]*b[3])
}

angle <-function(a,b){
	la = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
	lb = sqrt(b[1]^2 + b[2]^2 + b[3]^2)
	a = a/la
	b = b/lb
	
	acos(a[1]*b[1]+a[2]*b[2]+a[3]*b[3])
}

dot <- function(a,b){
	a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
}
*/



double Integrator::csc(double a){
	return 1.0/sin(a);
}



Vector Integrator::halfSphereIntersectionPoint(Vector &integrationOrigin, CircularRegion &c){
	double lambda;
	double psi;
	Vector v;
	double rho;
	
	lambda = c.openingAngle;
	psi = angle(integrationOrigin, c.vector);
	
	rho = acos(cos(lambda)*csc(psi));
	v(0) = 0;
	v(1) = cos(rho);
	v(2) = sin(rho);
	
	return v;
}



void Integrator::splitSASA(list<IntersectionPoint*> &sasa, vector<CircularRegion> &circles, int c, Vector &integrationOrigin, list<IntersectionPoint*>** frontHemisphere, list<IntersectionPoint*>** backHemisphere ){
	list<IntersectionPoint*>::iterator it;
	IntersectionPoint *first, *second;
	*frontHemisphere = new list<IntersectionPoint*>();
	*backHemisphere = new list<IntersectionPoint*>();
	
	IntersectionPoint ip;
	IntersectionPoint *ipp;
	
	
	
	
	
	double angleFirst, angleSecond;
	int locFirst, locSecond;
	
	second = *(--sasa.end());
	angleSecond = angle(integrationOrigin, second->vector);
	if(angleSecond<=M_PI/2){
		locSecond=0;
		(*frontHemisphere)->push_back(second);
	}
	else{
		locSecond=1;
		(*backHemisphere)->push_back(second);
	}
	
	for(it = sasa.begin(); it!=sasa.end(); ++it){
		first = second;
		locFirst = locSecond;
		second = *it;
		angleSecond = angle(integrationOrigin, second->vector);
		if(angleSecond<=M_PI/2) locSecond=0;
		else locSecond=1;

		if(locSecond != locFirst){
			ip.vector = halfSphereIntersectionPoint(integrationOrigin, circles[first->with]);
			ip.from = first->with;
			ip.with = c;
			ipp=&*(circles[c].forwardIntersections.insert(circles[c].forwardIntersections.end(),ip));
			
			(*frontHemisphere)->push_back(ipp);
			(*backHemisphere)->push_back(ipp);
		}
		
		
		if(locSecond==0)
			(*frontHemisphere)->push_back(second);
		else
			(*backHemisphere)->push_back(second);
		
		
	}
	
	
}

