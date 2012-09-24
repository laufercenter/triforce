#include "integrator.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;


Integrator::Integrator(){
	
}



Integrator::Integrator(Tessellation *tessellation, Interpolation *interpolation, Molecule *m){
	this->tessellation = tessellation;
	this->interpolation = interpolation;
	this->molecule = m;
	
}



double Integrator::integrate(){
	SASAs sasas;
	SASANodeList sasa;
	vector<double> *radii;
	Vector integrationOrigin;
	double radius;
	double area=0;
	radii = molecule->fetchRadii();
	
	for(int i=0; i<radii->size();i++){
		radius = radii->at(i);
		printf("RADIUS[%d]: %f\n",i,radius);
	}
		
	
	
	sasas = tessellation->sasas();
 
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		radius = sasas[i].radius;
		sasa=sasas[i].sasa;
		//iterate over all sasas
		printf("size: %d\n",sasa.size());
		for(int j=0;j<sasa.size();++j){
			area += integrateSASA(sasas[i]);
		}
		
	}
	
	return area;
	
	
}

double Integrator::integrateSASA(SASA &sasa){
	SASANodeList::iterator it;
	SASANode x0, x1;
	double area=0;
	
	x0 = *(--sasa.sasa.end());
	for(it = sasa.sasa.begin(); it!=sasa.sasa.end(); ++it){
		x1 = *it;
		
		area += integrateTriangle(x0, x1, sasa.tessellationOrigin);
		
		printf("area: %f\n",area);
		x0=x1;
		
	}
	
	return area;
}



double Integrator::integrateTriangle(SASANode &x0, SASANode &x1, Vector integrationOrigin){
	double psi;
	double lambda;
	double PHI0;
	double PHI1;
	double v0,v1;
	double area = 0;
	Vector vec0,vec1;
	Vector n(3), o;
	double maxArea;
	
	
	
	//calculate psi
	psi = angle(x1.normalForCircularRegion, integrationOrigin);
	//calculate lambda
	lambda = x1.lambda;
	
	
	printf("SEGMENT (%d,%d) - (%d,%d)\n",x0.id0,x0.id1,x1.id0,x1.id1);
	
	PHI0 = x0.angle1;
	PHI1 = x1.angle0;
	
	
	
	
	maxArea = interpolation->interpolate(0, psi, lambda);
	printf("MAXAREA: %f, psi %f, lambda %f, PHI0 %f, PHI1 %f, CIRCLE: %d\n",maxArea,psi,lambda,PHI0,PHI1,x0.id1);
	
	if(PHI0 >= 0){
		if(PHI1 >= 0){
			if(PHI0 <= PHI1){
				area += interpolation->interpolate(PHI0, psi, lambda);
				printf("CASE 0 0 %f\n",area);
				area -= interpolation->interpolate(PHI1, psi, lambda);
				printf("CASE 0 1 %f\n",area);
			}
			else{
				area += interpolation->interpolate(PHI0, psi, lambda);
				printf("CASE 1 0 %f\n",area);
				area -= interpolation->interpolate(PHI1, psi, lambda);
				printf("CASE 1 1 %f\n",area);
				
				area += maxArea;
				area += maxArea;
				
			}
		}
		else{
			area += interpolation->interpolate(PHI0, psi, lambda);
				printf("CASE 2 0 %f\n",area);
			area += interpolation->interpolate(PHI1, psi, lambda);
				printf("CASE 2 1 %f\n",area);
		}
	}
	else{
		if(PHI1 >= 0){
			area -= interpolation->interpolate(PHI0, psi, lambda);
				printf("CASE 3 0 %f\n",area);
			area -= interpolation->interpolate(PHI1, psi, lambda);
				printf("CASE 3 1 %f\n",area);
			area += maxArea;
			area += maxArea;
			
		}
		else{
			if(PHI0 <= PHI1){
				area += interpolation->interpolate(PHI1, psi, lambda);
				printf("CASE 4 0 %f\n",area);
				area -= interpolation->interpolate(PHI0, psi, lambda);
				printf("CASE 4 1 %f\n",area);
			}
			else{
				area -= interpolation->interpolate(PHI0, psi, lambda);
				printf("CASE 5 0 %f\n",area);
				area += interpolation->interpolate(PHI1, psi, lambda);
				printf("CASE 5 1 %f\n",area);
				
				area += maxArea;
				area += maxArea;
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
 * 
 * 
complAngle <-function(nij, nik, vi){
	v = acos(nij[1]*nik[1]+nij[2]*nik[2]+nij[3]*nik[3])
	ni = cross(vi,nij)
	s = sign(ni[1]*nik[1]+ni[2]*nik[2]+ni[3]*nik[3])
	
	if(s<0) v =  -v
	   
	v
		
}


n <-function(a){
	la = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
	a=a/la
	a
	}



aii <- function(a,n){
	
asin(n[1]*a[1]+n[2]*a[2]+n[3]*a[3])
}

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

cross <- function(a,b){
	x=c(0,0,0)
	x[1]=a[2]*b[3] - a[3]*b[2]
	x[2]=a[3]*b[1] - a[1]*b[3]
	x[3]=a[1]*b[2] - a[2]*b[1]
	
	x
}
*/



double Integrator::csc(double a){
	return 1.0/sin(a);
}


/*
Vector Integrator::halfSphereIntersectionPoint(Vector &integrationOrigin, CircularRegion &c, double radius, int sign){
	double lambda;
	double psi;
	Vector v(3);
	double rho;
	Vector orthoIntegrationOrigin(3);
	double l2;
	Vector n(3);
	double nu;
	
	orthoIntegrationOrigin(0) = 0;
	orthoIntegrationOrigin(1) = 1;
	orthoIntegrationOrigin(2) = 0;
	
	lambda = c.openingAngle;
	psi = angle(integrationOrigin, c.normal);
	rho = sign * acos(cos(lambda)*csc(psi));
	
	n = cross(c.normal, integrationOrigin);
	nu = complAngle(n,orthoIntegrationOrigin);
	

	
	v(0) = 0;
	v(1) = cos(rho + nu);
	v(2) = sin(rho + nu);
	
	v*=radius;
	
	
	
	
	return v;
}

*/
/*
 * sasa 1 size 4
intersectionpoint 0 vector 1.256455 0.354529 1.150268
intersectionpoint 1 vector 1.269964 0.264348 1.159653
intersectionpoint 2 vector 1.371810 0.136769 -1.061563
intersectionpoint 3 vector -0.114692 1.732367 -0.115048

 * 
 * 
 * 
 * */

/*
void Integrator::splitSASA(list<IntersectionNode*> &sasa, vector<CircularRegion> &circles, int c, Vector &integrationOrigin, double radius, list<IntersectionNode*>** frontHemisphere, list<IntersectionNode*>** backHemisphere, IntersectionGraph &intersectionGraph){
	list<IntersectionNode*>::iterator it;
	IntersectionNode *first, *second;
	*frontHemisphere = new list<IntersectionNode*>();
	*backHemisphere = new list<IntersectionNode*>();
	Vector mirrorIntegrationOrigin(3);
	Vector orthoIntegrationOrigin(3);
	double angleOrtho;
	Interfaces interfaces;
	Vector p0(3), p1(3);
	Vector v2(3);
	CircularRegion C, I;
	
	IntersectionNode ip;
	IntersectionNode *ipp;
	IntersectionAddress a;
	
	double angleFirst, angleSecond;
	int locFirst, locSecond;
	int sign;

	
	C = circles[c];
	
	
	mirrorIntegrationOrigin(0) = -integrationOrigin(0);
	mirrorIntegrationOrigin(1) = -integrationOrigin(1);
	mirrorIntegrationOrigin(2) = -integrationOrigin(2);
	
	orthoIntegrationOrigin(0) = 0;
	orthoIntegrationOrigin(1) = 1;
	orthoIntegrationOrigin(2) = 0;
	
	
	
	
	second = *(--sasa.end());
	angleSecond = angle(integrationOrigin, second->vector);
	if(angleSecond<=M_PI/2){
		locSecond=0;
		//(*frontHemisphere)->push_back(second);
	}
	else{
		locSecond=1;
		//(*backHemisphere)->push_back(second);
	}
	
	for(it = sasa.begin(); it!=sasa.end(); ++it){
		first = second;
		locFirst = locSecond;
		second = *it;
		angleSecond = angle(integrationOrigin, second->vector);
		angleOrtho = angle(orthoIntegrationOrigin, second->vector);
		if(angleSecond<=M_PI/2){
			locSecond=0;
			
			sign = -1;
			
		}
		else{
			locSecond=1;
			
			sign = 1;
		}
		
		

		if(locSecond != locFirst){
			I = circles[first->id1];
			ip.vector = halfSphereIntersectionPoint(integrationOrigin, I, radius, sign);
			
			//tessellation->measurementPoints(p0, p1, integrationOrigin, I);
			//v2 = I.normal * I.g;

			//interface = tessellation->angularInterface(ip.vector, v2, p0, p1);
			
			if(locSecond==1){
				a.id0 = first->id1;
				a.id1 = c;

				ip.id0 = first->id1;
				ip.id1 = c;
			}
			else{
				a.id0 = c;
				a.id1 = first->id1;
				
				ip.id0 = c;
				ip.id1=first->id1;
			}
			
			intersectionGraph[a] = ip;
			ipp = &intersectionGraph[a];
			
			(*frontHemisphere)->push_back(ipp);
			(*backHemisphere)->push_back(ipp);
		}
		
		
		if(locSecond==0)
			(*frontHemisphere)->push_back(second);
		else
			(*backHemisphere)->push_back(second);
		
		
	}
	
	
}
*/




