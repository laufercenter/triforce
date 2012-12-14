#include "integrator.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;


IntegratorTriforce::IntegratorTriforce(){
	
}



IntegratorTriforce::IntegratorTriforce(	Interpolation *dataConvex, Interpolation *forcesConvex0, Interpolation *forcesConvex1, Interpolation *forcesConvex2,
					Interpolation *dataConcave, Interpolation *forcesConcave1, Interpolation *forcesConcave2, Interpolation *forcesConcave3)
{
	
	convex.push_back(dataConvex);
	convex.push_back(forcesConvex0);
	convex.push_back(forcesConvex1);
	convex.push_back(forcesConvex2);

	convex.push_back(dataConcave);
	concave.push_back(forcesConcave0);
	concave.push_back(forcesConcave1);
	concave.push_back(forcesConcave2);
	
}



double IntegratorTriforce::integrate(Molecule *m, Tessellation *tessellation){
	SASAsForMolecule sasas;
	SASANodeList sasa;
	vector<double> *radii;
	Vector integrationOrigin;
	double radius;
	double area,result;
	
	this->molecule = m;
	this->tessellation = tessellation;
	
	
	
	radii = molecule->fetchRadii();
	forces = molecule->fetchForcePointers();
	areas = molecule->fetchAreaPointers();
	
	for(int i=0; i<radii->size();i++){
		radius = radii->at(i);
		printf("RADIUS[%d]: %f\n",i,radius);
	}
		
	
	
	sasas = tessellation->sasas();
	area = 0;
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		radius = sasas[i].radius;
		result = integrateAtomicSASA(sasas[i]);
		a = radius*radius * result(0);
		
		*(areas[i]) = a;
		area += a;
		
		for(j=0;j<3;++j)
			*(forces[i][j]) = result(j+1);
			
		
		
	}
	
	return area;
	
	
}


double IntegratorTriforce::integrateAtomicSASA(SASAsForAtom sasasForAtom){
	double radius;
	Vector result0(4), result1(4), a(4), result(4);
	
	radius = sasasForAtom.radius;
	result0 = 0;
	result1 = 0;
	for(int i=0;i<sasasForAtom.sasas.size();++i){
		//int i=1;
		
		a = integrateSASA(sasasForAtom.sasas[i]);
		
		if(sasasForAtom.sasas[i].hemisphere == FRONTHEMISPHERE)
			result0 += a;
		else result1+=a;
		
	}
	
	if(result0(0)<0) result0 = -result0;
	else result0(0)=2*M_PI - result0(0);

	if(result1(0)<0) result1 = -result1;
	else result1(0)=2*M_PI - result1(0);
	
	result = result0 + result1;
	
	return result;
	
}

double IntegratorTriforce::integrateSASA(SASA &sasa){
	SASANodeList::iterator it;
	SASANode x0, x1;
	Vector result(4);
	result(0)=0;
	result(1)=0;
	result(2)=0;
	result(3)=0;
	
	double totalAngle=0;
	
	x0 = *(--sasa.sasa.end());
	for(it = sasa.sasa.begin(); it!=sasa.sasa.end(); ++it){
		x1 = *it;
		
		result += integrateTriangle(x0, x1, sasa.tessellationOrigin, totalAngle);
		
		printf("AREA: %f totalAngle: %f\n",result(0),totalAngle);
		x0=x1;
		
	}
	
	printf("+++++ SAA END +++++\n\n");
	
	return result;
}


double IntegratorTriforce::PHI2phi(double PHI, double psi, double lambda){

	return acos( 	(-cos(PHI)*cos(psi)*sin(lambda)+cos(lambda)*sin(psi)) /
			(sqrt(pow(abs(cos(psi)*sin(lambda)*sin(PHI)),2) + pow(abs(sin(lambda)*sin(PHI)*sin(psi)),2)
			+ pow(abs(cos(PHI)*cos(psi)*sin(lambda) - cos(lambda) * sin(psi)),2))));


}


mat33 IntegratorTriforce::rotz(double theta){
	mat33 m;
	m(0,0) = cos(theta);
	m(0,1) = -sin(theta);
	m(0,2) = 0;
	
	m(1,0) = sin(theta);
	m(1,1) = cos(theta);
	m(1,2) = 0;
	
	m(2,0) = 0;
	m(2,1) = 0;
	m(2,2) = 1;
	
	return m;
		
}

double IntegratorTriforce::PHI2phi2(Vector integrationOrigin, double PHI, double psi, double lambda){
	mat33 T;
	mat33 r0, r1;
	Vector n(3);
	Vector n0(3);
	Vector n1(3);
	Vector v(3),v2(2);
	Vector ex(3);
	Vector p(3);
	double ux,uy,uz;
	double C,S,t;
	double g;
	
	ex(0) = 1;
	ex(1) = 0;
	ex(2) = 0;

	r0 = rotz(psi);
	r1 = rotz(psi-lambda);
	
	g = 1-cos(lambda);

	n = r0 * ex;
	n0=g*n;
	n1 = r1 * ex;
	
	v = n1-n0;
	
	ux=n(0);
	uy=n(1);
	uz=n(2);
	C=cos(PHI);
	S=sin(PHI);
	t=1-cos(PHI);
	
	
	T(0,0) = t*ux*ux + C;
	T(0,1) = t*ux*uy - S*uz;
	T(0,2) = t*ux*uz + S*uy;
	
	T(1,0) = t*ux*uy + S*uz;
	T(1,1) = t*uy*uy + C;
	T(1,2) = t*uy*uz - S*ux;
	
	T(2,0) = t*ux*uz - S*uy;
	T(2,1) = t*uy*uz + S*ux;
	T(2,2) = t*uz*uz + C;
	
	v2 = T * v;
	
	p= n0 + v2;
	
	return V2phi(integrationOrigin, n0, p);
	
	


}

bool IntegratorTriforce::isInPositiveEpsilonRange(double v, double eps){
	if(eps-(v+THRESHOLD_NUMERICAL) <= 0) return true;
	else return false;
}


double IntegratorTriforce::V2phi(Vector &integrationOrigin, Vector cv, Vector &v){
	Vector up(3);
	Vector n_origin(3);
	Vector n_v(3);
	Vector n_cv(3);

	
	up(0) = 0;
	up(1) = 1;
	up(2) = 0;
	
	
	if(isInPositiveEpsilonRange(fabs(norm_dot(integrationOrigin,cv)),1.0) ){
		cv = up;
	}
	
	
	
	n_origin = cross(up,integrationOrigin);
	n_v = cross(v, integrationOrigin);
	n_cv = cross(cv, integrationOrigin);
	
	return abs(acos(norm_dot(n_cv,n_v)));
}


Vector convexeLookUp(double PHI, double psi, double lambda){
	Vector r(4);
	for(i=0; i<4; ++i)
		r(i) = convex[i]->interpolate(PHI,psi,lambda);		
	return r;
}


Vector concaveLookUp(double PHI, double psi, double lambda){
	Vector r(4);
	for(i=0; i<4; ++i)
		r(i) = concave[i]->interpolate(PHI,psi,lambda);		
	return r;
}


Vector IntegratorTriforce::integrateTriangle(SASANode &x0, SASANode &x1, Vector integrationOrigin, double &totalAngle){
	double psi;
	double lambda;
	double PHI0;
	double PHI1;
	double aPHI0;
	double aPHI1;
	double v0,v1;
	double area;
	Vector maxResult(4);
	int form;
	double phi0, phi1, phi0a, phi1a, phi0b, phi1b;
	Vector a(4);
	Vector result(4);
	Vector force(3);
	
	
	area = 0;
	
	//calculate psi
	n = x1.normalForCircularRegion;
	if(x1.form!=CONVEX)
		n = -n;
	//psi = angle(n, integrationOrigin);
	psi = x1.psi.rotation;
	
	//calculate lambda
	lambda = x1.lambda.rotation;
	
	
	//printf("SEGMENT (%d,%d) - (%d,%d)\n",x0.id0,x0.id1,x1.id0,x1.id1);
	
	PHI0 = x0.rotation1.rotation;
	PHI1 = x1.rotation0.rotation;
	
	aPHI0 = abs(PHI0);
	aPHI1 = abs(PHI1);
	
	phi0 = abs(PHI2phi(PHI0,psi,lambda));
	phi0a = V2phi(integrationOrigin, x1.normalForCircularRegion, x0.vector);
	phi0b = PHI2phi2(integrationOrigin, PHI0, psi, lambda);
	if(PHI0<0){
		phi0=-phi0;
		phi0a=-phi0a;
		phi0b=-phi0b;
	}
	
	phi1 = abs(PHI2phi(PHI1,psi,lambda));
	phi1a = V2phi(integrationOrigin, x1.normalForCircularRegion, x1.vector);
	phi1b = PHI2phi2(integrationOrigin, PHI1, psi, lambda);
	if(PHI1<0) {
		phi1=-phi1;
		phi1a=-phi1a;
		phi1b=-phi1b;
	}
	
	
	if(x1.form!=CONVEX)
		totalAngle -= phi1 - phi0;
	else
		totalAngle += phi1 - phi0;
	
	
	
		
	if(psi+lambda>=M_PI){
		printf("MAXAREA 0\n");
		maxResult = concaveLookUp(M_PI, psi, lambda);	
	}
	else{
		printf("MAXAREA 1\n");
		maxResult = convexLookUp(0, psi, lambda);
	}
	printf("MAXAREA: %f, psi %f, lambda %f, PHI0 %f, PHI1 %f, phi0: %f (%f) [%f], phi1: %f (%f) [%f],totalAngle: %f, CIRCLE: %d\n",maxResult,psi,lambda,PHI0,PHI1,phi0,phi0a,phi0b,phi1,phi1a,phi1b,totalAngle,x0.id1);
	
	if(PHI0 >= 0){
		if(PHI1 >= 0){
			if(PHI0 <= PHI1){
				if(psi<=lambda){
					printf("CASE 0 A\n");
					a = maxResult;
					a -= convexLookUp(aPHI1, psi, lambda);
					
					result += a;
					printf("CASE 0 0 A %f\n",result);
					
					a = maxResult;
					a -= convexLookUp(aPHI0, psi, lambda);
					result -= a;
					
					printf("CASE 0 1 A %f\n",result);
				}
				else{
					printf("CASE 0 B\n");
					result += concaveLookUp(aPHI1, psi, lambda);
					printf("CASE 0 0 B %f\n",result);
					result -= concaveLookUp(aPHI0, psi, lambda);
					printf("CASE 0 1 B% f\n",result);
					
				}
			}
			else{
				if(psi<=lambda){
					printf("CASE 1 A\n");
					a = maxResult;
					a -= convexLookUp(aPHI1, psi, lambda);
					result += a;
					printf("CASE 1 0 A %f\n",result);
					
					result += convexLookUp(aPHI0, psi, lambda);
					printf("CASE 1 1 A %f\n",result);
					result += maxResult;
					printf("CASE 1 2 A %f\n",result);
					
				}
				else if(psi+lambda>=M_PI){
					printf("CASE 1 B\n");
					
					a = maxResult;
					a -= concaveLookUp(aPHI0, psi, lambda);
					result += a;
					printf("CASE 1 0 B%f\n",result);

					result += concaveLookUp(aPHI1, psi, lambda);
					printf("CASE 1 1 B %f\n",result);
					
					result += maxResult;
					printf("CASE 1 2 B %f\n",result);
					
				}
				else{
					printf("CASE 1 C\n");
					
					result += concaveLookUp(aPHI1, psi, lambda);
					printf("CASE 1 0 C %f\n",result);
					
					result += convexLookUp(aPHI0, psi, lambda);
					printf("CASE 1 1 C %f\n",result);
					
					result += maxResult;
					printf("CASE 1 2 C %f\n",result);
					
				}
				
			}
		}
		else{
			if(psi <= lambda){
				printf("CASE 2 A\n");
				result += convexLookUp(aPHI0, psi, lambda);
					printf("CASE 2 0 A %f\n",result);
				result += convexLookUp(aPHI1, psi, lambda);
					printf("CASE 2 1 A %f\n",result);
			}
			else if(psi+lambda>=M_PI){
					printf("CASE 2 B\n");
					a = maxResult;
					a -= concaveLookUp(aPHI0, psi, lambda);
					result += a;
					printf("CASE 2 0 B %f\n",result);
					
					a = maxResult;
					a -= concaveLookUp(aPHI1, psi, lambda);
					result += a;
					printf("CASE 2 1 B %f\n",result);
					
			}
			else{
				printf("CASE 2 C\n");
				result += convexLookUp(aPHI0, psi, lambda);
				printf("CASE 2 0 C %f\n",result);
				result += convexLookUp(aPHI1, psi, lambda);
				printf("CASE 2 1 C %f\n",result);
			}
				
		}
	}
	else{
		if(PHI1 >= 0){
			if(psi<=lambda){
				printf("CASE 3 A\n");
				a = maxResult;
				a -= convexLookUp(aPHI0, psi, lambda);
				result += a;
				printf("CASE 3 0 A %f\n",result);
				
				a = maxResult;
				a -= convexLookUp(aPHI1, psi, lambda);
				result += a;
				printf("CASE 3 1 A %f\n",result);
			}
			else if(psi+lambda>=M_PI){
				printf("CASE 3 B\n");
				result += concaveLookUp(aPHI0, psi, lambda);
				printf("CASE 3 0 B %f\n",result);
				result += concaveLookUp(aPHI1, psi, lambda);
				printf("CASE 3 1 B %f\n",result);
				
			}
			else{
				printf("CASE 3 C\n");
				result += concaveLookUp(aPHI0, psi, lambda);
				printf("CASE 3 0 C %f\n",result);
				result += concaveLookUp(aPHI1, psi, lambda);
				printf("CASE 3 1 C %f\n",result);
			}
			
		}
		else{
			if(PHI0 <= PHI1){
				if(psi<=lambda){
					printf("CASE 4 A\n");
					a = maxResult;
					a -= convexLookUp(aPHI0, psi, lambda);
					result += a;
					printf("CASE 4 0 A %f\n",result);
					
					a = maxResult;
					a -= convexLookUp(aPHI1, psi, lambda);
					result -= a;
					printf("CASE 4 1 A %f\n",result);
				}
				else if(psi+lambda>=M_PI){
					printf("CASE 4 B\n");
					result += concaveLookUp(aPHI0, psi, lambda);
					printf("CASE 4 0 B %f\n",result);
					result -= concaveLookUp(aPHI1, psi, lambda);
					printf("CASE 4 1 B %f\n",result);
					
				}
				else{
					printf("CASE 4 C\n");
					result += concaveLookUp(aPHI0, psi, lambda);
					printf("CASE 4 0 C %f\n",result);
					result -= concaveLookUp(aPHI1, psi, lambda);
					printf("CASE 4 1 C %f\n",result);
				}
			}
			else{
				if(psi<=lambda){
					printf("CASE 5 A\n");
					a = maxResult;
					a -= convexLookUp(aPHI0, psi, lambda);
					result += a;
					printf("CASE 5 0 A %f\n",result);
					
					result += maxResult;
					printf("CASE 5 1 A %f\n",result);
					
					result += convexLookUp(aPHI1, psi, lambda);
					printf("CASE 5 2 A %f\n",result);
					
				}
				else if(psi+lambda>=M_PI){
					printf("CASE 5 B\n");
					result += concaveLookUp(aPHI0, psi, lambda);
					printf("CASE 5 0 B %f\n",result);
					result += maxResult;
					printf("CASE 5 1 B %f\n",result);
					
					a = maxResult;
					a -= concaveLookUp(aPHI1, psi, lambda);
					result += a;
					printf("CASE 5 2 B %f\n",result);
					
					
				}
				else{
					printf("CASE 5 C\n");
					result += concaveLookUp(aPHI0, psi, lambda);
					printf("CASE 5 0 C %f\n",result);
					result += maxResult;
					printf("CASE 5 1 C %f\n",result);
					result += convexLookUp(aPHI1, psi, lambda);
					printf("CASE 5 2 C %f\n",result);
					
				}
			}
		}
		
	}
	
	
	
	
	if(x1.form!=CONVEX) s = -1;
	else s=1;
	
	
	result= s*result;
	
	force(0)
	
	
	
	
	return result;
	
	
	
	
}




double IntegratorTriforce::angle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}

double IntegratorTriforce::complAngle(Vector &a, Vector &b){
	return asin(norm_dot(a,b));
}

int IntegratorTriforce::sgn(double d){
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



double IntegratorTriforce::complLongAngle(Vector &n, Vector &o, Vector &a){
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



double IntegratorTriforce::csc(double a){
	return 1.0/sin(a);
}


/*
Vector IntegratorTriforce::halfSphereIntersectionPoint(Vector &integrationOrigin, CircularRegion &c, double radius, int sign){
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
void IntegratorTriforce::splitSASA(list<IntersectionNode*> &sasa, vector<CircularRegion> &circles, int c, Vector &integrationOrigin, double radius, list<IntersectionNode*>** frontHemisphere, list<IntersectionNode*>** backHemisphere, IntersectionGraph &intersectionGraph){
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




