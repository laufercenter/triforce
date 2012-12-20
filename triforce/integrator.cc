#include "integrator.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;


IntegratorTriforce::IntegratorTriforce(){
	
}



IntegratorTriforce::IntegratorTriforce(Interpolation *data){
	this->data = data;
	
}



double IntegratorTriforce::integrate(Molecule *m, Tessellation *tessellation){
	SASAsForMolecule sasas;
	SASANodeList sasa;
	vector<double> radii;
	Vector integrationOrigin;
	double radius;
	double area,a;
	
	this->molecule = m;
	this->tessellation = tessellation;
	
	
	
	radii = molecule->fetchRadii();
	
	for(int i=0; i<radii.size();i++){
		radius = radii[i];
		printf("RADIUS[%d]: %f\n",i,radius);
	}
		
	
	
	sasas = tessellation->sasas();
	area = 0;
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		radius = sasas[i].radius;
		a = radius*radius * integrateAtomicSASA(sasas[i]);
		area += a;
		
		
	}
	
	return area;
	
	
}


double IntegratorTriforce::integrateAtomicSASA(SASAsForAtom sasasForAtom){
	double radius;
	double area0, area1, a, area;
	
	radius = sasasForAtom.radius;
	area0 = 0;
	area1 = 0;
	for(int i=0;i<sasasForAtom.sasas.size();++i){
		//int i=1;
		
		a = integrateSASA(sasasForAtom.sasas[i]);
		
		if(sasasForAtom.sasas[i].hemisphere == FRONTHEMISPHERE)
			area0 += a;
		else area1+=a;
		
	}
	
	if(area0<0) area0 = -area0;
	else area0=2*M_PI - area0;

	if(area1<0) area1 = -area1;
	else area1=2*M_PI - area1;
	
	area = area0 + area1;
	
	return area;
	
}

double IntegratorTriforce::integrateSASA(SASA &sasa){
	SASANodeList::iterator it;
	SASANode x0, x1;
	double area=0;
	double totalAngle=0;
	
	x0 = *(--sasa.sasa.end());
	for(it = sasa.sasa.begin(); it!=sasa.sasa.end(); ++it){
		x1 = *it;
		
		area += integrateTriangle(x0, x1, sasa.tessellationOrigin, totalAngle);
		
		printf("AREA: %f totalAngle: %f\n",area,totalAngle);
		x0=x1;
		
	}
	
	printf("+++++ SAA END +++++\n\n");
	
	return area;
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




double IntegratorTriforce::calculateArea(double PHI, double psi, double lambda){
	double area;
	double aPHI;
	
	aPHI=abs(PHI);
	area = data->interpolate(aPHI, psi, lambda);
	
	if(PHI <0) area=-area;
	
	return area;
}
	
	


double IntegratorTriforce::integrateTriangle(SASANode &x0, SASANode &x1, Vector integrationOrigin, double &totalAngle){
	double psi;
	double lambda;
	double PHI0;
	double PHI1;
	double aPHI0;
	double aPHI1;
	double area,area0,area1;
	double maxArea;
	int form;
	double phi0, phi1, phi0a, phi1a, phi0b, phi1b;
	Vector n(3);
	double a2;
	
	
	area = 0;
	
	//calculate psi
	n = x1.normalForCircularInterface;
	if(x1.form!=CONVEX)
		n = -n;
	psi = angle(n, integrationOrigin);
	lambda = x1.lambda.rotation;
	
	
	
	PHI0 = x0.rotation1.rotation;
	PHI1 = x1.rotation0.rotation;
	
	aPHI0 = abs(PHI0);
	aPHI1 = abs(PHI1);

		
	maxArea = data->interpolate(M_PI, psi, lambda);	
	
	
	area0 = calculateArea(PHI0, psi, lambda);
	area1 = calculateArea(PHI1, psi, lambda);
	
	area = area1-area0;
	
	printf("PHIs: %f, %f\n",PHI0,PHI1);
	
	if(PHI1 < PHI0){
		area = 2*maxArea + area;
	}
	
	if(x1.form!=CONVEX) area=-area;

	
	return area;
	
	


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




double IntegratorTriforce::csc(double a){
	return 1.0/sin(a);
}


