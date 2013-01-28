#include "integrator.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;


IntegratorTriforce::IntegratorTriforce(){
	
}


/*
 * 409060
409379
409496


409060
*/


IntegratorTriforce::IntegratorTriforce(Interpolation *dataConcave, Interpolation *forcesConcave0, Interpolation *forcesConcave1, Interpolation *forcesConcave2){
	data.clear();
	data.push_back(dataConcave);
	data.push_back(forcesConcave0);
	data.push_back(forcesConcave1);
	data.push_back(forcesConcave2);
}
	

	
void IntegratorTriforce::clearForces(){
	int i,j;
	for(i=0; i<forces.size(); ++i)
		for(j=0; j<3; ++j){
			*(forces[i][j]) = 0;
		}
}



double IntegratorTriforce::integrate(Molecule *m, Tessellation *tessellation){
	SASAsForMolecule sasas;
	SASANodeList sasa;
	Vector integrationOrigin;
	double radius;
	double area,a;
	
	this->molecule = m;
	this->tessellation = tessellation;
	
	
	
	radii = molecule->fetchRadii();
	atoms = molecule->fetchCoordinates();
	forces = molecule->fetchForcePointers();
	areas = molecule->fetchAreaPointers();
	
	clearForces();
	
	
	for(int i=0; i<radii.size();i++){
		radius = radii[i];
	}
		
	
	
	sasas = tessellation->sasas();
	area = 0;
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		radius = sasas[i].radius;
		a = integrateAtomicSASA(i, sasas[i], radius);
		
		
		a = radius*radius * a;
		
		*(areas[i]) = a;
		area += a;
		
		
	}
	
	return area;
	
	
}


double IntegratorTriforce::integrateAtomicSASA(int l, SASAsForAtom sasasForAtom, double radius){
	double area0, area1, a, area;
	
	//radius = sasasForAtom.radius;
	area0 = 0;
	area1 = 0;
	for(int i=0;i<sasasForAtom.sasas.size();++i){
		//int i=1;
		
		a = integrateSASA(l, sasasForAtom.sasas[i], radius);
		

		
		if(sasasForAtom.sasas[i].hemisphere == FRONTHEMISPHERE)
			area0 += a;
		else area1+=a;
		
	}
	
	
	
	if(area0 < 0 && abs(area0) < THRESHOLD_NEGATIVE) area0 = 0;
	if(area0 < 0) area0 = 2*M_PI+area0;

	if(area1 < 0 && abs(area1) < THRESHOLD_NEGATIVE) area1 = 0;
	if(area1 < 0) area1 = 2*M_PI+area1;
	
	
	
	area = area0 + area1;
	
	return area;
	
}


void IntegratorTriforce::addForce(int i, Vector force){
	int j;
	if(i>=0){
		for(j=0; j<3; ++j){
			*(forces[i][j]) += force(j);
		}
	}
}



bool IntegratorTriforce::isWithinNumericalLimits(double x, double l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) return true;
	else return false;
}


Vector  IntegratorTriforce::recoverCircularInterface(Vector p, double psi_b, double lambda_b, double PHI_b0){
	Vector c(2);
	if(isWithinNumericalLimits(psi_b,0)) psi_b=0.0001;
	double denominator = (lambda_b*lambda_b + psi_b*psi_b - 2*lambda_b*psi_b*cos(PHI_b0));
	
	c(0) = (psi_b * (p(0) * psi_b - lambda_b * (p(0) * cos(PHI_b0) + p(1) * sin(PHI_b0)))) / denominator;
	c(1) = (psi_b * (p(1) * psi_b - p(1) * lambda_b * cos(PHI_b0) + p(0) * lambda_b * sin(PHI_b0))) / denominator;
	
	return c;
	
}


Vector  IntegratorTriforce::recoverCircularInterface(double psi_a, double lambda_a, double PHI_a1, double psi_b, double lambda_b, double PHI_b0){
	Vector c(2);
	
	double denominator = (lambda_b*lambda_b + psi_b*psi_b - 2*lambda_b*psi_b*cos(PHI_b0));
	
	c(0) = -(psi_b*(lambda_a*psi_b*sin(PHI_a1)-lambda_a*lambda_b*sin(PHI_a1-PHI_b0)+lambda_b*psi_a*sin(PHI_b0))) / denominator;
	c(1) = (psi_b*(psi_a*psi_b+lambda_a*psi_b*cos(PHI_a1)-lambda_b*(lambda_a*cos(PHI_a1-PHI_b0)+psi_a*cos(PHI_b0)))) / denominator;
	
	return c;
	

}

Vector  IntegratorTriforce::intersectionPoint(Vector c, double psi, double lambda, double PHI){
	Vector r(2);
	double d;
	d = norm(c,2);
	
		
	if(isWithinNumericalLimits(d,0)) {
		c(0) = 0;
		c(1) = 1;
	}
	else{
		c = c/d;
	}
	
	r(0) = c(0) - (c(0) * lambda * cos(PHI)) + (c(1)*lambda*sin(PHI));
	r(1) = c(1) - (c(1) * lambda * cos(PHI)) - (c(0)*lambda*sin(PHI));
	
	return r;
	
}

double IntegratorTriforce::integrateSASA(int l, SASA &sasa, double radius){
	SASANodeList::iterator it;
	SASANodeList::iterator it2;
	SASANode x0, x1, x2;
	double area=0;
	Area integral;
	double r_square;
	double actphi;
	double phi;
	double sign_prephi;
	Vector c(2),c2(2);
	double totalAngle,a;
	Vector p0, p1, p2, p01, p12, n01(2);
	double s;
	
	totalAngle=0;
	r_square = radius*radius;
	phi = 0;
	sign_prephi = 1;
	x0 = *(--sasa.sasa.end());
	for(it = sasa.sasa.begin(); it!=sasa.sasa.end(); ++it){
		x1 = *it;
		
		integral = integrateTriangle(l, x0, x1, sasa.tessellationOrigin, actphi);
		phi += sign_prephi * actphi;
		sign_prephi = sgn(actphi);
		

		area += integral.area;
		
		
		addForce(x0.index0, integral.force_i * r_square);
		addForce(x0.index1, integral.force_j * r_square);
		addForce(x1.index1, integral.force_k * r_square);
		addForce(l, integral.force_l * r_square);
		
		/*
		it2=it;
		++it2;
		if(it2==sasa.sasa.end()) it2 = sasa.sasa.begin();
		x2 = *it2;
		
		printf("X: (%d %d) (%d %d) (%d %d)\n",x0.index0, x0.index1, x1.index0, x1.index1, x2.index0, x2.index1);
		
		c(0) = 0;
		c(1) = x1.psi.rotation;
		if(isWithinNumericalLimits(x1.psi.rotation,0))
			c(1) = 0.0001;
		
		p0 = intersectionPoint(c, x1.psi.rotation, x1.lambda.rotation, x0.rotation1.rotation);
		printf("p0: %f %f\n",p0(0),p0(1));
		p1 = intersectionPoint(c, x1.psi.rotation, x1.lambda.rotation, x1.rotation0.rotation);
		printf("p1: %f %f\n",p1(0),p1(1));
		c2 = recoverCircularInterface(p1, x2.psi.rotation, x2.lambda.rotation, x1.rotation1.rotation);
		p2 = intersectionPoint(c2, x2.psi.rotation, x2.lambda.rotation, x2.rotation0.rotation);
		printf("p2: %f %f\n",p2(0),p2(1));
		
		p01 = p1 - p0;
		p12 = p2 - p1;
		n01(0) = p01(1);
		n01(1) = -p01(0);
		
		a = norm_dot(c,p12);
		printf("A1: %f\n",a);
		a = acos(a);
		printf("A: %f\n",a);
		s = sgn(dot(p12, n01));
		
		if(s < 0) a = -a;
		totalAngle += a;
		*/

		
		
		x0=x1;
		
	}
	
	//printf("TANGLE: %f %f\n",totalAngle, area);
	
	//s = sgn(totalAngle);
	//area = abs(area)*s;
	
	return area;
}



Vector IntegratorTriforce::lookUp(double PHI, double psi, double lambda, double &phi, CircularInterfaceForm form){
	Vector res(4);
	double aPHI, apsi, alambda;
	int i;
	
	aPHI=abs(PHI);
	apsi = max(0.0,psi);
	alambda= max(0.0,lambda);
	
	
	
	
	for(i=0;i<4;++i){
		//printf("start lookup %d (%f %f %f)\n",i,PHI,psi,lambda);
		res(i) = -data[i]->interpolate(aPHI, apsi, alambda, phi);
		//printf("end lookup %d\n",i);
	}
	
	if(PHI <0){
		res=-res;
		res(1) = -res(1);
		phi=-phi;
	}
	
	
	if(form != CONVEX){
		res(2) = -res(2);
		res(3) = -res(3);
	}
	
	
	
	
	return res;
}
	

double IntegratorTriforce::PHI2phi(double PHI, double psi, double lambda){

	return acos( 	(-cos(PHI)*cos(psi)*sin(lambda)+cos(lambda)*sin(psi)) /
			(sqrt(pow(abs(cos(psi)*sin(lambda)*sin(PHI)),2) + pow(abs(sin(lambda)*sin(PHI)*sin(psi)),2)
			+ pow(abs(cos(PHI)*cos(psi)*sin(lambda) - cos(lambda) * sin(psi)),2))));


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




Area IntegratorTriforce::integrateTriangle(int l, SASANode &x0, SASANode &x1, Vector integrationOrigin, double &phi){
	Rotation psi;
	Rotation lambda;
	Rotation PHIij;
	Rotation PHIjk;
	double area,area0,area1;
	Vector M;
	Area a;
	Vector Tij(4), Tjk(4);
	Vector Tij2(4), Tjk2(4);
	Vector force_i(3), force_j(3), force_k(3), force_l(3);
	double phi0,phi1;
	double q, q0, q1;
	CircularInterfaceForm form;
	CircularInterfaceForm formi;
	CircularInterfaceForm formj;
	CircularInterfaceForm formk;
	area = 0;
	double s_convex, s_complementation, s_direction;
	
	psi = x1.psi;
	lambda = x1.lambda;
	form = x1.form;
	
	formi = x0.form0;
	formj = x0.form1;
	formk = x1.form1;
	

	
	q0 = 1;
	if(formi != CONVEX) q0*=-1;
	if(formj != CONVEX) q0*=-1;
	
	q1 = 1;
	if(formj != CONVEX) q1*=-1;
	if(formk != CONVEX) q1*=-1;
	
	
	
	
	
	PHIij = x0.rotation1;
	PHIjk = x1.rotation0;
	
	
	
	
	
	s_convex=1;
	if(form!=CONVEX){
		PHIij.rotation *=-1;
		PHIjk.rotation *=-1;
		s_convex=-1;
		
	}
	
		
	M = lookUp(M_PI, psi.rotation, lambda.rotation, q0, form);
	
	
	Tij = lookUp(PHIij.rotation, psi.rotation, lambda.rotation, phi0, form);
	Tjk = lookUp(PHIjk.rotation, psi.rotation, lambda.rotation, phi1, form);
	
	
	if(PHIjk.rotation >= PHIij.rotation){
		s_direction=1;
	}
	else{
		s_direction=-1;
	}
	
	
	area = s_direction*Tjk(0) - s_direction*Tij(0);
	
	q= s_convex * s_direction;
	
	force_i = -q*(Tij(1) * PHIij.drotation_dxj);
	force_j = q*(Tjk(1) * PHIjk.drotation_dxi + Tjk(2) * psi.drotation_dxi + Tjk(3) * lambda.drotation_dxi) - q*(Tij(1) * PHIij.drotation_dxi + Tij(2) * psi.drotation_dxi + Tij(3) * lambda.drotation_dxi);
	force_k = q*(Tjk(1) * PHIjk.drotation_dxj);
	force_l = q*((Tjk(1) * PHIjk.drotation_dxl + Tjk(2) * psi.drotation_dxl + Tjk(3) * lambda.drotation_dxl) - (Tij(1) * PHIij.drotation_dxl + Tij(2) * psi.drotation_dxl + Tij(3) * lambda.drotation_dxl));
	
	
	if(PHIjk.rotation < PHIij.rotation){
		area = 2*M(0) - area;
		force_i = -force_i;
		force_j = 2*(M(2)*psi.drotation_dxi + M(3)*lambda.drotation_dxi) - force_j;
		force_k = -force_k;
		force_l = 2*(M(2)*psi.drotation_dxl + M(3)*lambda.drotation_dxl) - force_l;
	}
	
	
	if(x1.form==CONVEX) s_convex=1;
	else s_convex=-1;
	
	q = s_convex;
	
	a.area=s_convex * area;
	a.force_i = q* force_i;
	a.force_j = q* force_j;
	a.force_k = q* force_k;
	a.force_l = q* force_l;
	
	


	
	return a;

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


