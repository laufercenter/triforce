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


	data.clear();
	data.push_back(dataConcave);
	data.push_back(forcesConcave0);
	data.push_back(forcesConcave1);
	data.push_back(forcesConcave2);

*/


IntegratorTriforce::IntegratorTriforce(vector<Interpolator*> data){
	this->data=data;
	benchmark=Benchmark("integrator");
	
}
	

	
void IntegratorTriforce::clearForces(){
	purgeForces();
	
	for(unsigned int i=0; i<forces.size(); ++i)
		for(unsigned int j=0; j<3; ++j){
			*(forces[i][j]) = 0;
		}
}


void IntegratorTriforce::pushForces(){
	for(unsigned int i=0; i<forcesDelayed.size(); ++i)
		for(unsigned int j=0; j<3; ++j){
			*(forces[forcesDelayed[i].i][j]) += static_cast<ForcesDT>(forcesDelayed[i].force(j));
		}
}

void IntegratorTriforce::purgeForces(){
	forcesDelayed.clear();
}



float IntegratorTriforce::integrate(Molecule *m, Tessellation *tessellation){
	SASAs sasas;
	Vector integrationOrigin;
	float radius;
	float area,a;
	
	
	this->molecule = m;
	this->tessellation = tessellation;
	
	
	
	radii = molecule->fetchRadii();
	atoms = molecule->fetchCoordinates();
	forces = molecule->fetchForcePointers();
	areas = molecule->fetchAreaPointers();
	
	clearForces();
	
	
	
	
	sasas = tessellation->sasas();
	area = 0;
	
	//iterate over all atoms
	
	benchmark.start(string("integrating"));
	
	for(unsigned int i=0;i<sasas.size();++i){
		radius = radii[i];
		tradius=radius;
		a = integrateSASA(i, sasas[i], radius);
		
		
		a = radius*radius * a;
		
		*(areas[i]) = static_cast<AreasDT>(a);
		area += a;
		
		if(a>0) pushForces();
		purgeForces();
		
	}
	
	benchmark.stop();
	
	return area;
	
	
}




void IntegratorTriforce::addForce(int i, Vector force){
	ForceElement fe;
	if(i>=0){
// 		if(i==749) printf("FORCE: %f %f %f\n",force(0),force(1),force(2));
		fe.i=i;
		fe.force=force;
		forcesDelayed.push_back(fe);
	}
}



bool IntegratorTriforce::isWithinNumericalLimits(float x, float l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) return true;
	else return false;
}


Vector  IntegratorTriforce::recoverCircularInterface(Vector p, float psi_b, float lambda_b, float PHI_b0){
	Vector c(2);
	if(isWithinNumericalLimits(psi_b,0)) psi_b=0.0001;
	float denominator = (lambda_b*lambda_b + psi_b*psi_b - 2*lambda_b*psi_b*cos(PHI_b0));
	
	c(0) = (psi_b * (p(0) * psi_b - lambda_b * (p(0) * cos(PHI_b0) + p(1) * sin(PHI_b0)))) / denominator;
	c(1) = (psi_b * (p(1) * psi_b - p(1) * lambda_b * cos(PHI_b0) + p(0) * lambda_b * sin(PHI_b0))) / denominator;
	
	return c;
	
}


Vector  IntegratorTriforce::recoverCircularInterface(float psi_a, float lambda_a, float PHI_a1, float psi_b, float lambda_b, float PHI_b0){
	Vector c(2);
	
	float denominator = (lambda_b*lambda_b + psi_b*psi_b - 2*lambda_b*psi_b*cos(PHI_b0));
	
	c(0) = -(psi_b*(lambda_a*psi_b*sin(PHI_a1)-lambda_a*lambda_b*sin(PHI_a1-PHI_b0)+lambda_b*psi_a*sin(PHI_b0))) / denominator;
	c(1) = (psi_b*(psi_a*psi_b+lambda_a*psi_b*cos(PHI_a1)-lambda_b*(lambda_a*cos(PHI_a1-PHI_b0)+psi_a*cos(PHI_b0)))) / denominator;
	
	return c;
	

}

Vector  IntegratorTriforce::intersectionPoint(Vector c, float psi, float lambda, float PHI){
	Vector r(2);
	float d;
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

float IntegratorTriforce::integrateSASA(int l, SASASegmentList &sasa, float radius){
	SASASegmentList::iterator it;
	SASASegment x;
	float area0,area1,area;
	Integral integral;
	float r_square;
	Vector c(2),c2(2);
	float a;
	Vector p0, p1, p2, p01, p12, n01(2);
	
	r_square = radius*radius;
	
	area0=0;
	area1=0;
	for(it = sasa.begin(); it!=sasa.end(); ++it){
		x=*it;
		integral = integrateTriangle(l, x, x.tessellationAxis.v);
		
		a = integral.integral;
		
		addForce(x.index0, integral.force_i * r_square);
		addForce(x.index1, integral.force_j * r_square);
		addForce(x.index2, integral.force_k * r_square);
		addForce(l, integral.force_l * r_square);
		if(x.tessellationAxis.mode==ALIGNED){
			addForce(x.tessellationAxis.index, integral.force_t * r_square);
		}
		
		

		if(x.hemisphere == FRONTHEMISPHERE)
			area0 += a;
		else area1+=a;
		
		
		
		
		
		
	}
	
	if(area0 < 0 && abs(area0) < THRESHOLD_NEGATIVE) area0 = 0;
	if(area0 < 0) area0 = 2*M_PI+area0;

	if(area1 < 0 && abs(area1) < THRESHOLD_NEGATIVE) area1 = 0;
	if(area1 < 0) area1 = 2*M_PI+area1;
	
	
	
	
	area = area0 + area1;
	
	//if we have information from the depthbuffer, check if the calculated size coarsely matches the one from the buffer
	//and if not, disregard the sasa (!)
// 	if(sasa.size()>0){
// 		if(sasa[0].depthBufferEstimatedArea>=0)
// 			if(area> 2*sasa[0].depthBufferEstimatedArea*4*M_PI) area=0;
// 	}
	
	
	
	return area;
}



Vector IntegratorTriforce::lookUp(float PHI, float psi, float lambda, CircularInterfaceForm form){
	Vector res(4);
	float aPHI, apsi, alambda;
	int i;
	Vector v(3);

	aPHI=abs(PHI);
	apsi = max(0.0f,psi);
	alambda= max(0.0f,lambda);
	
	v(0)=aPHI;
	v(1)=apsi;
	v(2)=alambda;
	
	
	for(i=0;i<4;++i){
		//printf("start lookup %d (%f %f %f)\n",i,PHI,psi,lambda);
		res(i) = -(data[i]->interpolate(v));
		//printf("end lookup %d\n",i);
	}
	
	if(PHI <0){
		res=-res;
		res(1) = -res(1);
	}
	
	
	if(form != CONVEX){
		res(2) = -res(2);
		res(3) = -res(3);
	}
	
	
	
	
	return res;
}
	

float IntegratorTriforce::PHI2phi(float PHI, float psi, float lambda){

	return acos( 	(-cos(PHI)*cos(psi)*sin(lambda)+cos(lambda)*sin(psi)) /
			(sqrt(pow(abs(cos(psi)*sin(lambda)*sin(PHI)),2) + pow(abs(sin(lambda)*sin(PHI)*sin(psi)),2)
			+ pow(abs(cos(PHI)*cos(psi)*sin(lambda) - cos(lambda) * sin(psi)),2))));


}


float IntegratorTriforce::PHI2phi2(Vector integrationOrigin, float PHI, float psi, float lambda){
	fmat33 T;
	fmat33 r0, r1;
	Vector n(3);
	Vector n0(3);
	Vector n1(3);
	Vector v(3),v2(2);
	Vector ex(3);
	Vector p(3);
	float ux,uy,uz;
	float C,S,t;
	float g;
	
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

bool IntegratorTriforce::isInPositiveEpsilonRange(float v, float eps){
	if(eps-(v+THRESHOLD_NUMERICAL) <= 0) return true;
	else return false;
}


float IntegratorTriforce::V2phi(Vector &integrationOrigin, Vector cv, Vector &v){
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


fmat33 IntegratorTriforce::rotz(float theta){
	fmat33 m;
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


float IntegratorTriforce::logisticSmoother(float x){
	return 1.0/(1.0+exp(-((x*2*LOGISTIC_SMOOTHER_PARAMETER)-LOGISTIC_SMOOTHER_PARAMETER)));
}



float IntegratorTriforce::sech(float x){
	return 1.0/(0.5 * (exp(x)+exp(-x)));
}





float IntegratorTriforce::dlogisticSmoother(float x){
	return LOGISTIC_SMOOTHER_PARAMETER / (1+cosh(LOGISTIC_SMOOTHER_PARAMETER - 2*LOGISTIC_SMOOTHER_PARAMETER*x));

}



Vector IntegratorTriforce::areaSmoother(Vector &x, float area, float radius){
	return x*area/(2*M_PI*radius*radius);
}



Integral IntegratorTriforce::integrateTriangle(int l, SASASegment &x, Vector integrationOrigin){
	Rotation psi;
	Rotation lambda;
	Rotation PHIij;
	Rotation PHIjk;
	float area;
	Vector M;
	Integral a;
	Vector Tij(4), Tjk(4);
	Vector Tij2(4), Tjk2(4);
	Vector force_i(3), force_j(3), force_k(3), force_l(3), force_t(3);
	float q, q0, q1;
	CircularInterfaceForm form;
	CircularInterfaceForm formi;
	CircularInterfaceForm formj;
	CircularInterfaceForm formk;
	float s_convex, s_direction;
	
	
	area = 0;
	
	psi = x.psi;
	lambda = x.lambda;
	form = x.form1;
	
	formi = x.form0;
	formj = x.form1;
	formk = x.form2;
	

	
	q0 = 1;
	if(formi != CONVEX) q0*=-1;
	if(formj != CONVEX) q0*=-1;
	
	q1 = 1;
	if(formj != CONVEX) q1*=-1;
	if(formk != CONVEX) q1*=-1;
	
	
	
	
	
	PHIij = x.rotation0;
	PHIjk = x.rotation1;
	
	
		int i0,i1,i2;
	i0=749;
	i1=748;
	i2=746;

	
	
	s_convex=1;
	if(form!=CONVEX){
		PHIij.rotation *=-1;
		PHIjk.rotation *=-1;
		s_convex=-1;
		
	}
	
		
	M = lookUp(M_PI, psi.rotation, lambda.rotation, form);
	
	
	Tij = lookUp(PHIij.rotation, psi.rotation, lambda.rotation, form);
	Tjk = lookUp(PHIjk.rotation, psi.rotation, lambda.rotation, form);
	

	
	/*
	if(x1.form == CONVEX && psi.rotation < lambda.rotation){
		Tij*=-1;
		Tjk*=-1;
	}
	*/
	
	//if(l==0) fprintf(stderr,"T: %f %f %f %f (%f %f)\n",PHIij.rotation, PHIjk.rotation, psi.rotation, lambda.rotation, Tij(0), Tjk(0));
	
	if(PHIjk.rotation >= PHIij.rotation){
		s_direction=1;
	}
	else{
		s_direction=-1;
	}
	
	
	area = s_direction*Tjk(0) - s_direction*Tij(0);
	
	
	//if(l==0) fprintf(stderr,"T1: %f\n",area);
	
	
	q= s_convex * s_direction;
	
	
	force_i = -q*(Tij(1) * PHIij.drotation_dxj);
	force_j = q*(Tjk(1) * PHIjk.drotation_dxi + Tjk(2) * psi.drotation_dxi + Tjk(3) * lambda.drotation_dxi) - q*(Tij(1) * PHIij.drotation_dxi + Tij(2) * psi.drotation_dxi + Tij(3) * lambda.drotation_dxi);
	force_k = q*(Tjk(1) * PHIjk.drotation_dxj);
	force_l = q*((Tjk(1) * PHIjk.drotation_dxl + Tjk(2) * psi.drotation_dxl + Tjk(3) * lambda.drotation_dxl) - (Tij(1) * PHIij.drotation_dxl + Tij(2) * psi.drotation_dxl + Tij(3) * lambda.drotation_dxl));

	
	//printf("SIZES: %d %d %d %d %d %d\n",Tjk.size(), Tij.size(), PHIjk.drotation_dxt.size(), psi.drotation_dxt.size(), PHIij.drotation_dxt.size(), psi.drotation_dxt.size());
	
	force_t = q*((Tjk(1) * PHIjk.drotation_dxt + Tjk(2) * psi.drotation_dxt) - (Tij(1) * PHIij.drotation_dxt + Tij(2) * psi.drotation_dxt));
	
	//printf("Tij: %f, Tjk: %f, PSI: (%f %f %f)\n",Tij(2),Tjk(2),psi.drotation_dxj(0),psi.drotation_dxj(1),psi.drotation_dxj(2));
	//printf("FORCE_T: (%f %f %f)\n",force_t(0),force_t(1),force_t(2));
	
	
  	if(//(l==i0 || x.index0==i0 || x.index1==i0 || x.index2==i0) ||
 		(l==i0 && ((x.index0==i1 && x.index1==i2) || (x.index0==i2 && x.index1==i1) || (x.index1==i1 && x.index2==i2) || (x.index1==i2 && x.index2==i1))) ||
 		(l==i1 && ((x.index0==i0 && x.index1==i2) || (x.index0==i2 && x.index1==i0) || (x.index1==i0 && x.index2==i2) || (x.index1==i2 && x.index2==i0))) ||
 		(l==i2 && ((x.index0==i1 && x.index1==i0) || (x.index0==i0 && x.index1==i1) || (x.index1==i1 && x.index2==i0) || (x.index1==i0 && x.index2==i1))))
 	{
// 		printf("ROT: %f %f\n",PHIjk.rotation,PHIij.rotation);
	}
 	
// 		printf("X: %d (%d %d %d)\n",l,x.index0,x.index1,x.index2);
// 	}
// 		printf("INDEX: %d %d %d\n",x.index0,x.index1,x.index2);
//   		printf("PHIij dxi : (%f %f %f) \t PHIij dxj : (%f %f %f) \t PHIij dxl : (%f %f %f)\n",PHIij.drotation_dxi(0),PHIij.drotation_dxi(1),PHIij.drotation_dxi(2),PHIij.drotation_dxj(0),PHIij.drotation_dxj(1),PHIij.drotation_dxj(2),PHIij.drotation_dxl(0),PHIij.drotation_dxl(1),PHIij.drotation_dxl(2));
//   		printf("PHIjk dxi : (%f %f %f) \t PHIjk dxj : (%f %f %f) \t PHIjk dxl : (%f %f %f)\n",PHIjk.drotation_dxi(0),PHIjk.drotation_dxi(1),PHIjk.drotation_dxi(2),PHIjk.drotation_dxj(0),PHIjk.drotation_dxj(1),PHIjk.drotation_dxj(2),PHIjk.drotation_dxl(0),PHIjk.drotation_dxl(1),PHIjk.drotation_dxl(2));
//   		printf("psi dxi   : (%f %f %f) \t psi dxj   : (%f %f %f) \t psi dxl   : (%f %f %f)\n",psi.drotation_dxi(0),psi.drotation_dxi(1),psi.drotation_dxi(2),psi.drotation_dxj(0),psi.drotation_dxj(1),psi.drotation_dxj(2),psi.drotation_dxl(0),psi.drotation_dxl(1),psi.drotation_dxl(2));
//   		printf("lambda dxi: (%f %f %f) \t lambda dxj: (%f %f %f) \t lambda dxl: (%f %f %f)\n",lambda.drotation_dxi(0),lambda.drotation_dxi(1),lambda.drotation_dxi(2),lambda.drotation_dxj(0),lambda.drotation_dxj(1),lambda.drotation_dxj(2),lambda.drotation_dxl(0),lambda.drotation_dxl(1),lambda.drotation_dxl(2));
//  	}
	
	/*
	int ti=12;
	if(x0.index0==ti || x0.index1==ti || x1.index0==ti || x1.index1==ti || l==ti){
	
		fprintf(stderr, "INDEXES: (%d) %d %d %d %d\n",l,x0.index0, x0.index1, x1.index0, x1.index1);
		fprintf(stderr, "BEFORE FORCES: %f %f %f\n",*(forces[ti][0]),*(forces[ti][1]),*(forces[ti][2]));
		fprintf(stderr, "Tij: %f %f %f %f\n",Tij(0),Tij(1),Tij(2),Tij(3));
		fprintf(stderr, "Tjk: %f %f %f %f\n",Tjk(0),Tjk(1),Tjk(2),Tjk(3));
		fprintf(stderr, "M: %f %f %f %f\n",M(0),M(1),M(2),M(3));
		
		fprintf(stderr, "PHIij dxi : (%f %f %f) \t PHIij dxj : (%f %f %f) \t PHIij dxl : (%f %f %f)\n",PHIij.drotation_dxi(0),PHIij.drotation_dxi(1),PHIij.drotation_dxi(2),PHIij.drotation_dxj(0),PHIij.drotation_dxj(1),PHIij.drotation_dxj(2),PHIij.drotation_dxl(0),PHIij.drotation_dxl(1),PHIij.drotation_dxl(2));
		fprintf(stderr, "PHIjk dxi : (%f %f %f) \t PHIjk dxj : (%f %f %f) \t PHIjk dxl : (%f %f %f)\n",PHIjk.drotation_dxi(0),PHIjk.drotation_dxi(1),PHIjk.drotation_dxi(2),PHIjk.drotation_dxj(0),PHIjk.drotation_dxj(1),PHIjk.drotation_dxj(2),PHIjk.drotation_dxl(0),PHIjk.drotation_dxl(1),PHIjk.drotation_dxl(2));
		fprintf(stderr, "psi dxi   : (%f %f %f) \t psi dxj   : (%f %f %f) \t psi dxl   : (%f %f %f)\n",psi.drotation_dxi(0),psi.drotation_dxi(1),psi.drotation_dxi(2),psi.drotation_dxj(0),psi.drotation_dxj(1),psi.drotation_dxj(2),psi.drotation_dxl(0),psi.drotation_dxl(1),psi.drotation_dxl(2));
		fprintf(stderr, "lambda dxi: (%f %f %f) \t lambda dxj: (%f %f %f) \t lambda dxl: (%f %f %f)\n",lambda.drotation_dxi(0),lambda.drotation_dxi(1),lambda.drotation_dxi(2),lambda.drotation_dxj(0),lambda.drotation_dxj(1),lambda.drotation_dxj(2),lambda.drotation_dxl(0),lambda.drotation_dxl(1),lambda.drotation_dxl(2));
	}
	*/
	q=s_convex;
	
	if(PHIjk.rotation < PHIij.rotation){
		area = 2*M(0) - area;
		force_i = -force_i;
		force_j = q* (2*(M(2)*psi.drotation_dxi + M(3)*lambda.drotation_dxi)) - force_j;
		force_k = -force_k;
		force_l = q*(2*(M(2)*psi.drotation_dxl + M(3)*lambda.drotation_dxl)) - force_l;
		force_t = q*(2*(M(2)*psi.drotation_dxt)) - force_t;
		//force_t =  -force_t;
	}
	
	//if(l==0) fprintf(stderr,"T2: %f\n",area);

	
	if(form==CONVEX) s_convex=1;
	else s_convex=-1;
	
	q = s_convex;
	
	
	a.integral=q * area;
	a.force_i = q* force_i;
	a.force_j = q* force_j;
	a.force_k = q* force_k;
	a.force_l = q* force_l;
	//a.force_t = Vector(3).zeros();
	a.force_t = q* force_t;
	
	
// 	if(	(l==i0 && ((x.index0==i1 && x.index1==i2) || (x.index0==i2 && x.index1==i1) || (x.index1==i1 && x.index2==i2) || (x.index1==i2 && x.index2==i1))) ||
// 		(l==i1 && ((x.index0==i0 && x.index1==i2) || (x.index0==i2 && x.index1==i0) || (x.index1==i0 && x.index2==i2) || (x.index1==i2 && x.index2==i0))) ||
// 		(l==i2 && ((x.index0==i1 && x.index1==i0) || (x.index0==i0 && x.index1==i1) || (x.index1==i1 && x.index2==i0) || (x.index1==i0 && x.index2==i1))))
// 	
// 	{
// 		printf("forcei: %f %f %f\n",force_i(0)*tradius*tradius,force_i(1)*tradius*tradius,force_i(2)*tradius*tradius);
// 		printf("forcej: %f %f %f\n",force_j(0)*tradius*tradius,force_j(1)*tradius*tradius,force_j(2)*tradius*tradius);
// 		printf("forcek: %f %f %f\n",force_k(0)*tradius*tradius,force_k(1)*tradius*tradius,force_k(2)*tradius*tradius);
// 		printf("forcel: %f %f %f\n",force_l(0)*tradius*tradius,force_l(1)*tradius*tradius,force_l(2)*tradius*tradius);
// 		printf("---\n");
// 	}
// 	
	
	
	/*
	
	//if(l==0) fprintf(stderr,"T3: %f\n",a.area);
	
	
	//here, we apply a logistic function to smooth the derivative for small lambdas

	int index_j, index_l, index_i, index_k;
	Vector v_j,mu_j,mu_l;
	float d_j,r_j,r_l, t_j, dd;
	
	index_i=x0.index0;
	index_j=x0.index1;
	index_k=x1.index1;
	index_l=l;
	
	if(index_j>=0){
	
		v_j = atoms[index_j] - atoms[index_l];
		d_j = norm(v_j,2);
		mu_j = v_j/d_j;
		r_j = radii[index_j];
		
		mu_l = -mu_j;
		r_l = radii[index_l];
		
		
		t_j = 1.0-d_j/(r_j+r_l);
		
		
		if(t_j <= LOGISTIC_LIMIT){
			
			area=a.area;
			
			
			
			ls = logisticSmoother(t_j/LOGISTIC_LIMIT);
			dls = dlogisticSmoother(t_j/LOGISTIC_LIMIT);
			
			dd = -LOGISTIC_LIMIT/(r_j+r_l);
			
			//fprintf(stderr,"TJ: %f %f %f\n",t_j,ls, dls);
			
			
			
			a.area = a.area * ls;
			if(index_i==index_j)
				a.force_i = a.force_i * ls + area * dls * dd * mu_j;
			else
				a.force_i = a.force_i * ls;
			a.force_j = a.force_j * ls + area * dls * dd * mu_j;
			if(index_k==index_j)
				a.force_k = a.force_k * ls + area * dls * dd * mu_j;
			else
				a.force_k = a.force_k * ls;
			a.force_l = a.force_l * ls + area * dls * dd * mu_l;
		}
	}
	
	*/
	
	return a;

}










float IntegratorTriforce::angle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}

float IntegratorTriforce::complAngle(Vector &a, Vector &b){
	return asin(norm_dot(a,b));
}

int IntegratorTriforce::sgn(float d){
	if(d>=0) return 1;
	else return -1;
}




float IntegratorTriforce::csc(float a){
	return 1.0/sin(a);
}




Benchmark IntegratorTriforce::getBenchmark(){
	return benchmark;
}


void IntegratorTriforce::outputPatches(FILE* outputfile, Molecule *m, Tessellation *tessellation){
	SASAs sasas;
	Vector integrationOrigin;
	vector<Vector> atoms;
	vector<float> radii;	
	SASASegmentList::iterator it;
	SASASegment x;
	float radius;
	
	
	radii = molecule->fetchRadii();
	atoms = molecule->fetchCoordinates();
	
	
	sasas = tessellation->sasas();
	
	fprintf(outputfile,"normal_for_circular_interface\tPHI0\tPHI1\tpsi\tlambda\tradius\themisphere\tform\tindexl\tindex0\tindex1\tindex2\n");
	for(unsigned int i=0;i<sasas.size();++i){
		radius = radii[i];
		for(it = sasas[i].begin(); it!=sasas[i].end(); ++it){
			x=*it;
		
			fprintf(outputfile,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\n",
								x.normalForCircularInterface(0),
								x.normalForCircularInterface(1),
								x.normalForCircularInterface(2),
								x.rotation0.rotation,
								x.rotation1.rotation,
								x.psi.rotation,
								x.lambda.rotation,
								radii[i],
								x.hemisphere,
								x.form1,
								i,
								x.index0,
								x.index1,
								x.index2
								);
		}
	}
}
