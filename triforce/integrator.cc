#include "integrator.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;


IntegratorTriforce::IntegratorTriforce(){
	
}



IntegratorTriforce::IntegratorTriforce(Interpolation *dataConcave, Interpolation *forcesConcave0, Interpolation *forcesConcave1, Interpolation *forcesConcave2){
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
		printf("RADIUS[%d]: %f\n",i,radius);
	}
		
	
	
	sasas = tessellation->sasas();
	area = 0;
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		radius = sasas[i].radius;
		a = integrateAtomicSASA(i, sasas[i], radius);
		
		printf("TOTAL SUBAREA: %f\n",a);
		
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
		
		printf("ACC SUB SUBAREA: %f\n",a);

		
		if(sasasForAtom.sasas[i].hemisphere == FRONTHEMISPHERE)
			area0 += a;
		else area1+=a;
		
	}
	
	printf("ACC AREA0: %f AREA1: %f\n",area0,area1);
	
	
	if(area0<=0) area0 = -area0;
	else area0=2*M_PI - area0;

	if(area1<=0) area1 = -area1;
	else area1=2*M_PI - area1;
	
	printf("ACC CORRECTED AREA0: %f AREA1: %f\n",area0,area1);
	
	
	area = area0 + area1;
	
	return area;
	
}


void IntegratorTriforce::addForce(int i, Vector force){
	int j;
	printf("ADDING FORCE TO %d\n",i);
	if(i>=0){
		for(j=0; j<3; ++j){
			*(forces[i][j]) += force(j);
		}
	}
	else printf("DISREGARDING SPLITTER FORCE %f %f %f\n",force(0),force(1),force(2));
}

double IntegratorTriforce::integrateSASA(int l, SASA &sasa, double radius){
	SASANodeList::iterator it;
	SASANode x0, x1;
	double area=0;
	double totalAngle=0;
	Area integral;
	double r_square;
	
	r_square = radius*radius;
	
	x0 = *(--sasa.sasa.end());
	for(it = sasa.sasa.begin(); it!=sasa.sasa.end(); ++it){
		x1 = *it;
		
		integral = integrateTriangle(l, x0, x1, sasa.tessellationOrigin);
		printf("SUBAREA: %f\n",integral.area);

		area += integral.area;
		
		if(sasa.hemisphere==BACKHEMISPHERE){
			integral.force_i(0) *= -1;
			integral.force_j(0) *= -1;
			integral.force_k(0) *= -1;
			integral.force_l(0) *= -1;
		}
		
		addForce(x0.index0, integral.force_i * r_square);
		addForce(x0.index1, integral.force_j * r_square);
		addForce(x1.index1, integral.force_k * r_square);
		addForce(l, integral.force_l * r_square);
		
		
		
		x0=x1;
		
	}
	
	printf("+++++ SAA END +++++\n\n");
	
	return area;
}



Vector IntegratorTriforce::lookUp(double PHI, double psi, double lambda){
	Vector res(4);
	double aPHI;
	int i;
	
	aPHI=abs(PHI);
	
	
	
	for(i=0;i<4;++i)
		res(i) = data[i]->interpolate(aPHI, psi, lambda);
	
	if(PHI <0) res=-1*res;
	
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


Area IntegratorTriforce::integrateTriangle(int l, SASANode &x0, SASANode &x1, Vector integrationOrigin){
	Rotation psi;
	Rotation lambda;
	Rotation PHIij;
	Rotation PHIjk;
	double area,area0,area1;
	Vector M;
	Area a;
	Vector Tij(4), Tjk(4);
	Vector force_i(3), force_j(3), force_k(3), force_l(3);
	double totalAngle;
	double phi0,phi1;
	double q;
	CircularInterfaceForm form;
	area = 0;
	totalAngle=0;
	double s_convex, s_complementation, s_direction;
	
	psi = x1.psi;
	lambda = x1.lambda;
	form = x1.form;
	
	printf("NORMAL: %f %f %f\n",x1.normalForCircularInterface(0),x1.normalForCircularInterface(1),x1.normalForCircularInterface(2));
	
	
	PHIij = x0.rotation1;
	PHIjk = x1.rotation0;
	
	
	
	printf("ROTATIONS: %f %f %f %f\n",x0.rotation0.rotation,x0.rotation1.rotation,x1.rotation0.rotation,x1.rotation1.rotation);
	
	
	s_convex=1;
	if(form!=CONVEX){
		PHIij.rotation *=-1;
		PHIjk.rotation *=-1;
		s_convex=-1;
		
		printf("REVERSING ANGLE\n");
	}
	
		
	M = lookUp(M_PI, psi.rotation, lambda.rotation);
	//if(form!=CONVEX) maxArea*=-1;
	
	
	Tij = lookUp(PHIij.rotation, psi.rotation, lambda.rotation);
	Tjk = lookUp(PHIjk.rotation, psi.rotation, lambda.rotation);
	
	
	//finite differences
	//dPHIij_dxi
	Tessellation t(*this->molecule);
	Vector ai(3), aj(3), al(3), x(3), xp(3), xn(3);
	Vector nip,nin,nj;
	Vector njp,njn,ni;
	int indexi, indexl, indexj;
	double ri, rl, rj;
	Rotation lambdaip;
	Rotation lambdain;
	Rotation lambdaj;
	Rotation lambdajp;
	Rotation lambdajn;
	Rotation lambdai;
	double dip,din,dj;
	double djp,djn,di;
	Vector fd_dPHIij_dxi_in(3);
	Vector fd_dPHIij_dxi_out(3);
	Vector fd_dPHIij_dxj_in(3);
	Vector fd_dPHIij_dxj_out(3);
	PHIContainer PHIp, PHIn;
	Vector fd_dPHIij_dxl_out(3);
	Vector fd_dPHIij_dxl_in(3);
	
	printf("FORM: %d %d\n",x0.form, x1.form);
	
	
	//dPHIij_dxi
	indexj = x0.index0;
	indexi = x0.index1;
	indexl = l;
	if(indexi > 0 && indexj > 0){
		ai = atoms[indexi];
		aj = atoms[indexj];
		al = atoms[indexl];
		ri = radii[indexi];
		rj = radii[indexj];
		rl = radii[indexl];
		
		
		x=ai;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			nip = t.calculateInterfaceNormal(al, xp, dip);
			nin = t.calculateInterfaceNormal(al, xn, din);
			nj = t.calculateInterfaceNormal(al, aj, dj);
			
			lambdaip = t.calculateLambda(dip, rl, ri, nip);
			lambdain = t.calculateLambda(din, rl, ri, nin);
			lambdaj = t.calculateLambda(dj, rl, rj, nj);
			
			PHIp = t.calculatePHI(integrationOrigin, nip, nj, lambdaip, lambdaj);
			PHIn = t.calculatePHI(integrationOrigin, nin, nj, lambdain, lambdaj);
			
			fd_dPHIij_dxi_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
			fd_dPHIij_dxi_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
		}
		printf("DRV: dPHIij_dxi: %f %f %f fd_dPHIij_dxi_in: %f %f %f fd_dPHIij_dxi_out: %f %f %f\n",PHIij.drotation_dxi(0),PHIij.drotation_dxi(1),PHIij.drotation_dxi(2), fd_dPHIij_dxi_in(0), fd_dPHIij_dxi_in(1), fd_dPHIij_dxi_in(2), fd_dPHIij_dxi_out(0), fd_dPHIij_dxi_out(1), fd_dPHIij_dxi_out(2));
	}
	
	//dPHIij_dxj
	indexj = x0.index0;
	indexi = x0.index1;
	indexl = l;
	if(indexi > 0 && indexj > 0){
		ai = atoms[indexi];
		aj = atoms[indexj];
		al = atoms[indexl];
		ri = radii[indexi];
		rj = radii[indexj];
		rl = radii[indexl];
		
		
		x=aj;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			njp = t.calculateInterfaceNormal(al, xp, djp);
			njn = t.calculateInterfaceNormal(al, xn, djn);
			ni = t.calculateInterfaceNormal(al, ai, di);
			
			lambdajp = t.calculateLambda(djp, rl, rj, njp);
			lambdajn = t.calculateLambda(djn, rl, rj, njn);
			lambdai = t.calculateLambda(di, rl, ri, ni);
			
			PHIp = t.calculatePHI(integrationOrigin, ni, njp, lambdai, lambdajp);
			PHIn = t.calculatePHI(integrationOrigin, ni, njn, lambdai, lambdajn);
			
			fd_dPHIij_dxj_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
			fd_dPHIij_dxj_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
		}
		printf("DRV: dPHIij_dxj: %f %f %f fd_dPHIij_dxj_in: %f %f %f fd_dPHIij_dxj_out: %f %f %f\n",PHIij.drotation_dxj(0),PHIij.drotation_dxj(1),PHIij.drotation_dxj(2), fd_dPHIij_dxj_in(0), fd_dPHIij_dxj_in(1), fd_dPHIij_dxj_in(2), fd_dPHIij_dxj_out(0), fd_dPHIij_dxj_out(1), fd_dPHIij_dxj_out(2));
	}
	
	

	
	
	//dPHIij_dxl
	indexj = x0.index0;
	indexi = x0.index1;
	indexl = l;
	if(indexi > 0 && indexj > 0){
		ai = atoms[indexi];
		aj = atoms[indexj];
		al = atoms[indexl];
		ri = radii[indexi];
		rj = radii[indexj];
		rl = radii[indexl];
		
		
		x=al;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			nip = t.calculateInterfaceNormal(xp, ai, dip);
			nin = t.calculateInterfaceNormal(xn, ai, din);
			njp = t.calculateInterfaceNormal(xp, aj, djp);
			njn = t.calculateInterfaceNormal(xn, aj, djn);
			
			lambdaip = t.calculateLambda(dip, rl, ri, nip);
			lambdain = t.calculateLambda(din, rl, ri, nin);
			lambdajp = t.calculateLambda(djp, rl, rj, njp);
			lambdajn = t.calculateLambda(djn, rl, rj, njn);
			
			PHIp = t.calculatePHI(integrationOrigin, nip, njp, lambdaip, lambdajp);
			PHIn = t.calculatePHI(integrationOrigin, nin, njn, lambdain, lambdajn);
			
			fd_dPHIij_dxl_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
			fd_dPHIij_dxl_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
		}
		printf("DRV: dPHIij_dxl: %f %f %f fd_dPHIij_dxl_in: %f %f %f fd_dPHIij_dxl_out: %f %f %f\n",PHIij.drotation_dxl(0),PHIij.drotation_dxl(1),PHIij.drotation_dxl(2), fd_dPHIij_dxl_in(0), fd_dPHIij_dxl_in(1), fd_dPHIij_dxl_in(2), fd_dPHIij_dxl_out(0), fd_dPHIij_dxl_out(1), fd_dPHIij_dxl_out(2));
	}
		
		
	
	
	phi0 = PHI2phi2(integrationOrigin, PHIij.rotation,psi.rotation,lambda.rotation);
	phi1 = PHI2phi2(integrationOrigin, PHIjk.rotation,psi.rotation,lambda.rotation);
	
	
	if(PHIij.rotation<0){
		phi0=-phi0;
	}
	if(PHIjk.rotation<0) {
		phi1=-phi1;
	}
	if(x1.form!=CONVEX)
		totalAngle -= phi1 - phi0;
	else
		totalAngle += phi1 - phi0;
	
	printf("ID: (%d-%d) -> (%d-%d)\n",x0.id0, x0.id1,x1.id0,x1.id1);
	printf("CIRC: %d\n",x1.id0);
	printf("TOTAL ANGLE: %f phi0:%f (PHI0: %f) phi1:%f (PHI1: %f)\n",totalAngle,phi0, PHIij.rotation, phi1, PHIjk.rotation);
	printf("Tij: %f, Tjk: %f\n",Tij(0),Tjk(0));
	printf("MAXAREA: %f\n",M(0));
	
	
	if(PHIjk.rotation >= PHIij.rotation) s_direction=1;
	else s_direction=-1;
	
	area = s_direction*Tjk(0) - s_direction*Tij(0);
	
	q= s_convex * s_direction;
	
	force_i = -q*(Tij(1) * PHIij.drotation_dxj);
	force_j = q*(Tjk(1) * PHIjk.drotation_dxi + Tjk(2) * psi.drotation_dxi + Tjk(3) * lambda.drotation_dxi) - q*(Tij(1) * PHIij.drotation_dxi + Tij(2) * psi.drotation_dxi + Tij(3) * lambda.drotation_dxi);
	force_k = q*(Tjk(1) * PHIjk.drotation_dxj);
	force_l = q*(Tjk(1) * PHIjk.drotation_dxl + Tjk(2) * psi.drotation_dxl + Tjk(3) * lambda.drotation_dxl) - q*(Tij(1) * PHIij.drotation_dxl + Tij(2) * psi.drotation_dxl + Tij(3) * lambda.drotation_dxl);
	
	
	//printf("PHIij dxi : (%f %f %f) \t PHIij dxj : (%f %f %f) \t PHIij dxl : (%f %f %f)\n",PHIij.drotation_dxi(0),PHIij.drotation_dxi(1),PHIij.drotation_dxi(2),PHIij.drotation_dxj(0),PHIij.drotation_dxj(1),PHIij.drotation_dxj(2),PHIij.drotation_dxl(0),PHIij.drotation_dxl(1),PHIij.drotation_dxl(2));
	//printf("PHIjk dxi : (%f %f %f) \t PHIjk dxj : (%f %f %f) \t PHIjk dxl : (%f %f %f)\n",PHIjk.drotation_dxi(0),PHIjk.drotation_dxi(1),PHIjk.drotation_dxi(2),PHIjk.drotation_dxj(0),PHIjk.drotation_dxj(1),PHIjk.drotation_dxj(2),PHIjk.drotation_dxl(0),PHIjk.drotation_dxl(1),PHIjk.drotation_dxl(2));
	//printf("psi dxi   : (%f %f %f) \t psi dxj   : (%f %f %f) \t psi dxl   : (%f %f %f)\n",psi.drotation_dxi(0),psi.drotation_dxi(1),psi.drotation_dxi(2),psi.drotation_dxj(0),psi.drotation_dxj(1),psi.drotation_dxj(2),psi.drotation_dxl(0),psi.drotation_dxl(1),psi.drotation_dxl(2));
	//printf("lambda dxi: (%f %f %f) \t lambda dxj: (%f %f %f) \t lambda dxl: (%f %f %f)\n",lambda.drotation_dxi(0),lambda.drotation_dxi(1),lambda.drotation_dxi(2),lambda.drotation_dxj(0),lambda.drotation_dxj(1),lambda.drotation_dxj(2),lambda.drotation_dxl(0),lambda.drotation_dxl(1),lambda.drotation_dxl(2));
	
	
	
	
	
	
	
	
	Vector fd_dA_dxi(3);
	Vector nk;
	Rotation lambdak;
	Rotation psij;
	double dk;
	double Ap,An;
	PHIContainer PHI2;
	int indexk;
	double rk;
	Vector ak;
	
	
	
	//dAijk_dxi
	indexi = x0.index0;
	indexj = x0.index1;
	indexk = x1.index1;
	indexl = l;
	if(indexi >= 0 && indexj >= 0 && indexk >= 0){
		ai = atoms[indexi];
		aj = atoms[indexj];
		ak = atoms[indexk];
		al = atoms[indexl];
		ri = radii[indexi];
		rj = radii[indexj];
		rk = radii[indexk];
		rl = radii[indexl];
		
		
		x=ai;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			nip = t.calculateInterfaceNormal(al, xp, dip);
			nin = t.calculateInterfaceNormal(al, xn, din);
			nj = t.calculateInterfaceNormal(al, aj, dj);
			nk = t.calculateInterfaceNormal(al, ak, dk);
			
			lambdaip = t.calculateLambda(dip, rl, ri, nip);
			lambdain = t.calculateLambda(din, rl, ri, nin);
			lambdaj = t.calculateLambda(dj, rl, rj, nj);
			lambdak = t.calculateLambda(dk, rl, rk, nk);
			
			psij = t.calculatePsi(integrationOrigin, nj);
			
			PHIp = t.calculatePHI(integrationOrigin, nj, nip, lambdaj, lambdaip);
			PHIn = t.calculatePHI(integrationOrigin, nj, nin, lambdaj, lambdain);
			
			PHI2 = t.calculatePHI(integrationOrigin, nj, nk, lambdaj, lambdak);
			
			Ap = s_direction * (lookUp(PHI2.out.rotation,psij.rotation,lambdaj.rotation)(0) - lookUp(PHIp.in.rotation, psij.rotation, lambdaj.rotation)(0));
			An = s_direction * (lookUp(PHI2.out.rotation,psij.rotation,lambdaj.rotation)(0) - lookUp(PHIn.in.rotation, psij.rotation, lambdaj.rotation)(0));
			
			fd_dA_dxi(j) = (Ap-An) / (2*FD);
			
		}
		printf("DRV: dA_dxi: %f %f %f fd_dA_dxi: %f %f %f\n",force_i(0),force_i(1),force_i(2),fd_dA_dxi(0),fd_dA_dxi(1),fd_dA_dxi(2));
		
	}	
	
	
	
	
	Rotation psijp, psijn;
	PHIContainer PHI2p, PHI2n;
	Vector fd_dA_dxj(3);
	
	//dAijk_dxj
	indexi = x0.index0;
	indexj = x0.index1;
	indexk = x1.index1;
	indexl = l;
	if(indexi >= 0 && indexj >= 0 && indexk >= 0){
		ai = atoms[indexi];
		aj = atoms[indexj];
		ak = atoms[indexk];
		al = atoms[indexl];
		ri = radii[indexi];
		rj = radii[indexj];
		rk = radii[indexk];
		rl = radii[indexl];
		
		
		x=aj;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			njp = t.calculateInterfaceNormal(al, xp, djp);
			njn = t.calculateInterfaceNormal(al, xn, djn);
			ni = t.calculateInterfaceNormal(al, ai, di);
			nk = t.calculateInterfaceNormal(al, ak, dk);
			
			lambdajp = t.calculateLambda(djp, rl, rj, njp);
			lambdajn = t.calculateLambda(djn, rl, rj, njn);
			lambdai = t.calculateLambda(di, rl, ri, ni);
			lambdak = t.calculateLambda(dk, rl, rk, nk);
			
			psijp = t.calculatePsi(integrationOrigin, njp);
			psijn = t.calculatePsi(integrationOrigin, njn);
			
			PHIp = t.calculatePHI(integrationOrigin, njp, ni, lambdajp, lambdai);
			PHIn = t.calculatePHI(integrationOrigin, njn, ni, lambdajn, lambdai);
			
			PHI2p = t.calculatePHI(integrationOrigin, njp, nk, lambdajp, lambdak);
			PHI2n = t.calculatePHI(integrationOrigin, njn, nk, lambdajn, lambdak);
			
			Ap = s_direction * (lookUp(PHI2.out.rotation,psijp.rotation,lambdajp.rotation)(0) - lookUp(PHIp.in.rotation, psijp.rotation, lambdajp.rotation)(0));
			An = s_direction * (lookUp(PHI2.out.rotation,psijn.rotation,lambdajn.rotation)(0) - lookUp(PHIn.in.rotation, psijn.rotation, lambdajn.rotation)(0));
			
			fd_dA_dxj(j) = (Ap-An) / (2*FD);
			
		}
		printf("DRV: dA_dxj: %f %f %f fd_dA_dxj: %f %f %f\n",force_j(0),force_j(1),force_j(2),fd_dA_dxj(0),fd_dA_dxj(1),fd_dA_dxj(2));
		
	}	
	
	
		
	
	
	
	
	
	
	
	
	
	if(PHIjk.rotation < PHIij.rotation){
		area = 2*M(0) - area;
		force_j = 2*(M(2)*psi.drotation_dxi + M(3)*lambda.drotation_dxi) - force_j;
		force_l = 2*(M(2)*psi.drotation_dxl + M(3)*lambda.drotation_dxl) - force_l;
		
		printf("COMPLEMENTING\n");
	}
	

	
	
	/*
	if(psi.rotation < lambda.rotation){
		if(x1.form==CONVEX){
			q=1;
			printf("psi<lambda convex\n");
		}
		else{
			q=-1;
			printf("psi<lambda concave\n");
		}
	}
	else{
		if(x1.form==CONVEX){
			q=1;
			printf("psi>=lambda convex\n");
		}
		else{
			q=-1;
			printf("psi>=lambda concave\n");
		}
	}
	*/
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


