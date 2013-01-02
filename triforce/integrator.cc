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
			//*(forces[i][j]) = 0;
		}
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
		a = integrateAtomicSASA(i, sasas[i]);
		
		printf("TOTAL SUBAREA: %f\n",a);
		
		a = radius*radius * a;
		
		*(areas[i]) = a;
		area += a;
		
		
	}
	
	return area;
	
	
}


double IntegratorTriforce::integrateAtomicSASA(int l, SASAsForAtom sasasForAtom){
	double radius;
	double area0, area1, a, area;
	
	radius = sasasForAtom.radius;
	area0 = 0;
	area1 = 0;
	for(int i=0;i<sasasForAtom.sasas.size();++i){
		//int i=1;
		
		a = integrateSASA(l, sasasForAtom.sasas[i]);
		
		printf("ACC SUB SUBAREA: %f\n",a);

		
		if(sasasForAtom.sasas[i].hemisphere == FRONTHEMISPHERE)
			area0 += a;
		else area1+=a;
		
	}
	
	printf("ACC AREA0: %f AREA1: %f\n",area0,area1);
	
	
	if(area0<0) area0 = -area0;
	else area0=2*M_PI - area0;

	if(area1<0) area1 = -area1;
	else area1=2*M_PI - area1;
	
	printf("ACC CORRECTED AREA0: %f AREA1: %f\n",area0,area1);
	
	
	area = area0 + area1;
	
	return area;
	
}


void IntegratorTriforce::addForce(int i, Vector force){
	int j;
	for(j=0; j<3; ++j){
		//*(forces[i][j]) += force(j);
	}
}

double IntegratorTriforce::integrateSASA(int l, SASA &sasa){
	SASANodeList::iterator it;
	SASANode x0, x1;
	double area=0;
	double totalAngle=0;
	Area integral;
	
	x0 = *(--sasa.sasa.end());
	for(it = sasa.sasa.begin(); it!=sasa.sasa.end(); ++it){
		x1 = *it;
		
		integral = integrateTriangle(x0, x1, sasa.tessellationOrigin);
		printf("SUBAREA: %f\n",integral.area);

		area += integral.area;
		addForce(x0.index0, integral.force_i);
		addForce(x0.index1, integral.force_j);
		addForce(x1.index1, integral.force_k);
		addForce(l, integral.force_l);
		
		
		
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
	
	
	res(0) = data[0]->interpolate(aPHI, psi, lambda);
	res(1) = 0;
	res(2) = 0;
	res(3) = 0;
	
	//for(i=0;i<4;++i)
		//res(i) = data[i]->interpolate(aPHI, psi, lambda);
	
	if(PHI <0) res=-1*res;
	
	return res;
}
	
	


Area IntegratorTriforce::integrateTriangle(SASANode &x0, SASANode &x1, Vector integrationOrigin){
	Rotation psi;
	Rotation lambda;
	Rotation PHIij;
	Rotation PHIjk;
	double area,area0,area1;
	double maxArea;
	Area a;
	Vector Tij(4), Tjk(4);
	Vector force_i(3), force_j(3), force_k(3), force_l(3);
	
	area = 0;
	
	
	psi = x1.psi;
	lambda = x1.lambda;
	
	
	
	PHIij = x0.rotation1;
	PHIjk = x1.rotation0;
	
		
	maxArea = data[0]->interpolate(M_PI, psi.rotation, lambda.rotation);	
	
	
	Tij = lookUp(PHIij.rotation, psi.rotation, lambda.rotation);
	Tjk = lookUp(PHIjk.rotation, psi.rotation, lambda.rotation);
	
	
	
	area = Tjk(0)-Tij(0);
	
	force_i = -( Tij(1) * PHIij.drotation_dxi + Tij(2) * psi.drotation_dxi + Tij(3) * lambda.drotation_dxi);
	force_j = Tjk(1) * PHIjk.drotation_dxi + Tjk(2) * psi.drotation_dxi + Tjk(3) * lambda.drotation_dxi - (Tij(1) * PHIij.drotation_dxj + Tij(2) * psi.drotation_dxj + Tij(3) * lambda.drotation_dxj);
	force_k = Tjk(1) * PHIjk.drotation_dxj + Tjk(2) * psi.drotation_dxj + Tjk(3) * lambda.drotation_dxj;
	force_l = Tjk(1) * PHIjk.drotation_dxl + Tjk(2) * psi.drotation_dxl + Tjk(3) * lambda.drotation_dxl - (Tij(1) * PHIij.drotation_dxl + Tij(2) * psi.drotation_dxl + Tij(3) * lambda.drotation_dxl);
	
	
	if(PHIjk.rotation < PHIij.rotation){
		area = 2*maxArea + area;
	}
	
	if(x1.form!=CONVEX){
		area=-1 * area;
		force_i = -1* force_i;
		force_j = -1* force_j;
		force_k = -1* force_k;
		force_l = -1* force_l;
	}
	
	a.area = area;
	a.force_i = -1* force_i;
	a.force_j = -1* force_j;
	a.force_k = -1* force_k;
	a.force_l = -1* force_l;
	

	
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


