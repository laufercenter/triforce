#include "integratorGaussBonnet.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;




IntegratorGaussBonnet::IntegratorGaussBonnet(){
	
}





double IntegratorGaussBonnet::integrate(Molecule *molecule, Tessellation *tessellation){
	this->molecule = molecule;
	this->tessellation = tessellation;
	
	SASAsForMolecule sasas;
	double area;
	
	sasas = tessellation->sasas();
	
	area = 0;
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		area += integrateAtomicSASA(sasas[i]);
	}
	
	return area;
	
	
}

int IntegratorGaussBonnet::sgn(double d){
	if(d>=0) return 1;
	else return -1;
}


double IntegratorGaussBonnet::integrateAtomicSASA(SASAsForAtom sasasForAtom){
	int i;
	double radius;
	double area;
	
	radius = sasasForAtom.radius;
	area = 0;
	for(int i=0;i<sasasForAtom.sasas.size();++i){
		area += integrateSASA(sasasForAtom.sasas[i]);
		
	}
	
	return area;
	
}

double IntegratorGaussBonnet::integrateSASA(SASA &sasa){
	SASANodeList::iterator it;
	SASANode x0, x1;
	double area;
	
	x0 = *(--sasa.sasa.end());
	area = 0;
	for(it = sasa.sasa.begin(); it!=sasa.sasa.end(); ++it){
		x1 = *it;
		
		area += integrateArc(x0, x1);
		x0=x1;
		
	}
	printf("GBONNET: %f\n",area);
	area = 2*M_PI+area;
	
	
	return area;
}



double IntegratorGaussBonnet::integrateArc(SASANode &x0, SASANode &x1){
	double lambda_k, lambda_j;
	Vector mu_k(3), mu_j(3);
	double a_k, a_j;
	double g_k, g_j;
	Vector n_kj(3), m_kl(3), m_jk(3);
	double OMEGA;
	double S;
	double PHI;
	double cosTHETA;
	double res;
	double dot1, dot2;
	
	//calculate lambda
	lambda_k = x1.lambda;
	lambda_j = x0.lambda;
	mu_k=x1.normalForCircularRegion;
	mu_j=x0.normalForCircularRegion;
	
	//extract g and a from lambda
	a_k = sin(lambda_k);
	a_j = sin(lambda_j);
	
	g_k = 1-cos(lambda_k);
	g_j = 1-cos(lambda_j);
	
	printf("F: %f, %f, %f, %f\n",a_k,a_j,g_k,g_j);
	printf("lambdas: %f, %f\n",lambda_k,lambda_j);
	
	

	n_kj = cross(mu_k, x0.vector);
	n_kj /= norm(n_kj,2);

	m_kl = cross(mu_k, x1.vector);
	m_kl /= norm(m_kl,2);
	
	m_jk = cross(mu_j, x0.vector);
	m_jk /= norm(m_jk,2);
	
	dot1 = dot(n_kj ,m_jk);
	dot2 = dot(n_kj,m_kl);
	OMEGA = -acos(dot(n_kj ,m_jk));
	S = sgn(dot(mu_k,cross(n_kj,m_kl)));
	PHI = (1-S)*M_PI + S*acos(dot(n_kj,m_kl));
	
	cosTHETA = g_k;
	
	printf("V: %f, %f, %f, %f, %f, %f\n",OMEGA,S,PHI,cosTHETA, dot1, dot2);

	
	res = OMEGA + PHI * cosTHETA;
	printf("SUBTOTAL: %f\n",res);
	return res;
}
	
	
	
	
