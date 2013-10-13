#include "integratorGaussBonnet.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;




IntegratorGaussBonnet::IntegratorGaussBonnet(){
	
}





float IntegratorGaussBonnet::integrate(Molecule *molecule, Tessellation *tessellation){
	this->molecule = molecule;
	this->tessellation = tessellation;
	
	SASAsForMolecule sasas;
	float area;
	
	sasas = tessellation->sasas();
	
	area = 0;
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		area += integrateAtomicSASA(sasas[i]);
	}
	
	return area;
	
	
}

int IntegratorGaussBonnet::sgn(float d){
	if(d>=0) return 1;
	else return -1;
}


float IntegratorGaussBonnet::integrateAtomicSASA(SASAsForAtom sasasForAtom){
	int i;
	float radius;
	float area;
	
	radius = sasasForAtom.radius;
	area = 0;
	for(int i=0;i<sasasForAtom.sasas.size();++i){
		area += integrateSASA(sasasForAtom.sasas[i]);
		
	}
	
	return area;
	
}

float IntegratorGaussBonnet::integrateSASA(SASA &sasa){
	SASANodeList::iterator it;
	SASANode x0, x1;
	float area;
	
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



float IntegratorGaussBonnet::integrateArc(SASANode &x0, SASANode &x1){
	float lambda_k, lambda_j;
	Vector mu_k(3), mu_j(3);
	float a_k, a_j;
	float g_k, g_j;
	Vector n_kj(3), m_kl(3), m_jk(3);
	float OMEGA;
	float S;
	float PHI;
	float cosTHETA;
	float res;
	float dot1, dot2;
	
	//calculate lambda
	lambda_k = x1.lambda.rotation;
	lambda_j = x0.lambda.rotation;
	mu_k=x1.normalForCircularInterface;
	mu_j=x0.normalForCircularInterface;
	
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
	
	
	
	
