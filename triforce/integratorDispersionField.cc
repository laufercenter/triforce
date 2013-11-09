#include "integratorDispersionField.h"

#include <algorithm>
#include <string>
#include <limits>


using namespace std;
using namespace arma;
using namespace boost;

IntegratorDispersionField::IntegratorDispersionField(vector<Interpolation*> data){
	unsigned int j;
	j=0;
	
	this->dataSolvationFreeEnergy = data[j];
	++j;
	
	for(unsigned int i=0; i<4; ++i){
		this->data[i]=data[j];
		++j;
	}
	
	for(unsigned int parameter=0; parameter<2; ++parameter)
		for(unsigned int species=0; species<N_SPECIES; ++species)
			for(unsigned int hemisphere; hemisphere<2; ++hemisphere)
				for(unsigned int i=0; i<DIM_DISPERSIONFIELD; ++i){
					this->dataDispersion[parameter][species][hemisphere][i]=dataDispersion[j];
					++j;
				}
				
			
}


void IntegratorDispersionField::convertSegment(SASASegment &reference, SASASegment &container, Vector integrationAxis){
	//get the two integration plane normals (old and new)
	n_old = cross(reference.normalForCircularInterface, reference.tessellationAxis);
	n_new = cross(reference.normalForCircularInterface, integrationAxis);
	
	a = Tessellation::getAngle(n_old, n_new);
	
	container.PHIij.rotation += a;
	container.PHIjk.rotation += a;
	
	container.psi.rotation = Tessellation::getAngle(integrationAxis, reference.normalForCircularInterface);
	
	container.tessellationAxis = integrationAxis;
}


void IntegratorDispersionField::rotateSASA(SASASegmentList &reference, SASASegmentList &container, Vector &source, Vector &atom, float &d){
	Vector integrationAxis;
	Vector v;
	
	v=atom-source;
	d = norm(v,2);
		
	
	integrationAxis=-(v/d);
	
	for(unsigned int i=0; i<reference.size(); ++i){
		convertSegment(reference[i],container[i],integrationAxis);
	}
}


float IntegratorDispersionField::integrate(Molecule *m, Tessellation *tessellation){
	SASAs sasas;
	Vector integrationOrigin;
	float radius;
	Dispersion dispersion;
	vector<int> neighbourlist;
	SASASegmentList	rotatedSASA;
	set exposedBoundary;
	float d;
	this->molecule = m;
	this->tessellation = tessellation;
	float frontArea, backArea;
	float area,a,totalArea;
	DispersionElement M,M0,M1;
	Vector mLJP(2);
	
	radii = molecule->fetchRadii();
	atoms = molecule->fetchCoordinates();
	forces = molecule->fetchForcePointers();
	areas = molecule->fetchAreaPointers();
	species = molecule->fetchSpecies();
	nonpolar = molecule->fetchNonPolarPointers();
	epsilons = molecule->fetchEpsilons();
	sigmas = molecule->fetchSigmas();
	
	
	clearForces();
	
	
	
	
	sasas = tessellation->sasas();
	dispersion = 0;
	
	//iterate over all atoms
	
	benchmark.start(string("integrating"));
	
	for(unsigned int i=0;i<sasas.size();++i){
		radius = radii[i];
		
		//first we calculate areas
		area = integrateSASA(i, sasas[i], radius, frontArea, backArea);
		a = radius*radius * area;
		
		*(areas[i]) = a;
		totalArea += a;
		
		
		
		
		rotatedSASA = sasas[i];
		exposedBoundary.clear();
		
		for(unsigned int j=0; j<rotatedSASA[j]; ++j){
			exposedBoundary.insert(rotatedSASA[j].index1);
		}
		
		neighbourlist = molecule.getNeighborListFor(i);
		//we set the eps and sig estimates to the source atoms eps and sig values
		est_eps=epsilons[i];
		est_sig=sigmas[i];
		for(unsigned int j=0; j<neighbourlist.size(); ++j){
			M0 = lookUp(M_PI, 0, M_PI*0.5, d, eps, sig, CONVEX, FRONTHEMISPHERE, s);
			M1 = lookUp(M_PI, 0, M_PI*0.5, d, eps, sig, CONVEX, BACKHEMISPHERE, s);
			M.eps = M0.eps + M1.eps;
			M.sig = M0.sig + M1.sig;
			
			
			//if it is an exposed boundary atom, we treat it explicitely
			if(exposedBoundary.find(j)!=exposedBoundary.end()){
			
				rotateSASA(sasas[i],rotatedSASA, atoms[i], atoms[neighbourlist[j]], d);
				dispersion = integrateSASA(i, rotatedSASA, radius, d, est_eps, est_sig, M, frontArea, backArea);
				//average and get new estimates
				est_eps = dispersion.eps/=area;
				est_sig = dispersion.sig/=area;
				
			}
			//otherwise, we treat it implicitely
			else{
				//we take only the front hemisphere, and as such half the surface area of a unit sphere
				est_eps = M0.eps / (2*M_PI);
				est_sig = M0.sig / (2*M_PI);
			}
			
			
			//if(a>0) pushForces();
			//purgeForces();
		}
		
		//back out LJ parameter
		mLJP(0) = est_eps*est_eps;
		mLJP(1) = es_sig*2.0 - 0.82;
		
		
		//look up
		*(nonpolar[i]) = dataSolvationFreeEnergy->interpolate(mLJP);
		
		
	}
	
	benchmark.stop();
	
	return dispersion;
	
	
}



Dispersion IntegratorDispersionField::integrateSASA(int l, SASASegmentList &sasa, float radius, float d, float eps, float sig, float area, unsigned int s, DispersionElement M, float frontArea, float backArea){
	SASASegmentList::iterator it;
	SASASegment x;
	Dispersion dispersion0,dispersion1,dispersionAccumulated;
	DispersionIntegral dispersionIntegral;
	Integral areaIntegral;
	float r_square;
	float actphi;
	Vector c(2),c2(2);
	float dispersionEps;
	float dispersionSig;
	Vector p0, p1, p2, p01, p12, n01(2);
	
	M0 = lookUp(M_PI, 0, M_PI*0.5, d, eps, sig, CONVEX, FRONTHEMISPHERE, s);
	M1 = lookUp(M_PI, 0, M_PI*0.5, d, eps, sig, CONVEX, BACKHEMISPHERE, s);
	
	
	r_square = radius*radius;
	
	dispersion0=0;
	dispersion1=0;
	for(it = sasa.begin(); it!=sasa.end(); ++it){
		x=*it;
		dispersionIntegral = integrateTriangle(l, x, x.tessellationAxis, d, eps, sig, s);
		
		dispersionEps = dispersionIntegral.eps.integral;
		dispersionSig = dispersionIntegral.sig.integral;
		
		
		/*
		addForce(x.index0, integral.force_i * r_square);
		addForce(x.index1, integral.force_j * r_square);
		addForce(x.index2, integral.force_k * r_square);
		addForce(l, integral.force_l * r_square);
		*/
		

		if(x.hemisphere == FRONTHEMISPHERE){
			dispersion0.eps += dispersionEps;
			dispersion0.sig += dispersionSig;
		}
		else{
			dispersion1.eps += dispersionEps;
			dispersion1.sig += dispersionSig;
		}
		
		
		
		
		
		
	}
	
	if(frontArea < 0 && abs(frontArea) < THRESHOLD_NEGATIVE){
	}
	if(frontArea < 0){
		dispersion0.eps=M.eps(0) + dispersion0.eps;
		dispersion0.sig=M.sig(0) + dispersion0.sig;
	}

	if(backArea < 0 && abs(backArea) < THRESHOLD_NEGATIVE){
	}
	if(backArea < 0){
		dispersion1.eps=M.eps(0) + dispersion1.eps;
		dispersion1.sig=M.sig(0) + dispersion1.sig;
	}
	
	
	
	dispersion.eps = dispersion0.eps + dispersion1.eps;
	dispersion.sig = dispersion0.sig + dispersion1.sig;
	
	
	/*
	//if we have information from the depthbuffer, check if the calculated size coarsely matches the one from the buffer
	//and if not, disregard the sasa (!)
	if(sasa.size()>0){
		if(sasa[0].depthBufferEstimatedArea>=0)
			if(dispersion> 2*sasa[0].depthBufferEstimatedArea*4*M_PI) dispersion=0;
	}
	*/
	
	
	return dispersion;
}







float IntegratorTriforce::integrateSASA(int l, SASASegmentList &sasa, float radius, float &frontArea, float &backArea){
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
		integral = integrateTriangle(l, x, x.tessellationAxis);
		
		a = integral.integral;
		
/*		
		addForce(x.index0, integral.force_i * r_square);
		addForce(x.index1, integral.force_j * r_square);
		addForce(x.index2, integral.force_k * r_square);
		addForce(l, integral.force_l * r_square);
*/		
		

		if(x.hemisphere == FRONTHEMISPHERE)
			area0 += a;
		else area1+=a;
		
	}
	
	frontArea=area0;
	backArea=area1;
	
	
	if(area0 < 0 && abs(area0) < THRESHOLD_NEGATIVE) area0 = 0;
	if(area0 < 0) area0 = 2*M_PI+area0;

	if(area1 < 0 && abs(area1) < THRESHOLD_NEGATIVE) area1 = 0;
	if(area1 < 0) area1 = 2*M_PI+area1;
	
	
	
	area = area0 + area1;
	
	//if we have information from the depthbuffer, check if the calculated size coarsely matches the one from the buffer
	//and if not, disregard the sasa (!)
	if(sasa.size()>0){
		if(sasa[0].depthBufferEstimatedArea>=0)
			if(area> 2*sasa[0].depthBufferEstimatedArea*4*M_PI) area=0;
	}
	
	
	
	return area;
}




DispersionElement IntegratorDispersionField::lookUp(float PHI, float psi, float lambda, float d, float eps, float sig, CircularInterfaceForm form, Hemisphere hemisphere, unsigned int species){
	Vector resEps(7);
	Vector resSig(7);
	float aPHI, apsi, alambda, ad, aeps, asig;
	int i;
	Vector x(6);
	Vector p(2);
	vector<VectorInt> sp;
	vector<float> w;
	DispersionElement res;
	
	sp.clear();
	w.clear();
	
	aPHI=abs(PHI);
	apsi = max(0.0f,psi);
	alambda= max(0.0f,lambda);
	ad= max(0.0f,d);
	aeps= max(0.0f,eps);
	asig= max(0.0f,sig);
	
	x(0)=aPHI;
	x(1)=apsi;
	x(2)=alambda;
	x(3)=ad;
	x(4)=aeps;
	x(5)=asig;
	
	
	for(i=0;i<DIM_DISPERSIONFIELD;++i){
		//printf("start lookup %d (%f %f %f)\n",i,PHI,psi,lambda);
		resEps(i) = -dataDispersion[EPSILON][species][hemisphere][i]->interpolate(x, sp, w);
		resSig(i) = -dataDispersion[SIGMA][species][hemisphere][i]->interpolate(x, sp, w);
		//printf("end lookup %d\n",i);
	}
	
	if(PHI <0){
		resEps=-resEps;
		resEps(1) = -resEps(1);

		resSig=-resSig;
		resSig(1) = -resSig(1);
		
	}
	
	
	if(form != CONVEX){
		resEps(2) = -resEps(2);
		resEps(3) = -resEps(3);
		resEps(4) = -resEps(4);
		resEps(5) = -resEps(5);
		resEps(6) = -resEps(6);
		
		resSig(2) = -resSig(2);
		resSig(3) = -resSig(3);
		resSig(4) = -resSig(4);
		resSig(5) = -resSig(5);
		resSig(6) = -resSig(6);
		
	}
	
	
	res.eps=resEps;
	res.sig=resSig;
	
	return res;
}
	
	
	
	
	

Dispersion IntegratorDispersionField::integrateTriangle(int l, SASASegment &x, Vector integrationOrigin, float d, float eps, float sig, unsigned int s){
	Rotation psi;
	Rotation lambda;
	Rotation PHIij;
	Rotation PHIjk;
	float e_eps,e_sig;
	Integral a;
	Vector force_i(6), force_j(6), force_k(6), force_l(6);
	float phi0,phi1;
	float q, q0, q1;
	CircularInterfaceForm form;
	CircularInterfaceForm formi;
	CircularInterfaceForm formj;
	CircularInterfaceForm formk;
	float s_convex, s_direction;
	unsigned int hemisphere;
	unsigned int s;
	DispersionElement Tij, Tjk, M;
	
	
	e_eps = 0;
	e_sig = 0;
	
	psi = x.psi;
	lambda = x.lambda;
	form = x.form1;
	
	formi = x.form0;
	formj = x.form1;
	formk = x.form2;
	
	hemisphere=x.hemisphere;
	

	
	q0 = 1;
	if(formi != CONVEX) q0*=-1;
	if(formj != CONVEX) q0*=-1;
	
	q1 = 1;
	if(formj != CONVEX) q1*=-1;
	if(formk != CONVEX) q1*=-1;
	
	
	
	
	
	PHIij = x.rotation0;
	PHIjk = x.rotation1;
	
	
	
	
	
	s_convex=1;
	if(form!=CONVEX){
		PHIij.rotation *=-1;
		PHIjk.rotation *=-1;
		s_convex=-1;
		
	}
	
		
	M = lookUp(M_PI, psi.rotation, lambda.rotation, d, eps, sig, form, hemisphere, s);
	
	
	Tij = lookUp(PHIij.rotation, psi.rotation, lambda.rotation, d, eps, sig, form, hemisphere, s);
	Tjk = lookUp(PHIjk.rotation, psi.rotation, lambda.rotation, d, eps, sig, form, hemisphere, s);
	
	
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
	
	
	dispersionEps = s_direction*Tjk.eps(0) - s_direction*Tij.eps(0);
	dispersionSig = s_direction*Tjk.sig(0) - s_direction*Tij.sig(0);
	
	
	//if(l==0) fprintf(stderr,"T1: %f\n",dispersion);
	
	
	q= s_convex * s_direction;
	
	/*
	force_i = -q*(Tij(1) * PHIij.drotation_dxj);
	force_j = q*(Tjk(1) * PHIjk.drotation_dxi + Tjk(2) * psi.drotation_dxi + Tjk(3) * lambda.drotation_dxi) - q*(Tij(1) * PHIij.drotation_dxi + Tij(2) * psi.drotation_dxi + Tij(3) * lambda.drotation_dxi);
	force_k = q*(Tjk(1) * PHIjk.drotation_dxj);
	force_l = q*((Tjk(1) * PHIjk.drotation_dxl + Tjk(2) * psi.drotation_dxl + Tjk(3) * lambda.drotation_dxl) - (Tij(1) * PHIij.drotation_dxl + Tij(2) * psi.drotation_dxl + Tij(3) * lambda.drotation_dxl));
	*/

	q=s_convex;
	
	if(PHIjk.rotation < PHIij.rotation){
		dispersionEps = 2*M.eps(0) - dispersionEps;
		dispersionSig = 2*M.sig(0) - dispersionSig;
		
		/*
		force_i = -force_i;
		force_j = q* (2*(M(2)*psi.drotation_dxi + M(3)*lambda.drotation_dxi)) - force_j;
		force_k = -force_k;
		force_l = q*(2*(M(2)*psi.drotation_dxl + M(3)*lambda.drotation_dxl)) - force_l;
		*/
	}
	

	
	if(form==CONVEX) s_convex=1;
	else s_convex=-1;
	
	q = s_convex;
	
	
	a.eps=q * dispersionEps;
	a.sig=q * dispersionSig;
	
	/*
	a.force_i = q* force_i;
	a.force_j = q* force_j;
	a.force_k = q* force_k;
	a.force_l = q* force_l;
	*/
	
	
	
	
	
	
	
	
	/*
	
	//if(l==0) fprintf(stderr,"T3: %f\n",a.dispersion);
	
	
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
			
			dispersion=a.dispersion;
			
			
			
			ls = logisticSmoother(t_j/LOGISTIC_LIMIT);
			dls = dlogisticSmoother(t_j/LOGISTIC_LIMIT);
			
			dd = -LOGISTIC_LIMIT/(r_j+r_l);
			
			//fprintf(stderr,"TJ: %f %f %f\n",t_j,ls, dls);
			
			
			
			a.dispersion = a.dispersion * ls;
			if(index_i==index_j)
				a.force_i = a.force_i * ls + dispersion * dls * dd * mu_j;
			else
				a.force_i = a.force_i * ls;
			a.force_j = a.force_j * ls + dispersion * dls * dd * mu_j;
			if(index_k==index_j)
				a.force_k = a.force_k * ls + dispersion * dls * dd * mu_j;
			else
				a.force_k = a.force_k * ls;
			a.force_l = a.force_l * ls + dispersion * dls * dd * mu_l;
		}
	}
	
	*/
	
	return a;

}






