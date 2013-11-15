#include "integratorNumerical.h"

#include <algorithm>
#include <string>
#include <limits>


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>



using namespace std;
using namespace arma;
using namespace boost;



IntegratorNumerical::IntegratorNumerical(){
	IntegratorNumerical(100,0);
}



IntegratorNumerical::IntegratorNumerical(int trials, float fd){
	this->trials=trials;
	this->fd = fd;
	
}


float IntegratorNumerical::angle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}



float IntegratorNumerical::integrate(Molecule *molecule, Tessellation *tessellation){
	exit(-1);
}



float IntegratorNumerical::integrate(Molecule *molecule, int index, void (*feedback)(float)){
	Vector coordinate(3);
	Vector pp(3),pn(3);
	float areap,arean;
	vector<vector<ForcesDT*> > forces;	
	
	
	//numerical derivatives
	if(fd!=0){
		atoms = molecule->fetchCoordinates();
		forces = molecule->fetchForcePointers();
		
		
		
		for(unsigned int i=0; i<atoms.size(); ++i){
			coordinate = molecule->getInternallyStoredAtomCoordinates(i);
			for(int j=0; j<3; ++j){
				if(feedback!=NULL) feedback((float)(i*3+j)/((float)(atoms.size()*3)-1));

				pp = Vector(3).zeros();
				pn = Vector(3).zeros();
				pp(j) = fd;
				pn(j) = -fd;
				
				
				molecule->perturbInternallyStoredAtomCoordinates(i, pp);
				areap = integrateMolecule(molecule,-1);
				molecule->setInternallyStoredAtomCoordinates(i, coordinate);
				
				molecule->perturbInternallyStoredAtomCoordinates(i, pn);
				arean = integrateMolecule(molecule,-1);
				molecule->setInternallyStoredAtomCoordinates(i, coordinate);
				
				*(forces[i][j]) = (areap - arean)/(2*fd);
			}
			
		}		
	}
	
	return integrateMolecule(molecule, index);
}


/*
float IntegratorNumerical::partitionSearch(int p0, int p1, float phi0, float theta0, float phi1, float theta1, vector<float> &lambdas, vector<Vector> &mus, vector<CircularInterfaceForm> &forms){
	Vector v0 = Vector(2);
	Vector v1 = Vector(2);
	v0(0)=phi0;
	v0(1)=theta0;
	v1(0)=phi1;
	v1(1)=theta1;
	Vector l=v0;
	Vector r=v1;
	float length = norm(l-r,2);
	float h = 1;
	Vector c=r;
	int p;
	//printf("[%f %f - %f %f (%f) %d %d]\n",phi0,theta0,phi1,theta1,length,p0,p1);
	while(h>partitionThreshold){
		c = r+(l-r)*0.5;
		p = occludes(c(0),c(1), lambdas, mus, forms);
		if(p0==p) l=c;
		else r=c;
		
		h = norm(l-r,2)/length;
		//printf("%f %f (%f) [%d]\n",c(0),c(1),h,p);
		
	}
	//printf("[%f]\n",norm(v0-c,2)/length);
	
	return norm(v0-c,2)/length;
}
*/


Vector IntegratorNumerical::sphericalVector(float phi, float theta){
	Vector v(3);
	// Find the position of the current surface segment.
	v(0) = 1;
	v(1) = theta;
	v(2) = phi;
	v = spherical2cartesian(v);
	
	return v;
}

	
bool IntegratorNumerical::occludes(Vector v, vector<float> &lambdas, vector<Vector> &mus, vector<CircularInterfaceForm> &forms, float &conflict){
	bool occluded;
	CircularInterfaceForm form;
	float l;
	Vector n;
	float conflict_sasa, conflict_sesa;
	
	conflict_sasa=99999999999;
	
	occluded = false;
	
	for(unsigned int j=0; j<mus.size() && !occluded; ++j){
		n = mus[j];
		l = lambdas[j];
		form = forms[j];
		
		float a=angle(v,n);
	
		if(form==CONVEX){
			if(a <= l){
				occluded=true; 
				conflict_sesa=l-a;
			}
			else{
				if(a-l<conflict_sasa) conflict_sasa=a-l;
			}
		}
		else{
			if(a >= l){
				occluded=true; 
				conflict_sesa=a-l;
			}
			else{
				if(l-a<conflict_sasa) conflict_sasa=l-a;
			}
				
		}
	}
	
	if(occluded) conflict=conflict_sesa;
	else conflict=conflict_sasa;
	
	

	
	
	return occluded;
	
}


float IntegratorNumerical::integrateMolecule(Molecule *molecule, int index){
	Vector origin;
	float r_k;
	float lenv;
	Vector normal;
	Vector mu;
	float lambda;
	CircularInterfaceForm form;
	
	vector<float> lambdas;
	vector<Vector> mus;
	vector<CircularInterfaceForm> forms;
	
	float area;
	Vector v(3);
	float sasaCount;
	float sesaCount;
	bool occluded;
	int detail;
	boost::mt19937 randomizer;
	static boost::uniform_real<> dist01(0, 1);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > gen01(randomizer, dist01);
	float phi_step;
	float dtheta;
	float dtheta2;
	float dots_in_ring;
	float theta_step;
	float temp_dtheta;
	float g,g_normalised;
	Vector n;
	float a;
	float d_i;
	float r_i,r_l,radius;
	bool internalAtom;
	Vector v0,v1,v2,v3;
	float x;
	float conflict,conflict1;
	
	Vector lineNew, lineOld;
	
	
	area = 0;
	srand(120);
	randomizer.seed(0);
	detail = trials;
	phi_step = M_PI / detail;
	
	radii = molecule->fetchRadii();
	atoms = molecule->fetchCoordinates();
	areas = molecule->fetchAreaPointers();
	
	
	if(index>=0)
		file = fopen ("rays.csv","w");
	
	
	for(unsigned int i=0; i<atoms.size(); ++i){
		//i=12;	{
		origin = atoms[i];
		radius = radii[i];
		mus.clear();
		lambdas.clear();
		forms.clear();
		internalAtom = false;
		
		if(index == (signed int) i) fprintf(file, "atom radius %f\n", radius);
		
		
		for(unsigned int j=0; j<atoms.size(); ++j){
			if(i!=j){
				
				v=atoms[j] - origin;
				r_k = radii[j];
				r_i = radius;
				
				if(index==(signed int) i) fprintf(file, "intersector %d radius %f vector %f %f %f\n", j, r_k, v(0), v(1), v(2));
				
				
				
				lenv = norm(v,2);
				

				//reject, if no intersection
				if(lenv < r_i + r_k && lenv+r_k > r_i && lenv+r_i > r_k){
					normal = v / lenv;
					d_i = norm(v,2);
					mu = normal;
					
					
					r_l = radius;
					r_i = r_k;
		
		
					g = (d_i * d_i + r_l * r_l - r_i * r_i ) / (2 * d_i);
		
					if(g<0){
						g = abs(g);
						normal=-normal;
						form=CONCAVE;
					}
					else form = CONVEX;

					g_normalised = g/radius;
					lambda = acos(g_normalised);
					
					mus.push_back(normal);
					lambdas.push_back(lambda);
					forms.push_back(form);
				}
				else if(r_i + lenv < r_k){
					internalAtom = true;
					break;
				}
	
					
					
					
			}
				
				
				
			
		}
		
		
		
		
		
		
		
		if(internalAtom){
			*(areas[i]) = 0;
			continue;
		}
		
		
		sasaCount=0;
		sesaCount=0;
		
		// Set up a random reference rotation.
		//  static boost::mt19937 rand_gen(time(0));

		// Integrate over the surface as a set of rings.
		dtheta = gen01();
		if (dtheta < 0.5) {
			dtheta2 = 0.5 + dtheta;
		}
		else{
			dtheta2 = dtheta - 0.5;
		}
		
		for (float phi = phi_step * gen01(); phi < M_PI; phi += phi_step) {
			// Use an integer number of dots for uniform surface coverage.
			dots_in_ring = int(2*detail*sin(phi));
			theta_step = 2*M_PI / dots_in_ring;
			// To avoid a seam on the sphere, we alternate between dtheta and dtheta2
			temp_dtheta = dtheta;
			dtheta = dtheta2;
			dtheta2 = temp_dtheta;
			
			for (float theta = theta_step * dtheta;	theta < 2 * M_PI; theta += theta_step) {
				
				v = sphericalVector(phi,theta);

				
				if(occludes(v,lambdas,mus,forms,conflict)){
					v0 = sphericalVector(phi-phi_step,theta);
					v1 = sphericalVector(phi+phi_step,theta);
					v2 = sphericalVector(phi,theta-theta_step);
					v3 = sphericalVector(phi,theta+theta_step);
					
					x=0;
					
					if(!occludes(v0,lambdas,mus,forms,conflict1))
						x += conflict/(conflict+conflict1);
					else x+=1;
					if(!occludes(v1,lambdas,mus,forms,conflict1))
						x += conflict/(conflict+conflict1);
					else x+=1;
					if(!occludes(v2,lambdas,mus,forms,conflict1))
						x += conflict/(conflict+conflict1);
					else x+=1;
					if(!occludes(v3,lambdas,mus,forms,conflict1))
						x += conflict/(conflict+conflict1);
					else x+=1;

					sasaCount+=1.0-(x/4.0);
					sesaCount+=x/4.0;
					
					occluded=true;
					
					
					//sesaCount+=1.0;
						
						
						
				}
				else{
					
					v0 = sphericalVector(phi-phi_step,theta);
					v1 = sphericalVector(phi+phi_step,theta);
					v2 = sphericalVector(phi,theta-theta_step);
					v3 = sphericalVector(phi,theta+theta_step);
					
					x=0;
					
					if(occludes(v0,lambdas,mus,forms,conflict1))
						x += conflict/(conflict+conflict1);
					else x+=1;
					if(occludes(v1,lambdas,mus,forms,conflict1))
						x += conflict/(conflict+conflict1);
					else x+=1;
					if(occludes(v2,lambdas,mus,forms,conflict1))
						x += conflict/(conflict+conflict1);
					else x+=1;
					if(occludes(v3,lambdas,mus,forms,conflict1))
						x += conflict/(conflict+conflict1);
					else x+=1;

					sesaCount+=1.0-(x/4.0);
					sasaCount+=x/4.0;
					
					
					occluded=false;
					
					
					
					//sasaCount+=1.0;
					
					
					
				}
				
			
				if(index==(signed int)i) fprintf(file, "ray %f %f %f occluded %d\n", v(0),v(1),v(2),occluded);
				
			}
		}
		
//		printf("SASA: %d, SESA: %d\n",sasaCount,sesaCount);
		
		a = 4*M_PI*(float)sasaCount / ((float)(sasaCount+sesaCount));
		
		a = radius*radius * a;
				
		*(areas[i]) = a;
		area += a;
	
	}
	
	if(index>=0)
		fclose(file);
	
	return area;
	
	
}




Vector IntegratorNumerical::spherical2cartesian(Vector s) {
	Vector v(3);
	v(0) = s(0) * cos(s(1)) * sin(s(2));
	v(1) = s(0) * sin(s(1)) * sin(s(2));
	v(2) = s(0) * cos(s(2));
	return v;
}



