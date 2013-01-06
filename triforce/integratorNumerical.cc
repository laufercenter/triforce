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
	IntegratorNumerical(100);
}



IntegratorNumerical::IntegratorNumerical(int trials){
	this->trials=trials;
	
}


double IntegratorNumerical::angle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}



double IntegratorNumerical::integrate(Molecule *molecule){
	Vector origin;
	double r_k;
	double lenv;
	Vector normal;
	Vector mu;
	double lambda;
	CircularInterfaceForm form;
	
	vector<double> lambdas;
	vector<Vector> mus;
	vector<CircularInterfaceForm> forms;
	
	double area;
	Vector v(3);
	int sasaCount;
	int sesaCount;
	bool occluded;
	int detail;
	boost::mt19937 randomizer;
	static boost::uniform_real<> dist01(0, 1);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > gen01(randomizer, dist01);
	float phi_step;
	int total_segs;
	float dtheta;
	float dtheta2;
	float dots_in_ring;
	float theta_step;
	float temp_dtheta;
	int i,j;
	double g,g_normalised;
	double l;
	Vector n;
	double a;
	double d_i;
	double r_i,r_l,radius;
	
	
	
	area = 0;
	srand(120);
	randomizer.seed(0);
	detail = trials;
	phi_step = M_PI / detail;
	
	radii = molecule->fetchRadii();
	atoms = molecule->fetchCoordinates();
	areas = molecule->fetchAreaPointers();
	
	file = fopen ("gbonnet1.csv","a");
	
	
	//for(i=0; i<atoms.size(); ++i){
		i=0;
		{
		origin = atoms[i];
		radius = radii[i];
		mus.clear();
		lambdas.clear();
		forms.clear();
		for(j=0; j<atoms.size(); ++j){
			if(i!=j){
				
				v=atoms[j] - origin;
				r_k = radii[j];
				r_i = radius;
				
				
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
			}
				
				
				
			
		}
		
		
		
		sasaCount=0;
		sesaCount=0;
		total_segs = 0;
		
		// Set up a random reference rotation.
		//  static boost::mt19937 rand_gen(time(0));

		// Integrate over the surface as a set of rings.
		dtheta = gen01();
		dtheta2;
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
		
				// Find the position of the current surface segment.
				v(0) = 1;
				v(1) = theta;
				v(2) = phi;
				v = spherical2cartesian(v);
				
				
				occluded = false;
				
				for(j=0; j<mus.size() && !occluded; ++j){
					n = mus[j];
					l = lambdas[j];
					form = forms[j];
					
				
					if(form==CONVEX){
						if(angle(v,n) <= l) occluded=true; 
					}
					else{
						//n = -n;
						if(angle(v,n) >= l) occluded=true; 
					}
				}
			
					
			
				fprintf(file, "ray %f %f %f occluded %d\n", v(0),v(1),v(2),occluded);
				
				if(!occluded) sasaCount++;
				else sesaCount++;
			}
		}
		
		printf("SASA: %d, SESA: %d\n",sasaCount,sesaCount);
		
		a = 4*M_PI*(double)sasaCount / ((double)(sasaCount+sesaCount));
		
		a = radius*radius * a;
				
		*(areas[i]) = a;
		area += a;
	
	}
	
	
	fclose(file);
	
	return area;
	
	
}


double IntegratorNumerical::integrate(Molecule *molecule, Tessellation *tessellation){
	this->molecule = molecule;
	this->tessellation = tessellation;
	
	SASAsForMolecule sasas;
	double area,a;
	double radius;
	vector<double> radii;
	
	radii = molecule->fetchRadii();
	areas = molecule->fetchAreaPointers();
	
	
	sasas = tessellation->sasas();
	
	
	file = fopen ("gbonnet0.csv","a");

	
	area = 0;
	//iterate over all atoms
	for(int i=0;i<sasas.size();++i){
		radius = sasas[i].radius;
		a = radius*radius * integrateAtomicSASA(sasas[i]);
		*(areas[i]) = a;
		area += a;
	}
	
	return area;
	
	fclose(file);
	
}



Vector IntegratorNumerical::spherical2cartesian(Vector s) {
	Vector v(3);
	v(0) = s(0) * cos(s(1)) * sin(s(2));
	v(1) = s(0) * sin(s(1)) * sin(s(2));
	v(2) = s(0) * cos(s(2));
	return v;
}



double IntegratorNumerical::integrateAtomicSASA(SASAsForAtom sasasForAtom){
	double radius;
	double area;
	Vector v(3);
	int form;
	int sasaCount;
	int sesaCount;
	bool occluded;
	Vector n(3);
	double l;
	int detail;
	boost::mt19937 randomizer;
	static boost::uniform_real<> dist01(0, 1);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > gen01(randomizer, dist01);
	float phi_step;
	int total_segs;
	float dtheta;
	float dtheta2;
	float dots_in_ring;
	float theta_step;
	float temp_dtheta;
	
	
	radius = sasasForAtom.radius;
	area = 0;
	sasaCount=0;
	sesaCount=0;
	srand(120);
	
	randomizer.seed(0);
		
	detail = trials;
		

	phi_step = M_PI / detail;
	total_segs = 0;
	
	// Set up a random reference rotation.
	//  static boost::mt19937 rand_gen(time(0));

	// Integrate over the surface as a set of rings.
	dtheta = gen01();
	dtheta2;
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
	
			// Find the position of the current surface segment.
			v(0) = 1;
			v(1) = theta;
			v(2) = phi;
			v = spherical2cartesian(v);
			
			
			occluded = false;
			for(int i=0; !occluded && i<sasasForAtom.sasas.size();++i){
				//int i =1;
				for(int j=0; !occluded && j<sasasForAtom.sasas[i].sasa.size(); ++j){
					n = sasasForAtom.sasas[i].sasa[j].normalForCircularInterface;
					l = sasasForAtom.sasas[i].sasa[j].lambda.rotation;
					form = sasasForAtom.sasas[i].sasa[j].form;
					
					
					//printf("C[%d,%d] (%f,%f,%f) {%f}\n",sasasForAtom.sasas[i].sasa[j].id0,sasasForAtom.sasas[i].sasa[j].id1,sasasForAtom.sasas[i].sasa[j].normalForCircularInterface(0),sasasForAtom.sasas[i].sasa[j].normalForCircularInterface(1),sasasForAtom.sasas[i].sasa[j].normalForCircularInterface(2),sasasForAtom.sasas[i].sasa[j].lambda);
					
					if(form != SPLITTER){
						
						if(form==CONVEX){
							if(angle(v,n) <= l) occluded=true; 
							
							//printf("ANGLE: %f, LAMBDA: %f\n",angle(v,n),l);
						}
						else{
							//n = -n;
							if(angle(v,n) >= l) occluded=true; 
						}
					}
					n(0)=1;
					n(1)=0;
					n(2)=0;
					//if(angle(v,n) <= M_PI/2) occluded=true; 
					
					

							
				}
				
			}

			fprintf(file, "ray %f %f %f occluded %d\n", v(0),v(1),v(2),occluded);
			
			if(!occluded) sasaCount++;
			else sesaCount++;
		}
	}
	
	printf("SASA: %d, SESA: %d\n",sasaCount,sesaCount);
	
	area = 4*M_PI*(double)sasaCount / ((double)(sasaCount+sesaCount));
	
	return area;
	
}
