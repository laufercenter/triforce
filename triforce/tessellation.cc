#include "tessellation.h"



Tessellation::Tessellation(Molecule &m){
	molecule = m;
}



void Tessellation::build(){
	//molecule.update();
	atoms = molecule.coordinates();
	radii = molecule.fetchRadii();
	emptyCircularRegions();
	emptyIntersections();
	
	vector<CircularRegion>* circles;
	vector<list<IntersectionPoint*>*>* intersections;
	
	
	//iterate over all atoms and build the tesselation for each of them
	//for(int i=0; i<atoms.size(); ++i){
		int i=0;
		buildGaussBonnetPath(atoms[i], radii->at(i), atoms, *radii, &circles, &intersections);
		circleSet.push_back(circles);
		intersectionSet.push_back(intersections);
		
		
		for(int i=0; i<intersectionSet[0]->size(); ++i){
			printf("SASA %d\n",i);
			outputGaussBonnetPath(*(intersectionSet[0]->at(i)));
		}


		
		//outputCircularRegions(*circles[circles.size()-1]);

	//}
}

void Tessellation::emptyIntersections(){
	vector<vector<list<IntersectionPoint*>*>*>::iterator it0;
	vector<list<IntersectionPoint*>*>::iterator it1;
	
	for(it0 = intersectionSet.begin(); it0 != intersectionSet.end(); ++it0){
		for(it1 = (*it0)->begin(); it1 != (*it0)->end(); ++it1){
			delete (*it1);
		}
		delete (*it0);
	}
	intersectionSet.clear();
}

void Tessellation::emptyCircularRegions(){
	vector<vector<CircularRegion> *>::iterator it0;
	for(it0 = circleSet.begin(); it0!= circleSet.end(); ++it0){
		delete (*it0);
	}
	circleSet.clear();
	
}

vector<vector<list<IntersectionPoint*>*>*>* Tessellation::intersectionPoints(){
	return &intersectionSet;
}
vector<vector<CircularRegion>* >* Tessellation::circularRegions(){
	return &circleSet;
}



double Tessellation::vsign(double v)
{
	if(v>=0) return 1.0;
	else return -1.0;
}

double Tessellation::cot(double a){
	return 1.0/tan(a);
}

double Tessellation::csc(double a){
	return 1.0/sin(a);
}




double Tessellation::getAngleBetweenNormals(Vector &a, Vector &b){
	return acos(dot(a,b));
}


double Tessellation::getAngle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}


bool Tessellation::isZero(double v){
	if(fabs(v) <= THRESHOLD_NUMERICAL) return true;
	else return false;
}



void Tessellation::determineProjection(Vector &origin, double radius, CircularRegion &circle){
	double d_k;
	double r_i, r_k;
	double g;
	Vector mu;
	double a;
	d_k = norm(circle.vector,2);
	mu = circle.normal;
	r_i = radius;
	r_k = circle.sphereRadius;
	g = (d_k * d_k + r_i * r_i - r_k * r_k ) / (2 * d_k);
	a = sqrt(r_i * r_i - g * g);
	
	if(g<0) circle.form=CONCAVE;
	else circle.form=CONVEX;

	circle.openingAngle = acos(g / r_i);
	circle.g = g;
	circle.a = a;

	//NSWarnLog(@"DET: [%f %f %f] %f %f %f %f",circle.normal.vector[0],circle.normal.vector[1],circle.normal.vector[2], r_i, r_k, g, circle.openingAngle);


	//NSWarnLog(@"%f %f %f %f %f",d_k, r_i, r_k, g, circle.openingAngle);
}


IntersectionPair Tessellation::determineIntersectionPoints(double radius, CircularRegion &K, CircularRegion &J){
	IntersectionPair res;
	double g_k, g_j;
	double r_i, r_k, r_j;
	double a_k, a_j;
	Vector mu_k, mu_j;
	double phi_kj, phi_jk;
	double tau_kj,tau_jk;
	Vector tmp1, tmp2;
	Vector eta_kj, eta_jk;
	Vector omega_kj, omega_jk;
	Vector p_kj, p_jk;
	double gamma_kj,gamma_jk;





	g_k = K.g;
	g_j = J.g;

	r_i = radius;
	r_k = K.sphereRadius;
	r_j = J.sphereRadius;


	//NSWarnLog(@"INT: %f [%f %f %f] %f %f, [%f %f %f] %f %f",r_i,K.vector.vector[0],K.vector.vector[1],K.vector.vector[2],K.g,r_k,J.vector.vector[0],J.vector.vector[1],J.vector.vector[2],J.g,r_j);


	//a_k = sqrt(r_i * r_i - g_k * g_k);
	//a_j = sqrt(r_i * r_i - g_j * g_j);

	mu_k = K.normal;
	mu_j = J.normal;

	phi_kj = acos(dot(mu_k,mu_j));
	phi_jk = acos(dot(mu_j,mu_k));

	tau_kj = (g_k - g_j * cos(phi_kj)) / (pow(sin(phi_kj),2));
	tau_jk = (g_j - g_k * cos(phi_jk)) / (pow(sin(phi_jk),2));


	eta_jk = eta_kj = (mu_k * tau_kj) + (mu_j * tau_jk);

	omega_kj = cross(mu_k, mu_j) / sin(phi_kj);
	omega_jk = cross(mu_j, mu_k) / sin(phi_jk);

	gamma_kj = sqrt(r_i * r_i - g_k*tau_kj - g_j*tau_jk);
	gamma_jk = sqrt(r_i * r_i - g_j*tau_jk - g_k*tau_kj);

	p_kj = (omega_kj * gamma_kj) + eta_kj;
	p_jk = (omega_jk * gamma_jk) + eta_jk;
	
	
	res.k_j = p_kj;
	res.j_k = p_jk;


	return res;
	

	
}



void Tessellation::makeCircularRegions(Vector &origin, double radius, vector<vec> &atoms, vector<double> &radii, vector<CircularRegion> &circles){
	int i;
	Vector res,normal,v;
	CircularRegion circle;
	double r_i = radius;
	double r_k;
	double lenv;

	for(i=0;i<atoms.size();i++){
		v=atoms[i] - origin;
		r_k = radii[i];
		
		
		lenv = norm(v,2);
		
		//printf("RADIUS: %f LEN: %f\n",r_k,lenv);
		

		//reject, if no intersection
		if(lenv < r_i + r_k && lenv+r_k > r_i){
			normal = v / lenv;
			circle.id = circles.size();
			circle.vector = v;
			circle.normal = normal;
			circle.sphereRadius = r_k;
			circle.intersect = false;
			circles.push_back(circle);
			
			//printf("CIRCLE[%d]: (%f, %f, %f) %f\n",circle.id, circle.vector(0), circle.vector(1), circle.vector(2), circle.sphereRadius);
		}
	}
}



int Tessellation::filterCircularRegions(double radius, vector<CircularRegion> &circles){
	vector<CircularRegion>::iterator it;
	double angle;
	int i;
	bool erased;
	it = circles.begin();
	while(it != circles.end()){
		erased=false;
		for(i=0;i<circles.size();i++){
			if(it->id != circles[i].id){			
				angle = getAngleBetweenNormals(it->normal, circles[i].normal);
//				NSWarnLog(@"ANGLE: %f %f %f",angle,it->openingAngle, circles[i].openingAngle);
				
				if(it->form==CONVEX){
					if(circles[i].form == CONVEX){
						//the circle is inside of convex circle i
						if(it->openingAngle + angle < circles[i].openingAngle){
							it = circles.erase(it);
							erased=true;
							break;
						}
					}
					else{
						//the circle is outside of concave circle i
						if(angle > it->openingAngle + circles[i].openingAngle){
							it = circles.erase(it);
							erased=true;
							break;
						}
					}
					
				}
				else{
					if(circles[i].form == CONVEX){
						//the sasa part of a concave circle is being covered by a convex circle
						if(it->openingAngle + angle < circles[i].openingAngle){
							//in this special case, the atom is completely buried and the sasa is 0
							return -1;
						}
					}
					else{
						//the concave circle is completely surrounding another concave circle
						if(it->openingAngle > circles[i].openingAngle + angle){
							it = circles.erase(it);
							erased=true;
							break;
						}
					}

				}		
						
					
			}
		}
		if(!erased) ++it;
	}
	
	return 0;

}

void Tessellation::filterEmptyCircularRegions(vector<CircularRegion> &circles){
	vector<CircularRegion>::iterator it;
	bool erased;
	it = circles.begin();
	while(it != circles.end()){
		erased=false;
		if(it->forwardIntersections.size()==0){
			it = circles.erase(it);
			erased=true;
		}
		if(!erased) ++it;
	}
	
}



void Tessellation::reindexCircularRegions(vector<CircularRegion> &circles){
	list<IntersectionPoint>::iterator it;
	int i,j;
	int currentID;
	for(i=0;i<circles.size();i++){
		currentID = circles[i].id;
		circles[i].id = i;
		for(it = circles[i].forwardIntersections.begin(); it != circles[i].forwardIntersections.end(); ++it){
			it->from = i;
		}

		for(j=0;j<circles.size();j++){
			if(j!=i){
				it = circles[j].forwardIntersections.begin();
				while(it != circles[j].forwardIntersections.end()){
					if(it->with == currentID && !it->visited){
						it->with = i;
						it->visited=true;
					}
					
					++it;
				}

				
			}
		}
		
	}
}


void Tessellation::outputCircularRegions(vector<CircularRegion> &circles){
	list<IntersectionPoint>::iterator it;
	int i,u;
	fprintf(stdout, "......... outputting circular regions ...........\n");
	for(i=0;i<circles.size();i++){
		fprintf(stdout,"- circular region %d with id: %d:\n",i,circles[i].id);
		u=0;
		for(it = circles[i].forwardIntersections.begin(); it != circles[i].forwardIntersections.end(); ++it){
			fprintf(stdout,"intersection %d from: %d with: %d (%f %f %f)\n",u,it->from,it->with,it->vector(0),it->vector(1),it->vector(2));
			++u;
		}
	}
	fprintf(stdout,"......... end of data  ...........\n");

}

/*
void outputSingleCircularRegion(CircularRegion *circle){
	
	NSWarnLog(@"...............................................");
	NSWarnLog(@"circular region [%d / %d]:",circle->id, circle->index);
	NSWarnLog(@"vector: (%f %f %f) normal: (%f %f %f)",circle->vector.vector[0],circle->vector.vector[1],circle->vector.vector[2],circle->normal.vector[0],circle->normal.vector[1],circle->normal.vector[2]);
	NSWarnLog(@"length: %f(%f), g: %f, a: %f, openingAngle: %f, sphereRadius: %f",circle->vector.length,circle->normal.length,circle->g, circle->a, circle->openingAngle, circle->sphereRadius); 
	
}
*/


void Tessellation::filterIntersectionPoints(vector<CircularRegion> &circles, int except){
	list<IntersectionPoint>::iterator it;
	int i,u;
	double angle;
	bool erase;
	for(i=0; i<circles.size(); i++){
		it = circles[i].forwardIntersections.begin();
		while(it != circles[i].forwardIntersections.end()){
			erase=false;
			for(u=0; u<circles.size(); u++){
				if(it->with != circles[u].id && circles[u].id!=circles[i].id && circles[u].id!=except){
					angle = getAngle(it->vector, circles[u].vector);
					if(angle < circles[u].openingAngle) erase = true;
					
					if(erase){
						it = circles[i].forwardIntersections.erase(it);
						break;
					}
				}
			}
			if(!erase) ++it;
		}

	}

}



void Tessellation::clearFlags(vector<CircularRegion> &circles){
	list<IntersectionPoint>::iterator it;
	int i;	
	for(i=0;i<circles.size();i++){
		it = circles[i].forwardIntersections.begin();
		for(it = circles[i].forwardIntersections.begin(); it != circles[i].forwardIntersections.end(); ++it){
			it->visited=false;
		}
	}
}


vector<CircularRegion>* Tessellation::deepCopy(vector<CircularRegion> &circles){
	int i,j;
	list<IntersectionPoint>::iterator it;
	vector<CircularRegion>* newCircles;
	newCircles = new vector<CircularRegion>();
	
	CircularRegion c;
	for(i=0;i<circles.size();i++){
		c.id=circles[i].id;
		c.vector=circles[i].vector;
		c.normal=circles[i].normal;
		c.openingAngle=circles[i].openingAngle;
		c.g=circles[i].g;
		c.a=circles[i].a;
		c.sphereRadius=circles[i].sphereRadius;
		c.form = circles[i].form;
		c.forwardIntersections.clear();
		for(it = circles[i].forwardIntersections.begin(); it != circles[i].forwardIntersections.end(); ++it){
			c.forwardIntersections.push_back(*it);	
		}
		newCircles->push_back(c);
	}
	return newCircles;
}


void Tessellation::outputGaussBonnetPath(list<IntersectionPoint*> &points){
	list<IntersectionPoint*>::iterator it;
	int i;
	
	for(it=points.begin(), i=0; it!=points.end(); ++it, ++i)
		printf("GBPATH[%d] %d - %d - %d\n", i, (*it)->from, (*it)->id, (*it)->with);
}



void Tessellation::prepareCircularRegions(vector<CircularRegion> &circles, vector<CircularRegion> **newCircles){
	*newCircles = deepCopy(circles);
	
	filterEmptyCircularRegions(**newCircles);
	clearFlags(**newCircles);
	reindexCircularRegions(**newCircles);
	clearFlags(**newCircles);
}




void Tessellation::insertFakeIntersectionPoints(vector<CircularRegion> &circles){
	Vector v,v2,o,o2;
	int k;
	IntersectionPoint p;
	
	//we have to find circles that did not intersect with other circles and we have to insert fake intersectionpoints for them
	for(k=0;k<circles.size();k++){
		if(circles[k].circularIntersections.size()==0){
			printf("INSERTING FAKE INTERSECTIONPOINTS\n");
			v = randu<vec>(3);
			v2 = circles[k].normal * circles[k].g;
			o = cross(v2,v);
			o = circles[k].a*o / norm(o,2);
			p.vector = v2+o;
			p.with=k;
			p.from=k;
			circles[k].forwardIntersections.push_back(p);
			
			o2 = cross(v2,o);
			o2 = circles[k].a*o2 / norm(o2,2);
			p.vector=v2+o2;
			circles[k].forwardIntersections.push_back(p);
			
			p.vector = v2 + (o * -1);
			circles[k].forwardIntersections.push_back(p);
			
			p.vector = v2 + (o2 * -1);
			circles[k].forwardIntersections.push_back(p);
			
		}
	}
	
}

int Tessellation::sgn(double d){
	if(d>=0) return 1;
	else return -1;
}

void Tessellation::determineCircularIntersections(vector<CircularRegion> &circles){
	int i,k,j;
	double angle;
	CircularIntersection c;
	
	for(k=0;k<circles.size();k++)
		for(j=k+1;j<circles.size();j++){
			if(k!=j){
				//determine if there will be intersections
				angle = getAngleBetweenNormals(circles[k].normal,circles[j].normal);
				if(angle < circles[k].openingAngle + circles[j].openingAngle)
					if(angle + circles[k].openingAngle > circles[j].openingAngle && angle + circles[j].openingAngle > circles[k].openingAngle){
						c.d=angle;
						c.visited=false;
						c.blocked=false;
						c.sasa=false;
						circles[j].circularIntersections[k]=c;
						circles[k].circularIntersections[j]=c;
						
						//printf("circle intersection %d-%d\n",k,j);
					}
			}
		}
	
}



/**
 * n is a normal vector to the plane connecting integration origin, origin and center of circular region
 * o is a vector pointing in direction of integration origin but intersects with the interface of circular region
 * a is one of the PHI vectors from the center of circular region to its interface
 * 
 * function will return the angle between a and the plane between 0 and PI (instead of 0 and PI/2)
 * 
 */
double Tessellation::complLongAngle(Vector &vi, Vector &vj, Vector &vk){
	double v;
	int s;
	Vector nij(3), nik(3), vij(3), vik(3);
	Vector ni(3);
	
	
	vij = vj - vi;
	nij = cross(vi, vij);
	nij = nij/norm(nij,2);
	vik = vk - vi;
	nik = cross(vi, vik);
	nik = nik/norm(nik,2);
	
	ni = cross(vi,nij);
	
	
	
	v = acos(dot(nij,nik));
	s = sgn(dot(ni,nik));
	
	if(s>0) v = -v;
	   
	return v;
		
}

/*
dot <- function(a,b){
a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
}
 * 
 * */



double Tessellation::intersectionDistance(CircularRegion I, CircularRegion J, CircularRegion &K, double dij, double djk){
	CircularRegion T;
	double gij, gik, dik;
	double eta, mu;
	double fmin;
	
	
	
	T=I;
	I=J;
	J=T;
	dik=djk;
	
	gij = (dij*dij - J.openingAngle*J.openingAngle + I.openingAngle*I.openingAngle) / (2*dij*I.openingAngle);
	gik = (dik*dik - K.openingAngle*K.openingAngle + I.openingAngle*I.openingAngle) / (2*dik*I.openingAngle);
	
	eta = complLongAngle(I.normal,J.normal,K.normal);
	
	fmin = acos(gij)+acos(gik);
	
	mu = eta-fmin;
	if(mu < 0) mu = 2*M_PI+mu;
	
	return mu;
		
		
	
}
/*
int Tessellation::selectNextIntersectingCluster(CircularRegion &I, CircularRegion &J, vector<CircularRegion> &circles){
	map<int,CircularIntersection>::iterator it;
	map<int,CircularIntersection>::iterator it_l;
	double minDist,d;
	double minClusterDist;
	int k;
	vector<int> clusters;
	map<int,bool> used;
	int n;
	
	
	for(it=J.circularIntersections.begin(); it != J.circularIntersections.end(); ++it){
		k = it->first;
		if(used.find(k)==J.circularIntersections.end()){
			clusters.push_back(k);
			for(it_l = circles[k].circularIntersections.begin(); it_l = circles[k].circularIntersections.end(); ++it_l){
				used[it_l->first] = true;
			}
		}
	}
	
	k=-1;
	
	minClusterDist = numeric_limits<double>::max();
	for(m=0; m<clusters.size(); ++m){
		minDist = numeric_limits<double>::max();
		for(it=circles[m].circularIntersections.begin(); it != circles[m].circularIntersections.end(); ++it){
			if(it->first != I.id) d = intersectionDistance(I,J,circles[it->first],I.circularIntersections[J.id].d, it->second.d);
			else d = 2*M_PI;
			if(d<minDist){
				minDist = d;
			}
			
		}
		
		if(minDist<minClusterDist){
			minClusterDist = minDist;
			k = m;
		}
		
	}
	
	
	printf("I: %d J: %d K: %d bestk: %d, d: %f, mind: %f\n",I.id,J.id,it->first,k,d,minDist);
	
	//if(k>=0 && (J.circularIntersections[k].blocked || J.circularIntersections[k].visited)) k = -1;
	
	return k;

}
*/

int Tessellation::selectNextIntersectingCluster(CircularRegion &I, CircularRegion &J, vector<CircularRegion> &circles){
	map<int,CircularIntersection>::iterator it;
	map<int,CircularIntersection>::iterator it_l;
	double minDist,d;
	double minClusterDist;
	int k;
	vector<int> clusters;
	map<int,bool> used;
	int n;
	
	k=-1;
	
	minDist = numeric_limits<double>::max();
	for(it=J.circularIntersections.begin(); it != J.circularIntersections.end(); ++it){
		if(it->first != I.id) d = intersectionDistance(I,J,circles[it->first],I.circularIntersections[J.id].d, it->second.d);
		else d = 2*M_PI;
		if(d<minDist){
			minDist = d;
			k = it->first;
		}
		printf("I: %d J: %d K: %d bestk: %d, d: %f, mind: %f\n",I.id,J.id,it->first,k,d,minDist);
		
	}
	
	
	
	//if(k>=0 && (J.circularIntersections[k].blocked || J.circularIntersections[k].visited)) k = -1;
	
	return k;

}

OcclusionType Tessellation::occludesForwardIntersectionPoint(CircularRegion &I, CircularRegion &J, CircularRegion &K, double dij, double dik){
	Vector vij(3);
	Vector vik(3);
	Vector ni(3);
	Vector nij(3);
	
	double gij, gik, eta;
	double fmin, fmax, fmax2;
	
	bool top, bottom, both, inside, outside, triforce;
	
	top=false;
	bottom=false;
	both=false;
	inside=false;
	outside=false;
	triforce=false;
	
	
	
	// g determines how far j (or k) invades i, where as 1 is no invasion, 0 is half invasion and -1 complete invasion
	gij = (dij*dij - J.openingAngle*J.openingAngle + I.openingAngle*I.openingAngle) / (2*dij*I.openingAngle);
	gik = (dik*dik - K.openingAngle*K.openingAngle + I.openingAngle*I.openingAngle) / (2*dik*I.openingAngle);
	
	eta = complLongAngle(I.normal,J.normal,K.normal);
	
	fmin = -(acos(gij)+acos(gik));
	fmax = -(acos(gij)-acos(gik));
	fmax2 = fmax+2*M_PI;
	
	//printf("dij: %f, dik: %f, iop: %f, jop: %f, kop: %f\n",dij,dik,I.openingAngle,J.openingAngle,K.openingAngle);
	//printf("gij: %f, gik: %f, fmin: %f, fmax: %f, fmax2: %f, eta: %f\n",gij,gik,fmin,fmax,fmax2,eta);
	
	if(K.form==CONVEX){
		//bottom occluded
		if((eta < fmin && eta > fmax) || eta > fmax2){
			bottom=true;
			//printf("bottom\n");
		}
		//top occluded
		if((eta > -fmin && eta < -fmax) || eta < -fmax2){
			top=true;
			//printf("top\n");
		}
		//total occlusion
		if(bottom && top){
			both=true;
			//printf("both\n");
		}
		//no occlusion, inside intersection
		if(!bottom && ! top && eta > fmin && eta < -fmin){
			//printf("inside\n");
			inside=true;
		}
		//no occlusion, outside intersection
		if(!bottom && !top && !inside){
			//printf("outside\n");
			outside=true;
		}

	}
	else{
		//todo this needs to be figured out
		return NA;
	}
	
	if(both) return BOTH;
	else if(bottom) return BOTTOM;
	else if(top) return TOP;
	else if(inside) return INSIDE;
	else if(outside) return OUTSIDE;
	else if(triforce) return TRIFORCE;
	
	return NA;
		
}



int Tessellation::checkIntegrityAndBlock(CircularRegion &I, CircularRegion &J){
	if(I.circularIntersections[J.id].sasa==true){
		//printf("SANITY CHECK FAILED !!!\n");
		return -1;
	}
	else{
		//printf("blocking %d-%d\n",I.id,J.id);
		I.circularIntersections[J.id].blocked=true;
	}
	return 0;

}

void Tessellation::setTertiaryIntersection(CircularRegion &I, CircularRegion &J, CircularRegion &K, OcclusionState state){
	I.circularIntersections[J.id].tertiaryIntersections[K.id] = state;
}



list<int>::iterator incIterator(list<int>::iterator it, list<int> &sasa){
	++it;
	if(it == sasa.end()) it = sasa.begin();
	
	return it;
}




int Tessellation::whoOccludes(CircularRegion *I, CircularRegion *J, vector<CircularRegion> &circles){
	map<int,OcclusionState>::iterator it_test;
	map<int,CircularIntersection>::iterator it;
	map<int,CircularIntersection>::iterator it_j;
	map<int,CircularIntersection>::iterator it_k;
	CircularRegion *K;
	IntersectionPair ip;
	double dij,dik,djk;
	bool foundOccludingCircle;
	IntersectionPoint p;
	list<IntersectionPoint>::iterator it_l;
	int occlusion,occlusion2;
	OcclusionState occlusionState;
	int i,j,k;
	
	i = I->id;
	j = J->id;
	
	for(it_k=I->circularIntersections.begin(); it_k != I->circularIntersections.end(); ++it_k){
		if(it_k!=it_j){
			k = (*it_k).first;
			it = J->circularIntersections.find(k);
			//3body intersection
			if(it != J->circularIntersections.end()){
				K = &circles[k];
				
				it_test = I->circularIntersections[J->id].tertiaryIntersections.find(k);
				if(it_test != I->circularIntersections[J->id].tertiaryIntersections.end()){
					occlusionState = (*it_test).second;
					
					if(occlusionState == UNOBSTRUCTED){
						continue;
					}
					
				}
				
				
				
				
				dij = I->circularIntersections[j].d;
				dik = I->circularIntersections[k].d;
				occlusion = occludesForwardIntersectionPoint(*I,*J,*K,dij,dik);
				//printf("3body intersection: %d-%d-%d (%d)\n",i,j,k,occlusion);
				
				if(occlusion == BOTTOM){
					//printf("BOTTOM\n");
					if(checkIntegrityAndBlock(*K,*I)<0) exit(-1);
					if(checkIntegrityAndBlock(*I,*J)<0) exit(-1);
					if(checkIntegrityAndBlock(*J,*K)<0) exit(-1);
					
					setTertiaryIntersection(*I,*K,*J,UNOBSTRUCTED);
					setTertiaryIntersection(*K,*J,*I,UNOBSTRUCTED);
					setTertiaryIntersection(*J,*I,*K,UNOBSTRUCTED);
					
					
				}
				if(occlusion == BOTH){
					//printf("BOTH\n");
					if(checkIntegrityAndBlock(*I,*J)<0) exit(-1);
					if(checkIntegrityAndBlock(*J,*I)<0) exit(-1);
					
					setTertiaryIntersection(*K,*J,*I,UNOBSTRUCTED);
					setTertiaryIntersection(*J,*K,*I,UNOBSTRUCTED);
					setTertiaryIntersection(*I,*K,*J,UNOBSTRUCTED);
					setTertiaryIntersection(*K,*I,*J,UNOBSTRUCTED);
				}
				if(occlusion == TOP){
					//printf("TOP\n");
					if(checkIntegrityAndBlock(*I,*K)<0) exit(-1);
					if(checkIntegrityAndBlock(*K,*J)<0) exit(-1);
					if(checkIntegrityAndBlock(*J,*I)<0) exit(-1);
					
					setTertiaryIntersection(*J,*K,*I,UNOBSTRUCTED);
					setTertiaryIntersection(*K,*I,*J,UNOBSTRUCTED);
					
				}
				if(occlusion == INSIDE){
					djk = circles[j].circularIntersections[k].d;
					occlusion2 = occludesForwardIntersectionPoint(*J,*I,*K,dij,djk);
					if(occlusion2 == INSIDE){
						if(checkIntegrityAndBlock(*K,*J)<0) exit(-1);
						if(checkIntegrityAndBlock(*J,*K)<0) exit(-1);
						
					}
					else if(occlusion2 == OUTSIDE){
						setTertiaryIntersection(*J,*K,*I,UNOBSTRUCTED);
						setTertiaryIntersection(*K,*J,*I,UNOBSTRUCTED);
					}
					else exit(-1);
					
					//printf("INSIDE\n");
					if(checkIntegrityAndBlock(*I,*K)<0) foundOccludingCircle=true;
					if(checkIntegrityAndBlock(*K,*I)<0) foundOccludingCircle=true;
				}
				if(occlusion == OUTSIDE){
					//check for triforce
					djk = circles[j].circularIntersections[k].d;
					occlusion2 = occludesForwardIntersectionPoint(*J,*I,*K,dij,djk);
					if(occlusion2 == OUTSIDE){
						//printf("TRIFORCE\n");
						
					}
					else{

						//printf("OUTSIDE\n");
						if(checkIntegrityAndBlock(*K,*J)<0) foundOccludingCircle=true;
						if(checkIntegrityAndBlock(*J,*K)<0) foundOccludingCircle=true;
						
						setTertiaryIntersection(*I,*K,*J,UNOBSTRUCTED);
						setTertiaryIntersection(*K,*I,*J,UNOBSTRUCTED);
						setTertiaryIntersection(*J,*I,*K,UNOBSTRUCTED);
						
					}
				}
				
				if(occlusion == BOTTOM || occlusion == BOTH){
					//printf("does occlude\n");
					return k;
					
				}
				
			}
		}	
	}
	
	return -1;
}


void Tessellation::printSasa(list<int> &sasa){
	list<int>::iterator it_i;
	
	for(it_i = sasa.begin(); it_i != sasa.end(); ++it_i){
		if(it_i==sasa.begin())
			printf("SASA: %d",*it_i);
		else printf("-%d",*it_i);
	}
	printf("\n");
	
}


bool Tessellation::mergeCutList(int k, list<int>::iterator start, list<int> &sasa){
	start = incIterator(start,sasa);
	list<int>::iterator it;
	for(it = start; it != sasa.end(); ++it){
		printf("merge comparing %d-%d\n",*it,k);
		if(*it == k){
			sasa.erase(start,it);
			return true;
		}
	}
	for(it = sasa.begin(); it != start; ++it){
		printf("merge comparing %d-%d\n",*it,k);
		if(*it == k){
			sasa.erase(start,sasa.end());
			sasa.erase(sasa.begin(),it);
		}
		return true;
	}
	return false;
	
}


/**
 * Start with first circle and assume it belongs to the boundary of a sasa.
 * Take one of its circular intersections and put it in the list that defines this sasa.
 * check for 3 body intersections. if found, correct sasa list.
 * Loop until correct sasas have been established
 * 
 */
void Tessellation::buildIntersectionGraph(double radius, vector<CircularRegion> &circles,vector<list<IntersectionPoint*>*> &intersections){
	list<int> sasa;
	map<int,CircularIntersection>::iterator it1;
	
	list<IntersectionPoint*>* path;
	IntersectionPoint* pp;
	
	CircularRegion *I1;
	
	list<int>::iterator it_i;
	list<int>::iterator it_j;
	list<int>::iterator it_k;
	list<int>::iterator it_l;
	IntersectionPair ip;
	IntersectionPoint p;
	int i1,j1;
	int i,j,k;
	bool sasaComplete;
	int s =0;
	int i_prev;
	
	
	
	
	//determine all circle-circle intersections, so that each circle knows which other circles are intersecting with it
	determineCircularIntersections(circles);

	
	
	for(i1=0; i1 < circles.size(); ++i1){
		I1 = &circles[i1];
		printf("i1: %d (%d)\n",i1,I1->circularIntersections.size());
		for(it1=I1->circularIntersections.begin(); it1 != I1->circularIntersections.end(); ++it1){
			
			j1=(*it1).first;
			printf("j1: %d\n",j1);
			
			if(!it1->second.visited && !it1->second.blocked){
				
				sasa.clear();
				it_i = sasa.insert(sasa.end(),i1);
				i = *it_i;
				sasa.push_back(j1);
				
				it_j = incIterator(it_i,sasa);
				j = *it_j;
				
				
				
				printSasa(sasa);

				sasaComplete = false;
				while(!sasaComplete){
					printf("I: %d J: %d\n",i,j);
					if(sasa.size()<2) break; // to be inspected, seems fishy
					++s;
					if(s>100) exit(-1);
					
					
					//in this iteration we inspect I and J
					if(circles[i].circularIntersections[j].visited){
						sasaComplete=true;
						break;
					}
					circles[i].circularIntersections[j].visited=true;
					
					k = whoOccludes(&circles[i], &circles[j], circles);
					
					//if there is an occlusion, we add the occluding circle to the sasa
					if(k>=0){
						//if I-K has already been visited, we delete I
						if(circles[i].circularIntersections[k].visited){
							printf("erasing totally occluded circle\n");
							it_i = sasa.erase(it_i);
							if(it_i == sasa.end()) it_i = incIterator(it_i, sasa);
							i = *it_i;
							
						}
						else{						
							printf("adding blocking circle\n");
							it_j = sasa.insert(it_j,k);
							j = *it_j;
							
						}
					}
					else{
						//I-J is free of occlusion. We will advance I to J and select a new J
						
						
						k = selectNextIntersectingCluster(circles[i],circles[j],circles);
						printf("selected next k: %d\n",k);
						
						if(k<0){
							exit(-1);
							sasaComplete=true;
						}
						else{
							it_l = incIterator(it_j,sasa);
							
							if(k != *it_l){
								//we have to check if the newly selected j is cutting a part of the established sasa
								
								if(!mergeCutList(k,it_j,sasa)){
									//adding circle
									printf("adding new circle\n");
									it_k = sasa.insert(it_l,k);
									sasa.insert(it_l,j);
									
									it_i = it_j;
									i = j;
									it_j = it_k;
									j = k;
									
								}
								else{
									it_i = it_j;
									i=j;
									it_j = incIterator(it_i,sasa);
									j = *it_j;
									printf("merge cut\n");
								}
							}
							else{
								//if the new selected J is the next item in our already established sasa,
								//we simply advance i and j
								
								it_i = it_j;
								i = j;
								it_j = it_l;
								j = k;
								
								printf("continuing with existing circle\n");
							}
							
						}
						
					}
					
					printSasa(sasa);
					
					
				}
			
				if(sasa.size()>1){
					path = new list<IntersectionPoint*>();
					intersections.push_back(path);
					for(it_i = sasa.begin(); it_i != sasa.end(); ++it_i){
						it_j = incIterator(it_i, sasa);
						i = *it_i;
						j = *it_j;
						
						printf("found sasa ip %d-%d\n",i,j);
						//here, we add a forwardIntersection to the circle
						ip = determineIntersectionPoints(radius, circles[i], circles[j]);
						p.vector = ip.k_j;
						p.with = j;
						p.from = i;
						pp = &*circles[i].forwardIntersections.insert(circles[i].forwardIntersections.end(),p);	
						circles[i].intersect = true;	
						path->push_back(pp);
					}
				}
				
				printf("-------------------------------------------\n");
			}
		}
	}
	
}




/*
void Tessellation::buildIntersectionGraph(double radius, vector<CircularRegion> &circles){
	int i,j,k;
	map<int,OcclusionState>::iterator it_test;
	map<int,CircularIntersection>::iterator it;
	map<int,CircularIntersection>::iterator it_j;
	map<int,CircularIntersection>::iterator it_k;
	CircularRegion *I,*J,*K;
	IntersectionPair ip;
	double dij,dik,djk;
	bool foundOccludingCircle;
	IntersectionPoint p;
	list<IntersectionPoint>::iterator it_l;
	int occlusion,occlusion2;
	OcclusionState occlusionState;

	//determine all circle-circle intersections, so that each circle knows which other circles are intersecting with it
	determineCircularIntersections(circles);
	
	
	//iterate over all circles
	for(i=0; i< circles.size(); i++){
		I = &circles[i];
		printf("I: %d, id %d, size: %d\n",i, I->id, I->circularIntersections.size());
		
		
		for(it_j=I->circularIntersections.begin(); it_j != I->circularIntersections.end(); ++it_j){
			
			j=(*it_j).first;
			J=&circles[j];
			
			if(I->circularIntersections[j].blocked){
				printf("BLOCKED\n");
				continue;
			}
			
			printf("J: %d, size: %d\n",j, J->circularIntersections.size());
			
			
			foundOccludingCircle = false;
			for(it_k=I->circularIntersections.begin(); it_k != I->circularIntersections.end(); ++it_k){
				if(it_k!=it_j){
					k = (*it_k).first;
					it = J->circularIntersections.find(k);
					//3body intersection
					if(it != J->circularIntersections.end()){
						K = &circles[k];
						
						it_test = I->circularIntersections[J->id].tertiaryIntersections.find(k);
						if(it_test != I->circularIntersections[J->id].tertiaryIntersections.end()){
							occlusionState = (*it_test).second;
							
							if(occlusionState == UNOBSTRUCTED){
								continue;
							}
							
						}
						
						
						
						
							dij = I->circularIntersections[j].d;
							dik = I->circularIntersections[k].d;
							occlusion = occludesForwardIntersectionPoint(*I,*J,*K,dij,dik);
							printf("3body intersection: %d-%d-%d (%d)\n",i,j,k,occlusion);
							if(occlusion == BOTTOM || occlusion == BOTH){
								foundOccludingCircle = true;
								printf("does occlude\n");
								
							}
							
							if(occlusion == BOTTOM){
								printf("BOTTOM\n");
								if(checkIntegrityAndBlock(*K,*I)<0) exit(-1);
								if(checkIntegrityAndBlock(*I,*J)<0) exit(-1);
								if(checkIntegrityAndBlock(*J,*K)<0) exit(-1);
								
								setTertiaryIntersection(*I,*K,*J,UNOBSTRUCTED);
								setTertiaryIntersection(*K,*J,*I,UNOBSTRUCTED);
								setTertiaryIntersection(*J,*I,*K,UNOBSTRUCTED);
								
								
							}
							if(occlusion == BOTH){
								printf("BOTH\n");
								if(checkIntegrityAndBlock(*I,*J)<0) exit(-1);
								if(checkIntegrityAndBlock(*J,*I)<0) exit(-1);
								
								setTertiaryIntersection(*K,*J,*I,UNOBSTRUCTED);
								setTertiaryIntersection(*J,*K,*I,UNOBSTRUCTED);
								setTertiaryIntersection(*I,*K,*J,UNOBSTRUCTED);
								setTertiaryIntersection(*K,*I,*J,UNOBSTRUCTED);
							}
							if(occlusion == TOP){
								printf("TOP\n");
								if(checkIntegrityAndBlock(*I,*K)<0) exit(-1);
								if(checkIntegrityAndBlock(*K,*J)<0) exit(-1);
								if(checkIntegrityAndBlock(*J,*I)<0) exit(-1);
								
								setTertiaryIntersection(*J,*K,*I,UNOBSTRUCTED);
								setTertiaryIntersection(*K,*I,*J,UNOBSTRUCTED);
								
							}
							if(occlusion == INSIDE){
								djk = circles[j].circularIntersections[k].d;
								occlusion2 = occludesForwardIntersectionPoint(*J,*I,*K,dij,djk);
								if(occlusion2 == INSIDE){
									if(checkIntegrityAndBlock(*K,*J)<0) exit(-1);
									if(checkIntegrityAndBlock(*J,*K)<0) exit(-1);
									
								}
								else if(occlusion2 == OUTSIDE){
									setTertiaryIntersection(*J,*K,*I,UNOBSTRUCTED);
									setTertiaryIntersection(*K,*J,*I,UNOBSTRUCTED);
								}
								else exit(-1);
								
								printf("INSIDE\n");
								if(checkIntegrityAndBlock(*I,*K)<0) foundOccludingCircle=true;
								if(checkIntegrityAndBlock(*K,*I)<0) foundOccludingCircle=true;
							}
							if(occlusion == OUTSIDE){
								//check for triforce
								djk = circles[j].circularIntersections[k].d;
								occlusion2 = occludesForwardIntersectionPoint(*J,*I,*K,dij,djk);
								if(occlusion2 == OUTSIDE){
									printf("TRIFORCE\n");
									
								}
								else{

									printf("OUTSIDE\n");
									if(checkIntegrityAndBlock(*K,*J)<0) foundOccludingCircle=true;
									if(checkIntegrityAndBlock(*J,*K)<0) foundOccludingCircle=true;
									
									setTertiaryIntersection(*I,*K,*J,UNOBSTRUCTED);
									setTertiaryIntersection(*K,*I,*J,UNOBSTRUCTED);
									setTertiaryIntersection(*J,*I,*K,UNOBSTRUCTED);
									
								}
							}

							if(occlusion == NA){
								printf("NA\n");
								exit(-1);
							}
							
							
						//}

						
					}
				}
				
				
			}
			
			
			if(!foundOccludingCircle){
				printf("found sasa ip %d-%d\n",i,j);
				//here, we add a forwardIntersection to the circle
				ip = determineIntersectionPoints(radius, *I, *J);
				p.vector = ip.k_j;
				p.with = j;
				p.from = i;
				circles[i].forwardIntersections.push_back(p);	
				circles[i].intersect = true;
				I->circularIntersections[j].sasa=true;
			}			
			
		}
	}
	
	
	for(i=0; i< circles.size(); i++){
		printf("circle %d\n",i);
		for(it_l=circles[i].forwardIntersections.begin(); it_l!=circles[i].forwardIntersections.end();++it_l){
			printf("ip: %d-%d\n",it_l->from, it_l->with);
		}
		
	}
	
			
	

}
*/

void Tessellation::buildGaussBonnetPath(Vector &origin, double radius, vector<Vector> &atoms, vector<double> &radii, vector<CircularRegion> **circles, vector<list<IntersectionPoint*>*>** intersections){
	IntersectionPoint p;
	int i,k,j;
	IntersectionPair points;
	double angle;
	Vector v,v2,o,o2;
	
	*circles = new vector<CircularRegion>();
	*intersections = new vector<list<IntersectionPoint*>*>();
	srand(2);


	makeCircularRegions(origin, radius, atoms, radii, **circles);
	for(i=0;i<(*circles)->size();i++){
		determineProjection(origin, radius, (*circles)->at(i));
	}
	filterCircularRegions(radius, **circles);
	
	clearFlags(**circles);
	reindexCircularRegions(**circles);
	
	buildIntersectionGraph(radius, **circles, **intersections);

	p.visited = false;
	clearFlags(**circles);
	
	//NSWarnLog(@"gauss bonnet path is ready...");
	
	insertFakeIntersectionPoints(**circles);
	
	
}




void Tessellation::harvestIntersectionPoints(vector<CircularRegion> &circles, vector<vec> &intersections){
	int i;
	list<IntersectionPoint>::iterator it;

	for(i=0;i<circles.size();i++)
		for(it = circles[i].forwardIntersections.begin(); it != circles[i].forwardIntersections.end(); ++it){
			//NSWarnLog(@"POINT: [%f %f %f]",it->vector.vector[0],it->vector.vector[1],it->vector.vector[2]);
			intersections.push_back(it->vector);
			}
}


bool Tessellation::hasUnflaggedIntersectionPoints(CircularRegion &circle, IntersectionPoint **ip){
	list<IntersectionPoint>::iterator it;
	for(it = circle.forwardIntersections.begin(); it != circle.forwardIntersections.end(); ++it){
		if(!it->visited){
			*ip=&*it;
			return true;
		}
	}
	return false;
}


list<IntersectionPoint*>* Tessellation::retrieveIntersections(CircularRegion &circle){
	list<IntersectionPoint*>* res;
	res = new list<IntersectionPoint*>();
	
	list<IntersectionPoint>::iterator it;
	for(it = circle.forwardIntersections.begin(); it != circle.forwardIntersections.end(); ++it){
		//if(!it->visited){
			//save a pointer to the intersectionpoint
			res->push_back(&(*it));
		//}
	}
	return res;
	
}


void Tessellation::showIntersections(list<IntersectionPoint*> &intersections){
	list<IntersectionPoint*>::iterator it;
	
	for(it=intersections.begin();it!=intersections.end();++it)
		fprintf(stdout,"INTERSECTION %d-%d",(*it)->from, (*it)->with);
}

/*
 * path and res needs to be deleted somewhere
 */
vector<list<IntersectionPoint*>*>*  Tessellation::harvestGaussBonnetPaths(vector<CircularRegion> &circles){
	vector<list<IntersectionPoint*>*>* res;
	list<IntersectionPoint*>* path;
	list<IntersectionPoint*>* intersectionsTest;
	list<IntersectionPoint*>* intersections;
	list<IntersectionPoint*>::iterator it_i, minIt;

	Vector v,ortho,test,reference;
	IntersectionPoint  *p, *p_prev;
	CircularRegion circle;
	bool done;
	double minAngle,angle;
	int i,j;
	
	res = new vector<list<IntersectionPoint*>*>();
	path = new list<IntersectionPoint*>();
	//the whole thing can actually be emtpy, in that case we better get outta here
	if(circles.size()==0) return res;
	//outputCircularRegions(circles);
	
	clearFlags(circles);
	
	
	for(j=0;j<circles.size();++j){
		if(hasUnflaggedIntersectionPoints(circles[j],&p_prev)){
			//printf("USING SEED %d, %d-%d\n",j, p_prev->from, p_prev->with);
			//p_prev->visited=true;
			done=false;
			while(!done){
				//NSWarnLog(@"STARTING ROUND %d",p_prev->with);
				circle = circles[p_prev->with];
				//NSWarnLog(@"GETTING INTERSECTIONS FOR: %d",p_prev->with);
				intersections = retrieveIntersections(circles[p_prev->with]);
				
				if(intersections->size() > 1){
					//showIntersections(*intersections);
					//NSWarnLog(@"MULTIINTERSECTIONS");
					//determine which point we need to take, using p_prev as a reference
					
					//first we draw an orthogonal vector to the plane of the circular region
					v = circle.normal * circle.g;
					//from there we draw a vector to the last visited intersectionpoint
					reference = p_prev->vector - v;
		
					//we create a plane that lies in the previously constructed vector and the plane-orthogonal vector.
					//ortho is then the normal vector of the plane
					ortho = cross(p_prev->vector,v);
					ortho = ortho / norm(ortho,2);
		
					minAngle = 2 * M_PI;
					for(it_i=intersections->begin(); it_i!=intersections->end(); ++it_i){
						
						//NSWarnLog(@"vector: (%f %f %f)",intersections[i]->vector.vector[0],intersections[i]->vector.vector[1],intersections[i]->vector.vector[2]);
						//we draw a vector from the plane-orthogonal vector to the currently inspected intersectionpoint
						test = (*it_i)->vector - v;
						angle = getAngle(reference,test);
						
						double oldAngle = angle;
						//here, we check on which side of the plane the intersection point is
						//(if it is counterclockwise or clockwise oriented regarding to the reference)
						if(circle.form == CONVEX){
							if(dot(ortho,test) < 0) angle = 2 * M_PI - angle;
						}
						else{
							if(dot(ortho,test) > 0) angle = 2 * M_PI - angle;
						}
						//NSWarnLog(@"EVAL: %f %f - %f",angle,minAngle,Ad3DDotProduct(&ortho,&test));
						
						if(isZero(angle)) angle = 2 * M_PI;
						
						if(angle < minAngle){
							minAngle = angle;
							minIt = it_i;
						}
						
						//NSWarnLog(@"INT: %d-%d ANG: %f / %f [%f] %f %d",(*it_i)->from,(*it_i)->with,angle,oldAngle,minAngle,Ad3DDotProduct(&ortho,&test),circle.form);
						
					}
		
					//we got the next point :)
					
					p = *minIt;
					if(!p->visited){
						p->visited = true;
						p->id=path->size();
						p_prev = p;
						path->push_back(p);
					}
					else done=true;
					
				}
				else{
					//NSWarnLog(@"NORMAL INTERSECTION");
					p = (*intersections->begin());
					if(!p->visited){
						p->visited = true;
						p->id = path->size();
						p_prev = p;
						path->push_back(p);
					}
					else done=true;
				}
				
			
				delete intersections;
			}
			res->push_back(path);
			path = new list<IntersectionPoint*>();
		}
		delete intersectionsTest;
	}
	return res;
}




/*
Vector Tessellation::getSphericalGreatArcGreatArcInterfacePoint(Vector &a, Vector &b, double radius){
	Vector f,ab,en,n;
	double s;
	
	//NSWarnLog(@"sn : (%f, %f, %f)",n->vector[0],n->vector[1],n->vector[2]);
	//NSWarnLog(@"sa : (%f, %f, %f)",a->vector[0],a->vector[1],a->vector[2]);
	//NSWarnLog(@"sb : (%f, %f, %f)",b->vector[0],b->vector[1],b->vector[2]);
	
	//create edge-plane
	ab = b-a;
	en = cross(a,ab);
	
	//the required interface is perpendicular to the two plane normals. Since the two planes also cross in the origin by definition, we
	//immediately have found the two interfacepoints. We select the right one by considering the dot product of the normals 
	n=cross(en,f);
	s = vsign(dot(f,a));

	flen = norm(f,2);
	f = f * (s*radius/flen);

	return f;
}
*/

	
/*
CircularRegion getGreatOpeningCircle(double radius, Vector3D *K, double r_k, int index){
	double g,a,side,alpha;
	CircularRegion res;
	Vector3D normal;

	Ad3DVectorLength(K);
	//check if neighbor is completely inside estimate:
	if(radius >= K->length + r_k){
		alpha = 0;
		a = 0;
		g = 0;
		//NSWarnLog(@"SPLITTER CASE 1 %f, %f, %f",radius, K->length, r_k);
	}
	else
	//check if our origin-sphere is completely inside K!!! this happens practically for all hydrogens and it's horrible!!!
	if(r_k >= K->length + radius)
	{
		alpha = M_M_PI;
		g = 0;
		a = 0;
		//NSWarnLog(@"SPLITTER CASE 2");
	}
	else
	{
		if(sqrt(pow(radius,2) + pow(r_k,2)) > K->length)
		{
			g = (K->length * K->length - r_k * r_k + radius * radius) /  (2 * K->length);
			a = (1.0/(2.0*K->length)) * sqrt((-K->length+r_k-radius)*(-K->length-r_k+radius)*(-K->length+r_k+radius)*(K->length+r_k+radius));
			//side = sqrt(g * g + a * a);

			if(g<0)	alpha = M_M_PI - asin(a/radius);
			else alpha = asin(a/radius);
		//NSWarnLog(@"SPLITTER CASE 3");
		}
		else
		{
			alpha = asin(r_k/K->length);
			g = cos(alpha) * radius;
			a = sin(alpha) * radius;

		//NSWarnLog(@"SPLITTER CASE 4");
		}
	}

	res.openingAngle=alpha;
	res.g = g;
	res.a = a;
	res.vector = *K;
	Ad3DVectorLength(&res.vector);
	Ad3DVectorScalarMultiply2(K, 1.0/K->length, &normal);
	res.normal = normal;
	Ad3DVectorLength(&res.normal);
	res.index=index;
	//probably not useful
	if(g<0) res.form=CONCAVE;
	else res.form=CONVEX;
	//this particular sphereradius is just a construct for later calculations and does represent anything physical 
	res.sphereRadius = sqrt((K->length-g)*(K->length-g) + a * a);

	return res;
}
*/


void Tessellation::outputGaussBonnetData(string filename){
	FILE* file;
	
	file = fopen (filename.c_str(),"w");
	
	vector<list<IntersectionPoint*>*>* sasas;
	list<IntersectionPoint*>* sasa;
	list<IntersectionPoint*>::iterator it;
	vector<CircularRegion>* crs;
	CircularRegion circle;
	Vector integrationOrigin;
	double area=0;
	int k;
	
	
	//iterate over all atoms
	for(int i=0;i<intersectionSet.size();++i){
		//iterate over all sasas
		sasas = intersectionSet.at(i);
		printf("size: %d\n",sasas->size());
		crs = circleSet.at(i);

		fprintf(file,"atom %d radius %f crd %f %f %f\n", i, radii->at(i), atoms[i](0), atoms[i](1), atoms[i](2));
		for(int j=0; j<atoms.size(); j++)
			if(j!=i) fprintf(file, "neighbor %d radius %f vector %f %f %f\n",i,radii->at(j), atoms[j](0)-atoms[i](0), atoms[j](1)-atoms[i](1), atoms[j](2)-atoms[i](2));
		
		
		for(int j=0;j<crs->size();++j){
			circle = crs->at(j);
			fprintf(file, "circularregion %d radius %f vector %f %f %f form %d\n", j, circle.a, circle.normal(0)*circle.g, circle.normal(1)*circle.g, circle.normal(2)*circle.g, circle.form);
			fprintf(file, "intersector %d radius %f vector %f %f %f\n", j, circle.sphereRadius, circle.vector(0), circle.vector(1), circle.vector(2));
		}
			
		for(int j=0;j<sasas->size();++j){
			sasa = sasas->at(j);
			fprintf(file, "sasa %d size %d\n", j, sasa->size());
			k=0;
			for(it=sasa->begin();it!=sasa->end();++it,++k){
				fprintf(file, "intersectionpoint %d vector %f %f %f\n", k, (*it)->vector(0), (*it)->vector(1), (*it)->vector(2));
			}
		}
		
	}
	
		
	

}
