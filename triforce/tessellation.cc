#include "tessellation.h"


Tessellation::Tessellation(Molecule &m){
	molecule = m;
}



void Tessellation::build(){
	CircularRegionsPerAtom circlesPerAtom;
	sasasForMolecule.clear();
	
	//molecule.update();
	atoms = molecule.coordinates();
	radii = molecule.fetchRadii();
	
	
	
	//iterate over all atoms and build the tessellation for each of them
	//for(int i=0; i<atoms.size(); ++i){
		int i=0;
		buildGaussBonnetPath(atoms[i], radii->at(i), atoms, *radii, sasasForMolecule);
		
		
		
		for(int i=0; i<sasasForMolecule[0].sasas.size(); ++i){
			printf("SASA %d\n",i);
			outputGaussBonnetPath(sasasForMolecule[0].sasas[i]);
		}


		
		//outputCircularRegions(*circles[circles.size()-1]);

	//}
	
	
	
}


CircularRegionsPerAtom Tessellation::coverHemisphere(Vector tessellationOrigin, double radius, CircularRegionsPerAtom circles){
	CircularRegion C;
	Vector v(3);
	
	v = -tessellationOrigin;
	
	C.openingAngle = M_PI * 0.5;
	C.a = radius;
	C.g=0;
	C.vector = v;
	C.normal = v;
	C.form = SPLITTER;
	C.id = circles.size();
	
	circles.push_back(C);
	
	return circles;
	
}



void Tessellation::buildGaussBonnetPath(Vector &origin, double radius, vector<Vector> &atoms, vector<double> &radii, SASAsForMolecule &sasas){
	int i,k,j;
	CircularRegionsPerAtom circles;
	CircularRegionsPerAtom circlesFrontHemisphere;
	CircularRegionsPerAtom circlesBackHemisphere;
	
	Vector frontTessellationOrigin(3);
	Vector backTessellationOrigin(3);
	
	
	frontTessellationOrigin(0) = 1;
	frontTessellationOrigin(1) = 0;
	frontTessellationOrigin(2) = 0;

	backTessellationOrigin = -frontTessellationOrigin;
	SASAs *newSasas;
	
	srand(2);
	
	SASAsForAtom sasasForAtom;
	sasasForAtom.radius = radius;
	sasasForAtom.vector = origin;
	
	newSasas = &(sasas.insert(sasas.end(),sasasForAtom)->sasas);
	


	makeCircularRegions(origin, radius, atoms, radii, circles);
	for(i=0;i<circles.size();i++){
		determineProjection(origin, radius, circles[i]);
	}
	
	//there is room for optimisation here...
	
	circlesFrontHemisphere = coverHemisphere(frontTessellationOrigin, radius, circles);
	circlesBackHemisphere = coverHemisphere(backTessellationOrigin, radius, circles);
	
	
	filterCircularRegions(radius, circlesFrontHemisphere);
	filterCircularRegions(radius, circlesBackHemisphere);
	reindexCircularRegions(circlesFrontHemisphere);
	reindexCircularRegions(circlesBackHemisphere);
	
	determineCircularIntersections(circlesFrontHemisphere);
	determineCircularIntersections(circlesBackHemisphere);
	
	printf("CIRC: %d, FRONT: %d, BACK: %d\n",circles.size(), circlesFrontHemisphere.size(), circlesBackHemisphere.size());
	
	buildIntersectionGraph(radius, frontTessellationOrigin, circlesFrontHemisphere, *newSasas,string("gbonnet0.csv"));
	buildIntersectionGraph(radius, backTessellationOrigin, circlesBackHemisphere, *newSasas,string("gbonnet1.csv"));
	
}






SASAsForMolecule &Tessellation::sasas(){
	return sasasForMolecule;
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

bool Tessellation::isInPositiveEpsilonRange(double v, double eps){
	printf("V: %f, EPS: %f, V+T: %f, EPS-V+T: %f\n",v,eps,v+THRESHOLD_NUMERICAL, eps-(v+THRESHOLD_NUMERICAL));
	if(eps-(v+THRESHOLD_NUMERICAL) <= 0) return true;
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
	
	if(g<0){
		circle.form=CONCAVE;
		g = abs(g);
	}
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




void Tessellation::outputGaussBonnetPath(SASA &points){
	SASANodeList::iterator it;
	int i;
	
	for(it=points.sasa.begin(), i=0; it!=points.sasa.end(); ++it, ++i)
		printf("GBPATH[%d] %d - %d\n", i, it->id0, it->id1);
}



void Tessellation::reindexCircularRegions(CircularRegionsPerAtom &circles){
	int i,j;
	int currentID;
	for(i=0;i<circles.size();i++){
		currentID = circles[i].id;
		circles[i].id = i;
	}
}




void Tessellation::insertFakeIntersectionPoints(vector<CircularRegion> &circles){
	Vector v,v2,o,o2;
	int k;
	//IntersectionPoint p;
/*	
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
*/	
}

int Tessellation::sgn(double d){
	if(d>=0) return 1;
	else return -1;
}

void Tessellation::determineCircularIntersections(CircularRegionsPerAtom &circles){
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
	
	//addition
	if(v<0) v = 2*M_PI+v;
	   
	return v;
		
}




/*
 * 
[7] (7,1) d: 0.725943 (OUT)
[7] (5,7) d: 0.769139 (IN)
[7] (7,3) d: 1.506200 (OUT)
[7] (7,4) d: 2.223748 (OUT)
[7] (4,7) d: 2.265115 (IN)
[7] (3,7) d: 3.938589 (IN)
[7] (1,7) d: 4.631484 (IN)
[7] (7,5) d: 4.641597 (OUT)

 * 
 * */


double Tessellation::angularInterface(Vector &x0, Vector &v, Vector &p0, Vector &p1){
	Vector vx(3);
	Vector vp0(3);
	Vector vp1(3);
	double eta;
	int s;
	
	
	
	
	vx = x0-v;
	vp0 = p0-v;
	vp1 = p1-v;

	eta = getAngle(vp1,vx);
	
	
	s = sgn(dot(vp0,vx));
	
	if(s > 0) eta = -eta;
	
	return eta;
	
}

void Tessellation::measurementPoints(Vector &p0, Vector &p1, Vector &tessellationOrigin, CircularRegion &I){
	Vector v(3);
	Vector o(3);
	Vector vx(3);
	Vector s(3), s2(3);
	
	o = tessellationOrigin;
	v = I.normal * I.g;
	
	if(isInPositiveEpsilonRange(fabs(dot(o,I.normal)),1.0) ){
		printf("++++++++++++++++++++++++DOTDOTDOT+++++++++++++++++++\n");
		p1 = Vector(3);
		p1(0) = 0;
		p1(1) = 1;
		p1(2) = 0;
		
		p0 = cross(p1,o);
	}
	else{
		s = cross(v,o);
		s = I.a*s / norm(s,2);
		p0 = v+s;
		
		s2 = -cross(v,s);
		s2 = I.a*s2 / norm(s2,2);
		p1 = v+s2;
	}
	
}

Interfaces Tessellation::angularInterfaces(Vector &x0, Vector &x1, Vector &tessellationOrigin, CircularRegion &I){
	Vector p0(3), p1(3);
	Interfaces res;
	Vector v(3);
	
	//measurement points will be written into p0 and p1
	measurementPoints(p0,p1,tessellationOrigin,I);
	
	
	v = I.normal * I.g;
	
	
	res.out = angularInterface(x0,v,p0,p1);
	res.vectorOut = x0;
	res.in = angularInterface(x1,v,p0,p1);
	res.vectorIn = x1;
	
	//if(I.form == CONCAVE){
	//	res.out = -res.out;
	//	res.in = -res.in;
	//}
	
	return res;
	
}

Interfaces Tessellation::retrieveInterfaces(Vector tessellationOrigin, CircularRegion &I, CircularRegion J, double dij, double radius){
	double gij;
	double eta;
	double mu;
	double entryPoint, exitPoint;
	IntersectionPair ip;
	Vector x0(3), x1(3);
	
	
	
	
	ip = determineIntersectionPoints(radius, I, J);
	x0 = ip.k_j;
	x1 = ip.j_k;
	
	
	return angularInterfaces(x0,x1, tessellationOrigin, I);
	
	
	
	
}






IntersectionBranches::iterator Tessellation::increaseBranchInterator(IntersectionBranches::iterator it){
	IntersectionBranches* p;
	p = it->second.body;
	++it;
	if(it == p->end()) it=p->begin();
	
	return it;
	
}

IntersectionBranches::iterator Tessellation::decreaseBranchInterator(IntersectionBranches::iterator it){
	IntersectionBranches* p;
	p = it->second.body;
	if(it == p->begin()) it=p->end();
	
	
	--it;
	
	return it;
}



IntersectionBranches::iterator Tessellation::increaseBranchInterator(multimap<double, IntersectionBranch>::iterator it, int ignore){
	IntersectionBranches* p;
	p = it->second.body;
	++it;
	if(it == p->end()) it=p->begin();
	
	if(it->second.it->second.id == ignore) it = increaseBranchInterator(it, ignore);
	
	return it;
	
}

IntersectionBranches::iterator Tessellation::decreaseBranchInterator(IntersectionBranches::iterator it, int ignore){
	IntersectionBranches* p;
	p = it->second.body;
	if(it == p->begin()) it=p->end();
	
	--it;
	
	//printf("chk %d against %d\n",it->second.it->second.id, ignore);
	if(it->second.it->second.id == ignore) it = decreaseBranchInterator(it, ignore);
	
	
	return it;
}



void Tessellation::disconnectIntersectionPoint(IntersectionNode &a){
	a.pointsTo0 = -1;
	a.pointsTo1 = -1;
}
void Tessellation::connectIntersectionPoints(IntersectionNode &a, IntersectionNode &b){
	a.pointsTo0 = b.id0;
	a.pointsTo1 = b.id1;
}



void Tessellation::createIntersectionNode(IntersectionAddress &address, IntersectionGraph &intersectionGraph){
	IntersectionNode node;
	
	node.id0 = address.id0;
	node.id1 = address.id1;
	node.pointsTo0=-1;
	node.pointsTo1=-1;
	node.visited=false;
	intersectionGraph[address] = node;
}


void Tessellation::createIntersectionNode(int id0, int id1, IntersectionGraph &intersectionGraph){
	IntersectionAddress address;
	
	address.id0 = id0;
	address.id1 = id1;
	
	createIntersectionNode(address, intersectionGraph);
	
}





void Tessellation::createIntersectionBranch(IntersectionAddress &address, Interfaces interfacesI, Interfaces interfacesJ, CircularRegion &I, CircularRegion &J, IntersectionGraph &intersectionGraph){
	
	pair<double, IntersectionBranch> x;
	IntersectionBranches::iterator it0;
	IntersectionBranches::iterator it1;
	
	x.second.visited = false;
	
	x.first = interfacesJ.in;
	x.second.node=&intersectionGraph[address];
	x.second.direction = IN;
	x.second.it = I.intersectionBranches.end();
	x.second.body = &J.intersectionBranches;
	x.second.id = address.id1;
	it0 = J.intersectionBranches.insert(x);
	
	printf("address %d-%d side %d IN\n",address.id0, address.id1, address.id1);
	
	x.first = interfacesI.out;
	x.second.node=&intersectionGraph[address];
	x.second.direction = OUT;
	x.second.it = it0;
	x.second.body = &I.intersectionBranches;
	x.second.id = address.id0;
	it1 = I.intersectionBranches.insert(x);
	it0->second.it=it1;

	printf("address %d-%d side OUT\n",address.id0, address.id1, address.id0);
	
	
	intersectionGraph[address].angle0=interfacesI.out;
	intersectionGraph[address].vector=interfacesI.vectorOut;
	intersectionGraph[address].angle1=interfacesJ.in;
	intersectionGraph[address].vector=interfacesJ.vectorIn;
	
}


void Tessellation::printBranch(const char* s, multimap<double, IntersectionBranch>::iterator &it){
	if(it->second.direction == OUT)
		printf("iterator %s: %d-%d (OUT)\n",s,it->second.node->id0, it->second.node->id1);
	else
		printf("iterator %s: %d-%d (IN)\n",s,it->second.node->id0, it->second.node->id1);
}

void Tessellation::printIntersectionGraph(IntersectionGraph &g){
	IntersectionGraph::iterator it;
	
	for(it = g.begin(); it != g.end(); ++it){
		printf("NODE[%d,%d]: (%d,%d) -> (%d,%d)\n", it->first.id0, it->first.id1, it->second.id0, it->second.id1, it->second.pointsTo0, it->second.pointsTo1);
	}
}


void Tessellation::buildIntersectionGraph(double radius, Vector &tessellationOrigin, CircularRegionsPerAtom &circles, SASAs &sasas, string filename){
	map<int, bool> processed;
	map<int, bool>::iterator it_p;
	IntersectionGraph intersectionGraph;
	IntersectionGraph::iterator it_g, it_x;
	
	SASA potentialSasa;
	SASANode sasaNode;
	
	int i,j;
	CircularRegion *I, *J;
	double tau0, tau1;
	map<int,CircularIntersection>::iterator it_j;
	
	bool empty;
	IntersectionAddress start, x, t;
	bool valid;
	
	IntersectionPair ip;
	
	
	Interfaces interfacesJ, interfacesI;
	IntersectionBranches::iterator it_main;
	IntersectionBranches::iterator it0, it1;
	IntersectionBranches::iterator it, it_mirror, it_next, it_prev, it_mirror_next, it_mirror_prev, it_mirror_next_ignore, it_mirror_prev_ignore;
	IntersectionAddress addressIJ, addressJI;
	
	
	
	
	//iterate through all circles and add them to the intersectiongraph, one by one
	for(i=0; i < circles.size(); ++i){
		
		printf("adding %d\n",i);
		I = &circles[i];
		processed[i]=true;
		
		
		
		
		
		//we have to sort this circle into the intersectiongraph
		//go through all intersecting circles, check whether they have been processed, and if so, calculate the intersections
		
		for(it_j=I->circularIntersections.begin(); it_j != I->circularIntersections.end(); ++it_j){
			j=(*it_j).first;
			J = &circles[j];
			
			
			
			
			
			it_p = processed.find(j);
			if(it_p != processed.end()){
				printf("processing intersection with %d\n",j);
				
				
				addressIJ.id0 = i;
				addressIJ.id1 = j;
				addressJI.id0 = j;
				addressJI.id1 = i;
				
				//push the intersection points
				createIntersectionNode(addressIJ,intersectionGraph);
				createIntersectionNode(addressJI,intersectionGraph);
				
				
				//retrieve external and internal interfaces (respectively)
				interfacesJ = retrieveInterfaces(tessellationOrigin, *J, *I, it_j->second.d, radius);
				interfacesI = retrieveInterfaces(tessellationOrigin, *I, *J, it_j->second.d, radius);
				
				createIntersectionBranch(addressIJ, interfacesI, interfacesJ, *I, *J, intersectionGraph);
				createIntersectionBranch(addressJI, interfacesJ, interfacesI, *J, *I, intersectionGraph);
				
			}
		}
		
		for(int m=0; m<circles.size(); ++m){
			printf("circle %d\n",m);
			for(it_main = circles[m].intersectionBranches.begin(); it_main != circles[m].intersectionBranches.end(); ++it_main){
				if(it_main->second.direction==OUT)
					printf("[%d] (%d,%d) d: %f - (%f,%f) (OUT)\n",it_main->second.id, it_main->second.node->id0, it_main->second.node->id1, it_main->first, it_main->second.node->angle0, it_main->second.node->angle1);
				else
					printf("[%d] (%d,%d) d: %f - (%f,%f) (IN)\n",it_main->second.id, it_main->second.node->id0, it_main->second.node->id1, it_main->first, it_main->second.node->angle0, it_main->second.node->angle1);
			}
		}
		
			
			
		//all intersectionpoints have been added, it is time to change topologies
		//start at one of I's intersectionbranches
		
		for(it_main = I->intersectionBranches.begin(); it_main != I->intersectionBranches.end(); ++it_main){
			it_main->second.visited=false;
		}
		
		for(it_main = I->intersectionBranches.begin(); it_main != I->intersectionBranches.end(); ++it_main){
			//if(!it_main->second.visited)
			{
				it = it_main;
				printf("investigating starts\n");
				
				//first, decide if it is an internal or external intersection point.
				//internal points have an IN anchestor and an OUT successor.
				//external points have an OUT anchestor and an IN successor
				
				it_mirror = it->second.it;
				if(it_mirror->second.body->size()<=2) empty=true;
				else empty=false;
				
				if(!empty){
					it_mirror_next_ignore = increaseBranchInterator(it_mirror, it->second.id);
					it_mirror_prev_ignore = decreaseBranchInterator(it_mirror, it->second.id);
				}
				else{
					printf("EMPTY\n");
					it_mirror_next_ignore = increaseBranchInterator(it_mirror);
					it_mirror_prev_ignore = decreaseBranchInterator(it_mirror);
				}

				it_mirror_next = increaseBranchInterator(it_mirror);
				it_mirror_prev = decreaseBranchInterator(it_mirror);
				
				it_next = increaseBranchInterator(it);
				it_prev = decreaseBranchInterator(it);
				
				
				printBranch("it",it);
				printBranch("it_mirror",it_mirror);
				printBranch("it_mirror_next",it_mirror_next);
				printBranch("it_mirror_prev",it_mirror_prev);
				printBranch("it_mirror_next_ignore",it_mirror_next_ignore);
				printBranch("it_mirror_prev_ignore",it_mirror_prev_ignore);
				printBranch("it_next",it_next);
				printBranch("it_prev",it_prev);
				
				
				//external
				if(it_mirror_prev_ignore->second.direction == IN && it_mirror_next_ignore->second.direction == OUT){
					//if((circles[it->second.id].form==CONVEX && it->second.direction == IN) || (circles[it->second.id].form==CONCAVE && it->second.direction == OUT)){
					if(it->second.direction == IN){
						printf("EXTERNAL IN\n");
						connectIntersectionPoints(*it_mirror_prev->second.node, *it_mirror->second.node);
						connectIntersectionPoints(*it_mirror->second.node, *it_next->second.node);
						
						it_mirror_prev->second.visited=true;
						
					}
					else{
						printf("EXTERNAL OUT\n");
						connectIntersectionPoints(*it_mirror->second.node, *it_mirror_next->second.node);
						connectIntersectionPoints(*it_prev->second.node, *it->second.node);
						if(!empty && !it_mirror_prev_ignore->second.visited) disconnectIntersectionPoint(*it_mirror_prev_ignore->second.node);
						
						it_mirror_prev_ignore->second.visited=true;
					}
				}
				//internal
				else{
					if(it->second.direction == IN){
						printf("INTERNAL IN\n");
						disconnectIntersectionPoint(*it_mirror_next->second.node);
					}
					else{
						printf("INTERNAL OUT\n");
						disconnectIntersectionPoint(*it_mirror_prev->second.node);
					}
						
				}
				
				
				printf("-------------------------------------------------\n");
				printIntersectionGraph(intersectionGraph);
				printf("-------------------------------------------------\n");
				
			}
			
		}
		
	}
	printf("GENERATING SASAS\n");
	
	int s = 0;
	
	for(it_g = intersectionGraph.begin(); it_g != intersectionGraph.end(); ++it_g){
		if(!it_g->second.visited){
			start.id0 = it_g->second.id0;
			start.id1 = it_g->second.id1;
			
			x.id0 = start.id0;
			x.id1 = start.id1;
			
			printf("STARTING WITH: %d-%d\n",x.id0,x.id1);
			
			
			valid=true;
			potentialSasa.sasa.clear();
			do{
				
				if(intersectionGraph[x].visited || intersectionGraph[x].pointsTo0<0 || intersectionGraph[x].pointsTo1<0){
					printf("BREAK\n");
					valid = false;
					break;
					
				}
				
				intersectionGraph[x].visited=true;
				//if(s > 20) exit(-1);
				++s;
				
				sasaNode.id0 = x.id0;
				sasaNode.id1 = x.id1;
				sasaNode.angle0 = intersectionGraph[x].angle0;
				sasaNode.angle1 = intersectionGraph[x].angle1;
				sasaNode.lambda = circles[x.id0].openingAngle;
				sasaNode.vector = intersectionGraph[x].vector;
				sasaNode.normalForCircularRegion = circles[x.id0].normal;
				sasaNode.form = circles[x.id0].form;
				
				potentialSasa.sasa.push_back(sasaNode);
				
				t.id0 = intersectionGraph[x].pointsTo0;
				t.id1 = intersectionGraph[x].pointsTo1;
				
				x = t;
					
				printf("NEXT: %d-%d\n",x.id0,x.id1);
				
				
			}
			while(!(intersectionGraph[x].id0 == start.id0 && intersectionGraph[x].id1 == start.id1));
			
			if(valid){
				potentialSasa.tessellationOrigin = tessellationOrigin;
				sasas.push_back(potentialSasa);
			}
		}
	}
	


	outputGaussBonnetData(filename, circles, sasas, intersectionGraph);

	
	
}












void Tessellation::outputGaussBonnetData(string filename, CircularRegionsPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph){
	FILE* file;
	
	file = fopen (filename.c_str(),"w");
	
	SASA sasa;
	SASANodeList::iterator it;
	CircularRegion circle;
	double area=0;
	int k;
	
	
	for(int i=0; i< circles.size(); ++i){
		circle = circles[i];
		fprintf(file, "circularregion %d radius %f vector %f %f %f form %d\n", i, circle.a, circle.normal(0)*circle.g, circle.normal(1)*circle.g, circle.normal(2)*circle.g, circle.form);
		fprintf(file, "intersector %d radius %f vector %f %f %f\n", i, circle.sphereRadius, circle.vector(0), circle.vector(1), circle.vector(2));
	}
	
	
	
	for(int i=0;i<sasas.size();++i){
		sasa = sasas[i];
		fprintf(file, "sasa %d size %d\n", i, sasa.sasa.size());
		k=0;
		for(int j=0; j<sasa.sasa.size(); ++j){
			fprintf(file, "intersectionpoint %d vector %f %f %f\n", j, sasa.sasa[j].vector(0), sasa.sasa[j].vector(1), sasa.sasa[j].vector(2));
		}
		
	}
	
	fclose(file);
	
		
	

}
