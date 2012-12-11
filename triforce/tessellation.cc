#include "tessellation.h"


Tessellation::Tessellation(Molecule &m){
	molecule = m;
}


vector<double*> &Tessellation::fetchForcePointers(){
	return forces;
}

vector<vector<double*> > &Tessellation::fetchAreaPointers(){
	return areas;
}



void Tessellation::build(bool split){
	CircularRegionsPerAtom circlesPerAtom;
	sasasForMolecule.clear();
	
	//molecule.update();
	atoms = molecule.fetchCoordinates();
	radii = molecule.fetchRadii();
	forces = molecule.fetchForcePointers();
	areas = molecule.fetchAreaPointers();
	
	for(int i=0; i<radii->size();i++){
		double radius = radii->at(i);
		printf("RADIUS[%d]: %f\n",i,radius);
	}
	
	
	
	//atoms.size()
	//iterate over all atoms and build the tessellation for each of them
	//for(int i=0; i<atoms.size(); ++i){
		int i=0;
		buildGaussBonnetPath(i, atoms, radii, sasasForMolecule, split);
		
		
		
		for(int i=0; i<sasasForMolecule[0].sasas.size(); ++i){
			printf("SASA %d\n",i);
			outputGaussBonnetPath(sasasForMolecule[0].sasas[i]);
		}


		
		//outputCircularRegions(*circles[circles.size()-1]);

	//}
	
	
	
}


CircularRegionsPerAtom Tessellation::coverHemisphere(Vector tessellationOrigin, double radius, CircularRegionsPerAtom circles, CircularInterfaceForm form){
	CircularRegion C;
	Vector v(3);
	
	v = -tessellationOrigin;
	
	C.openingAngle = M_PI * 0.5;
	C.a = radius;
	C.g=0;
	C.g_normalised=0;
	C.vector = v;
	C.normal = v;
	C.form = form;
	C.id = circles.size()+1;
	
	circles.push_back(C);
	
	return circles;
	
}



void Tessellation::buildGaussBonnetPath(int i, vector<Vector> &atoms, vector<double> &radii, SASAsForMolecule &sasas, bool split){
	int k,j;
	CircularRegionsPerAtom circles;
	CircularRegionsPerAtom circlesFrontHemisphere;
	CircularRegionsPerAtom circlesBackHemisphere;
	
	Vector frontTessellationOrigin(3);
	Vector backTessellationOrigin(3);
	Vector origin;
	double radius;
	
	
	origin = atoms[i];
	radius = radii[i];
	
	
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
	


	makeCircularRegions(i,origin, radius, atoms, radii, circles);
	for(j=0;j<circles.size();j++){
		determineProjection(origin, radius, circles[j]);
	}
	
	//there is room for optimisation here...
	
	if(split){
	
		circlesFrontHemisphere = coverHemisphere(frontTessellationOrigin, radius, circles, SPLITTER);
		circlesBackHemisphere = coverHemisphere(backTessellationOrigin, radius, circles, SPLITTER);
		
		
		filterCircularRegions(frontTessellationOrigin,radius, circlesFrontHemisphere);
		filterCircularRegions(backTessellationOrigin,radius, circlesBackHemisphere);
		reindexCircularRegions(circlesFrontHemisphere);
		reindexCircularRegions(circlesBackHemisphere);
		
		determineCircularIntersections(circlesFrontHemisphere);
		determineCircularIntersections(circlesBackHemisphere);
		
		printf("CIRC: %d, FRONT: %d, BACK: %d\n",circles.size(), circlesFrontHemisphere.size(), circlesBackHemisphere.size());
		
		buildIntersectionGraph(radius, frontTessellationOrigin, circlesFrontHemisphere, *newSasas, FRONTHEMISPHERE, string("gbonnet0.csv"));
		buildIntersectionGraph(radius, backTessellationOrigin, circlesBackHemisphere, *newSasas, BACKHEMISPHERE, string("gbonnet1.csv"));
	}
	else{
		
		
		filterCircularRegions(frontTessellationOrigin, radius, circles);
		reindexCircularRegions(circles);
		
		determineCircularIntersections(circles);
		
		buildIntersectionGraph(radius, frontTessellationOrigin, circles, *newSasas, FRONTHEMISPHERE, string("gbonnet0.csv"));
	}
	
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
	//printf("V: %f, EPS: %f, V+T: %f, EPS-V+T: %f\n",v,eps,v+THRESHOLD_NUMERICAL, eps-(v+THRESHOLD_NUMERICAL));
	if(eps-(v+THRESHOLD_NUMERICAL) <= 0) return true;
	else return false;
}

bool Tessellation::isWithinNumericalLimits(double x, double l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) return true;
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
		circle.normal=-circle.normal;
	}
	else circle.form=CONVEX;

	circle.openingAngle = acos(g / r_i);
	circle.g = g;
	circle.g_normalised = g/radius;
	circle.a = a;

	printf("DET: [%f %f %f] %f %f %f %f\n",circle.normal(0),circle.normal(1),circle.normal(2), r_i, r_k, g, circle.openingAngle);


	printf("%f %f %f %f %f\n",d_k, r_i, r_k, g, circle.openingAngle);
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



void Tessellation::makeCircularRegions(int i,Vector &origin, double radius, vector<vec> &atoms, vector<double> &radii, vector<CircularRegion> &circles){
	int j;
	Vector res,normal,v;
	CircularRegion circle;
	double r_i = radius;
	double r_k;
	double lenv;

	for(j=0;j<atoms.size();j++){
		if(i != j){
			v=atoms[j] - origin;
			r_k = radii[j];
			
			
			lenv = norm(v,2);
			
			printf("RADIUS: %f LEN: %f\n",r_k,lenv);
			

			//reject, if no intersection
			if(lenv < r_i + r_k && lenv+r_k > r_i && lenv+r_i > r_k){
				normal = v / lenv;
				circle.id = circles.size()+1;
				circle.vector = v;
				circle.normal = normal;
				circle.sphereRadius = r_k;
				circle.intersect = false;
				circles.push_back(circle);
				
				printf("CIRCLE[%d]: (%f, %f, %f) %f\n",circle.id, circle.vector(0), circle.vector(1), circle.vector(2), circle.sphereRadius);
			}
		}
	}
}



void Tessellation::depleteCircularRegions(Vector tessellationOrigin, double radius, vector<CircularRegion> &circles){
	circles.clear();
	Vector t;
	t = -tessellationOrigin;
	circles = coverHemisphere(t, radius, circles, CONCAVE);
}

int Tessellation::filterCircularRegions(Vector tessellationOrigin, double radius, vector<CircularRegion> &circles){
	vector<CircularRegion>::iterator it;
	double angle;
	int i;
	bool erased;
	Vector n0,n1;
	it = circles.begin();
	printf("----------------\n");
	while(it != circles.end()){
		erased=false;
		printf("I HAVE :%d\n",it->id);
	
		if(it->form==CONVEX) n0 = it->normal;
		else n0 = -it->normal;
		
		for(i=0;i<circles.size();i++){
			if(it->id != circles[i].id){
				printf("I LOOK AT :%d\n",circles[i].id);
				
				if(circles[i].form == CONVEX) n1 = circles[i].normal;
				else n1 = -circles[i].normal;
				
				
				
				angle = getAngleBetweenNormals(n0,n1);
				printf("ANGLE: %f %f %f\n",angle,it->openingAngle, circles[i].openingAngle);
						
						
				printf("V: (%f, %f, %f)   (%f %f %f)\n",n0(0),n0(1),n0(2),n1(0),n1(1),n1(2));
				
				
				//angle = getAngleBetweenNormals(it->normal, circles[i].normal);
				
				if(it->form==CONVEX){
					//convex circle IT is inside of convex circle i
					if(circles[i].form == CONVEX){
						
						if(it->openingAngle + angle < circles[i].openingAngle){
							printf("CASE 0\n");
							it = circles.erase(it);
							erased=true;
							break;
						}
							printf("CASE -0\n");
						
					}
					//convex circle is outside of concave circle i
					else{
						if(angle-it->openingAngle >  circles[i].openingAngle){
							printf("CASE 1\n");
							it = circles.erase(it);
							erased=true;
							break;
						}
							printf("CASE -1\n");
						
					}
					
				}
				else{
					//concave circle IT has a free area. This area is covered by convex circle i
					if(circles[i].form == CONVEX){
						if(it->openingAngle + angle < circles[i].openingAngle){
							printf("CASE 2\n");
							depleteCircularRegions(tessellationOrigin, radius, circles);
							return -1;
						}
						printf("CASE -2\n");
						
					}
					else{
						//concave circle IT is completely inside the occlusion of concave circle i
						if(it->openingAngle > circles[i].openingAngle + angle){
							printf("CASE 3\n");
							it = circles.erase(it);
							erased=true;
							break;
						}
						//concave circle IT has a free area. This area is covered by concave circle i
						else if(angle > it->openingAngle + circles[i].openingAngle){
							printf("CASE 4\n");
							depleteCircularRegions(tessellationOrigin, radius, circles);
							return -1;
							
						}
						printf("CASE -3\n");
						
					}

				}	
				
						
					
			}
		}
		if(!erased) ++it;
	}
	
	
	//if just the splitter is left, the hemisphere is completely uncovered.
	if(circles.size()==1 && circles[0].form==SPLITTER)
		circles.clear();
	
	return 0;

}




void Tessellation::outputGaussBonnetPath(SASA &points){
	SASANodeList::iterator it;
	int i;
	
	printf("ORIGIN: %f %f %f\n",points.tessellationOrigin(0),points.tessellationOrigin(1),points.tessellationOrigin(2));
	for(it=points.sasa.begin(), i=0; it!=points.sasa.end(); ++it, ++i)
		printf("GBPATH[%d] %d - %d\n", i, it->id0, it->id1);
}



void Tessellation::reindexCircularRegions(CircularRegionsPerAtom &circles){
	int i,j;
	for(i=0;i<circles.size();i++){
		circles[i].id = i+1;
	}
}




void Tessellation::insertFakeIntersectionPoints(CircularRegion &I, IntersectionGraph &intersectionGraph, Vector &tessellationOrigin){
	Vector v,v2,o,o2;
	int k;
	IntersectionNode a,b;
	IntersectionAddress addressA,addressB;
	Vector x0,x1;
	Interfaces f;
	
	
	
	
	//this will create two points on the border of the circular region on opposite sides.
	
	

	if(I.g==0){
		v = randu<vec>(3);
		v2 = I.normal;
		o = cross(v2,v);
		o = I.a*o / norm(o,2);
		x0 = o;
		
		x1 = (o * -1);
	}
	else{
		v = randu<vec>(3);
		v2 = I.normal * I.g;
		o = cross(v2,v);
		o = I.a*o / norm(o,2);
		x0 = v2+o;
		
		x1 = v2 + (o * -1);
	}
	
	f = angularInterfaces(x0,x1, tessellationOrigin, I);
	
	printf("F0: %f %f\n",f.in,f.out);
	
	a.id0 = addressA.id0 = I.id;
	a.id1 = addressA.id1 = -I.id;
	a.pointsTo0 = a.prev0 = -I.id;
	a.pointsTo1 = a.prev1 = I.id;
	a.vector = x0;
	a.angle1=f.in;
	a.angle0=f.out;
	a.visited=false;

	b.id0 = addressB.id0 = -I.id;
	b.id1 = addressB.id1 = I.id;
	b.pointsTo0 = a.prev0 = I.id;
	b.pointsTo1 = a.prev1 = -I.id;
	b.vector = x1;
	b.angle1=-f.in;
	b.angle0=-f.out;
	b.visited=false;
	
	intersectionGraph[addressA]=a;
	intersectionGraph[addressB]=b;
	
			
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
	Vector t;
	
	
	//measurement points will be written into p0 and p1
	if(I.form != CONVEX){
		t = -tessellationOrigin;
		measurementPoints(p0,p1,t,I);
	}
	else measurementPoints(p0,p1,tessellationOrigin,I);
	
	
	v = I.normal * I.g;
	
	
	res.out = angularInterface(x0,v,p0,p1);
	res.vectorOut = x0;
	res.in = angularInterface(x1,v,p0,p1);
	res.vectorIn = x1;
	
	if(I.form != CONVEX){
		res.out = -res.out;
		res.in = -res.in;
	}
	
	return res;
	
}




double Tessellation::rotationalAngle(Vector &tessellationOrigin, CircularRegion &I, CircularRegion &J){
	Vector ex(3);
	Vector ez(3);
	Vector n0(3);
	Vector n1(3);
	Vector n2(3);
	double sn;
	double a;
	double rho;
	
	ez(0)=0;
	ez(1)=0;
	ez(2)=1;
	double d;

	
	d = dot(tessellationOrigin,I.normal);
	
	if(isWithinNumericalLimits(d,1.0)){
		n0 = Vector(3);
		n0(0) = 0;
		n0(1) = 0;
		n0(2) = 1;
		printf("POSITIVE RANGE +1 %f\n",d);
	}
	else if(isWithinNumericalLimits(d,-1.0)){
		n0 = Vector(3);
		n0(0) = 0;
		n0(1) = 0;
		n0(2) = -1;
		printf("POSITIVE RANGE -1 %f\n",d);
	}
	else n0 = cross(tessellationOrigin,I.normal);
	
	n2 = cross(I.normal,n0);
	
	n1 = cross(I.normal, J.normal);
	
	s0 = sgn(norm_dot(n1,n2));
	a = asin(norm_dot(n0,n1));
	s1 = sgn(a);
	
	//if(sn<0 && a>=0) a = M_PI-a;
	//else if(sn<0 && a<0) a = -M_PI-a;
	
	a = -s0*(s1*(1-s0)*pi/2 - a);


	
	
	return a;
	
}



Interfaces Tessellation::retrieveInterfaces(Vector &tessellationOrigin, CircularRegion &I, CircularRegion J, double dij, double radius){
	double gij;
	double eta;
	double mu;
	double entryPoint, exitPoint;
	IntersectionPair ip;
	Vector x0(3), x1(3);
	Interfaces intf;
	Interfaces intf1;
	double rotationalPartIJ;
	double rotationalPartJI;
	double lambda_j, lambda_k, rho;
	
/*	
	if(PHIInterpolation){
		baseAngleIJ = dataPHI->interpolate(I.g_normalised,J.g_normalised,norm_dot(I.normal,J.normal));
		intf.out = baseAngleIJ + rotationalPartIJ;
		if(intf.out>=M_PI) intf.out = (baseAngleIJ + rotationalPartIJ) - 2*M_PI;

		
		printf("gi %f gj %f cosrho %f\n",I.g_normalised,J.g_normalised,norm_dot(I.normal,J.normal));
		printf("BASE ANGLE: %f %f ROTATIONAL PART: %f %f\n",baseAngleIJ,baseAngleJI,rotationalPartIJ,rotationalPartJI);
		
		//return intf;
		
		
		
		intf.out = baseAngleIJ + rotationalPartIJ;
		if(intf.out>=M_PI) intf.out = (baseAngleIJ + rotationalPartIJ) - 2*M_PI;

		intf.vectorOut = Vector(3);
		intf.vectorOut(0)=0;
		intf.vectorOut(1)=0;
		intf.vectorOut(2)=0;
		
		//baseAngleJI = dataPHI->interpolate(J.g_normalised,I.g_normalised,norm_dot(J.normal,I.normal));
		//rotationalPartJI = rotationalAngle(J,I);
		intf.in = (-M_PI-baseAngleIJ) + rotationalPartIJ;
		if(intf.in<=-M_PI) intf.in = (M_PI-baseAngleIJ) + rotationalPartIJ;
		
		
	}
*/	
	
	lambda_j = acos(I.g_normalised);
	lambda_k = acos(J.g_normalised);
	rho = acos(norm_dot(I.normal,J.normal));
	
	eta = M_PI/2 - acos(-(1/tan(lambda_j))*(1/tan(rho))+cos(lambda_k)*(1/sin(lambda_j))*(1/sin(rho)));
	
	printf("VECTORS (%f, %f, %f) (%f, %f, %f) (%f, %f, %f)\n",I.normal(0),I.normal(1),I.normal(2),J.normal(0),J.normal(1),J.normal(2), tessellationOrigin(0), tessellationOrigin(1), tessellationOrigin(2));
	
	rotationalPartIJ = rotationalAngle(tessellationOrigin,I,J);
	
	if(I.form != CONVEX){
		if(rotationalPartIJ>=0)
			rotationalPartIJ = -M_PI + rotationalPartIJ;
	}

	
	intf.vectorOut = Vector(3);
	intf.vectorOut(0)=0;
	intf.vectorOut(1)=0;
	intf.vectorOut(2)=0;
	intf.out = eta + rotationalPartIJ;
	if(intf.out>=M_PI) intf.out = (eta + rotationalPartIJ) - 2*M_PI;
	
	//baseAngleJI = dataPHI->interpolate(J.g_normalised,I.g_normalised,norm_dot(J.normal,I.normal));
	//rotationalPartJI = rotationalAngle(J,I);
	
	intf.in = (-M_PI-eta) + rotationalPartIJ;
	if(intf.in<=-M_PI) intf.in = (M_PI-eta) + rotationalPartIJ;
	intf.vectorIn = Vector(3);
	intf.vectorIn(0)=0;
	intf.vectorIn(1)=0;
	intf.vectorIn(2)=0;
	
	printf("gi %f gj %f cosrho %f\n",I.g_normalised,J.g_normalised,norm_dot(I.normal,J.normal));
	printf("BASE ANGLE: %f ROTATIONAL PART: %f\n",eta,rotationalPartIJ);
	
	
	
	
	
	
	ip = determineIntersectionPoints(radius, I, J);
	x0 = ip.k_j;
	x1 = ip.j_k;
	
	
	intf1 = angularInterfaces(x0,x1, tessellationOrigin, I);


	
	printf("TESTING: error:(%f %f) out: (%f %f) in: (%f %f) \n",abs(intf.out-intf1.out),abs(intf.in-intf1.in),intf.out,intf1.out,intf.in,intf1.in);
	double threshold=0.1;
	if(min(abs(intf.out-intf1.out),2*M_PI - abs(intf.out-intf1.out)) >= threshold || abs(intf.in-intf1.in) >= threshold){
		printf("ABORTING\n");
		exit(-2);
	}
	
	
	
	return intf;
	
	
	
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
	a.pointsTo0 = 0;
	a.pointsTo1 = 0;
}
void Tessellation::connectIntersectionPoints(IntersectionNode &a, IntersectionNode &b, IntersectionGraph &intersectionGraph){
	IntersectionAddress c;
	IntersectionAddress d;

	if(a.pointsTo0 != -1 && a.pointsTo1 != -1){
		c.id0=a.pointsTo0;
		c.id1=a.pointsTo1;
		
		intersectionGraph[c].prev0 = 0;
		intersectionGraph[c].prev1 = 0;
	}

	
	if(b.prev0 != -1 && b.prev1 != -1){
		d.id0=b.prev0;
		d.id1=b.prev1;
		
		intersectionGraph[d].pointsTo0 = 0;
		intersectionGraph[d].pointsTo1 = 0;
	}
	
	
	a.pointsTo0 = b.id0;
	a.pointsTo1 = b.id1;

	b.prev0 = a.id0;
	b.prev1 = a.id1;
	
}

void Tessellation::deleteIntersectionPoint(IntersectionBranches::iterator &it,IntersectionGraph &intersectionGraph){
	//delete point from intersectionGraph
	IntersectionAddress address;
	IntersectionBranches::iterator it_mirror;
	
	
	if(it->second.node->prev0 != 0 && it->second.node->prev1 != 0){
		address.id0 = it->second.node->prev0;
		address.id1 = it->second.node->prev1;
		
		intersectionGraph[address].pointsTo0 = 0;
		intersectionGraph[address].pointsTo1 = 0;
	}
	
	if(it->second.node->pointsTo0 != 0 && it->second.node->pointsTo1 != 0){
		address.id0 = it->second.node->pointsTo0;
		address.id1 = it->second.node->pointsTo1;
		
		intersectionGraph[address].prev0 = 0;
		intersectionGraph[address].prev1 = 0;
	}
	
	address.id0 = it->second.node->id0;
	address.id1 = it->second.node->id1;
	
	intersectionGraph.erase(address);
	
	//first destroy the mirror
	it_mirror = it->second.it;
	it_mirror->second.body->erase(it_mirror);
	
	//then destroy the point
	it->second.body->erase(it);
	
	//from now on "it" is invalid...
	
}



void Tessellation::createIntersectionNode(IntersectionAddress &address, IntersectionGraph &intersectionGraph){
	IntersectionNode node;
	
	node.id0 = address.id0;
	node.id1 = address.id1;
	node.pointsTo0=0;
	node.pointsTo1=0;
	node.prev0=0;
	node.prev1=0;
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
	
	x.second.visited = -1;
	
	x.first = interfacesJ.in;
	x.second.node=&intersectionGraph[address];
	x.second.direction = IN;
	x.second.it = I.intersectionBranches.end();
	x.second.body = &J.intersectionBranches;
	x.second.id = address.id1;
	it0 = J.intersectionBranches.insert(x);
	
	//printf("address %d-%d side %d IN\n",address.id0, address.id1, address.id1);
	
	x.first = interfacesI.out;
	x.second.node=&intersectionGraph[address];
	x.second.direction = OUT;
	x.second.it = it0;
	x.second.body = &I.intersectionBranches;
	x.second.id = address.id0;
	it1 = I.intersectionBranches.insert(x);
	it0->second.it=it1;

	//printf("address %d-%d side OUT\n",address.id0, address.id1, address.id0);
	
	
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


void Tessellation::buildIntersectionGraph(double radius, Vector &tessellationOrigin, CircularRegionsPerAtom &circles, SASAs &sasas, Hemisphere hemisphere, string filename){
	map<int, bool> processed;
	map<int, bool>::iterator it_p;
	IntersectionGraph intersectionGraph;
	IntersectionGraph::iterator it_g, it_x;
	int cid0;
	
	SASA potentialSasa;
	SASANode sasaNode;
	
	int i,j;
	CircularRegion *I, *J;
	double tau0, tau1;
	map<int,CircularIntersection>::iterator it_j;
	
	bool empty;
	IntersectionAddress start, x, t;
	bool valid;
	map<IntersectionBranches::iterator,bool,IteratorComparator> eraseList;
	map<IntersectionBranches::iterator,bool,IteratorComparator>::iterator it_e;
	
	IntersectionPair ip;
	
	
	Interfaces interfacesJ, interfacesI;
	IntersectionBranches::iterator it_main;
	IntersectionBranches::iterator it0, it1;
	IntersectionBranches::iterator it, it_mirror, it_next, it_prev, it_mirror_next, it_mirror_prev, it_mirror_next_ignore, it_mirror_prev_ignore;
	IntersectionAddress addressIJ, addressJI;
	
	
	
	
	//iterate through all circles and add them to the intersectiongraph, one by one
	for(i=0; i < circles.size(); ++i){
		
		//printf("adding %d\n",i);
		I = &circles[i];
		processed[i]=true;
		
		
		
		
		
		//we have to sort this circle into the intersectiongraph
		//go through all intersecting circles, check whether they have been processed, and if so, calculate the intersections
		
		if(I->circularIntersections.size()==0 && I->form!=SPLITTER){
			insertFakeIntersectionPoints(*I,intersectionGraph,tessellationOrigin);			
		}
		
		for(it_j=I->circularIntersections.begin(); it_j != I->circularIntersections.end(); ++it_j){
			j=(*it_j).first;
			J = &circles[j];
			
			
			
			
			
			it_p = processed.find(j);
			if(it_p != processed.end()){
				//printf("processing intersection with %d\n",j);
				
				
				addressIJ.id0 = I->id;
				addressIJ.id1 = J->id;
				addressJI.id0 = J->id;
				addressJI.id1 = I->id;
				
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
			it_main->second.visited=-1;
		}
		
		eraseList.clear();
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
					//printf("EMPTY\n");
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
				
				if(empty){
					connectIntersectionPoints(*it_mirror->second.node, *it_next->second.node, intersectionGraph);
				}
				else{
					//external
					if(it_mirror_prev_ignore->second.direction == IN){
						if(it->second.direction == IN){
							printf("EXTERNAL IN\n");
							
							if(it_mirror_prev->second.direction==IN && it_next->second.direction==OUT){
								connectIntersectionPoints(*it_mirror_prev->second.node, *it_mirror->second.node, intersectionGraph);
								connectIntersectionPoints(*it_mirror->second.node, *it_next->second.node, intersectionGraph);
								
								if(!empty && I->form!=CONVEX && it_mirror_next->second.it->second.id != it->second.id){
									//deleteIntersectionPoint(it_mirror_next, intersectionGraph);
									eraseList[it_mirror_next]=true;
									printf("DELETE %d-%d (NEXT)\n",it_mirror_next->second.node->id0,it_mirror_next->second.node->id1);

								}
								
							}
							else printf("REJECT\n");
							
						}
						else{
							printf("EXTERNAL OUT\n");
							if(it_prev->second.direction==IN && it_mirror_next->second.direction==OUT){
								connectIntersectionPoints(*it_mirror->second.node, *it_mirror_next->second.node, intersectionGraph);
								connectIntersectionPoints(*it_prev->second.node, *it->second.node, intersectionGraph);
								
								if(!empty && I->form==CONVEX && it_mirror_prev->second.it->second.id != it->second.id){
									//deleteIntersectionPoint(it_mirror_prev, intersectionGraph);
									eraseList[it_mirror_prev]=true;
									printf("DELETE %d-%d (PREV)\n",it_mirror_prev->second.node->id0,it_mirror_prev->second.node->id1);
									
								}
							
							}
							else printf("REJECT\n");
							
						}
					}
					//internal
					else{
							
							if(it->second.direction == IN){
								if(!empty && it_mirror_next->second.it->second.id != it->second.id){
									//deleteIntersectionPoint(it_mirror_next, intersectionGraph);
									eraseList[it_mirror_next]=true;
									
									printf("DELETE CONV IN  %d-%d (NEXT)\n",it_mirror_next->second.node->id0,it_mirror_next->second.node->id1);
								}
							}
							else{
								if(!empty && it_mirror_prev->second.it->second.id != it->second.id){
									//deleteIntersectionPoint(it_mirror_prev, intersectionGraph);
									eraseList[it_mirror_prev]=true;
									
									printf("DELETE CONV OUT %d-%d (PREV)\n",it_mirror_prev->second.node->id0,it_mirror_prev->second.node->id1);
									
								}
							}
							
						
						
						eraseList[it]=true;
						//deleteIntersectionPoint(it, intersectionGraph);
						printf("DELETE POINT\n");
					}
				}
				
				
				printf("-------------------------------------------------\n");
				printIntersectionGraph(intersectionGraph);
				printf("-------------------------------------------------\n");
				
			}
			
		}
		
		printf("REMOVING IPS\n");
		for(it_e=eraseList.begin(); it_e != eraseList.end(); ++it_e){
			it = it_e->first;
			
			printf("REMOVING IP: %d\n",&*it);
			printf("--: %d-%d\n",it->second.node->id0,it->second.node->id1);
			deleteIntersectionPoint(it,intersectionGraph);
		}
		printf("REMOVING IPS - DONE\n");
		
	}
	printf("GENERATING SASAS\n");
	
	int s = 0;
	
	for(it_g = intersectionGraph.begin(); it_g != intersectionGraph.end(); ++it_g){
		if(!it_g->second.visited){
			start.id0 = it_g->second.id0;
			start.id1 = it_g->second.id1;
			
			x.id0 = start.id0;
			x.id1 = start.id1;
			
			printf("STARTING WITH: %d - %d\n",x.id0,x.id1);
			
			
			valid=true;
			potentialSasa.sasa.clear();
			do{
				
				if(intersectionGraph[x].visited || (intersectionGraph[x].pointsTo0==0 && intersectionGraph[x].pointsTo1==0)){
					printf("BREAK\n");
					valid = false;
					break;
					
				}
				
				intersectionGraph[x].visited=true;
				//if(s > 20) exit(-1);
				++s;
				
				//all id's start at 1, so we have to decrease them by 1. 
				cid0=abs(x.id0)-1;
				
				sasaNode.id0 = x.id0;
				sasaNode.id1 = x.id1;
				sasaNode.angle0 = intersectionGraph[x].angle0;
				sasaNode.angle1 = intersectionGraph[x].angle1;
				sasaNode.lambda = circles[cid0].openingAngle;
				sasaNode.vector = intersectionGraph[x].vector;
				sasaNode.normalForCircularRegion = circles[cid0].normal;
				sasaNode.form = circles[cid0].form;
				
				int cid1=abs(x.id1)-1;
				double ank = getAngle(circles[cid0].normal, circles[cid1].normal);
				printf("VERIFICATION.. g0: %f g1: %f cosrho: %f PHI0: %f PHI1: %f rotational part: %f n0(%f, %f, %f) n1(%f, %f, %f)\n",circles[cid0].g/radius,circles[cid1].g/radius, cos(ank), sasaNode.angle0,  sasaNode.angle1, rotationalAngle(tessellationOrigin, circles[cid0],circles[cid1]),circles[cid0].normal(0),circles[cid0].normal(1),circles[cid0].normal(2),circles[cid1].normal(0),circles[cid1].normal(1),circles[cid1].normal(2)); 
				
				
				
				
				potentialSasa.sasa.push_back(sasaNode);
				
				t.id0 = intersectionGraph[x].pointsTo0;
				t.id1 = intersectionGraph[x].pointsTo1;
				
				x = t;
				
				
				
					
				printf("NEXT: %d - %d\n",x.id0,x.id1);
				
				
			}
			while(!(intersectionGraph[x].id0 == start.id0 && intersectionGraph[x].id1 == start.id1));
			
			if(valid){
				printf("ADDING SASA\n");
				potentialSasa.tessellationOrigin = tessellationOrigin;
				potentialSasa.hemisphere = hemisphere;
				sasas.push_back(potentialSasa);
			}
		}
	}
	


	outputGaussBonnetData(filename, radius, circles, sasas, intersectionGraph);

	
	
}












void Tessellation::outputGaussBonnetData(string filename, double radius, CircularRegionsPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph){
	FILE* file;
	
	file = fopen (filename.c_str(),"w");
	
	SASA sasa;
	SASANodeList::iterator it;
	CircularRegion circle;
	double area=0;
	int k;
	
	fprintf(file, "atom radius %f\n", radius);
	
	
	for(int i=0; i< circles.size(); ++i){
		circle = circles[i];
		fprintf(file, "circularregion %d radius %f vector %f %f %f form %d\n", circle.id, circle.a, circle.normal(0)*circle.g, circle.normal(1)*circle.g, circle.normal(2)*circle.g, circle.form);
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
