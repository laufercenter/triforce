#include "tessellation.h"


Tessellation::Tessellation(Molecule &m){
	molecule = m;
}



void Tessellation::build(bool split){
	CircularInterfacesPerAtom circlesPerAtom;
	sasasForMolecule.clear();
	
	//molecule.update();
	atoms = molecule.fetchCoordinates();
	radii = molecule.fetchRadii();
	
	for(int i=0; i<radii.size();i++){
		double radius = radii[i];
		printf("RADIUS[%d]: %f\n",i,radius);
	}
	
	
	
	//atoms.size()
	//iterate over all atoms and build the tessellation for each of them
	for(int i=0; i<atoms.size(); ++i){
		//int i=16;
		buildGaussBonnetPath(i, atoms, radii, sasasForMolecule, split);
		
		
		
		for(int i=0; i<sasasForMolecule[0].sasas.size(); ++i){
			printf("SASA %d\n",i);
			outputGaussBonnetPath(sasasForMolecule[0].sasas[i]);
		}


		
		//outputCircularInterfaces(*circles[circles.size()-1]);

	}
	
	
	
}



CircularInterfacesPerAtom Tessellation::coverHemisphere(Vector tessellationOrigin, double radius, CircularInterfacesPerAtom circles, CircularInterfaceForm form){
	CircularInterface C;
	Vector v(3);
	Matrix NullMatrix(3,3);
	NullMatrix.zeros();
	
	v = tessellationOrigin;
	
	C.lambda.rotation = M_PI * 0.5;
	C.g=0;
	C.g_normalised=0;
	C.vector = v;
	C.normal = v;
	C.form = form;
	C.id = circles.size()+1;
	C.a=radius;
	C.valid = true;
	
	C.dmu_dx = NullMatrix;
	C.lambda.drotation_dxi = Vector(3).zeros();
	C.lambda.drotation_dxj = Vector(3).zeros();
	C.lambda.drotation_dxl = Vector(3).zeros();
	
	
	circles.push_back(C);
	
	return circles;
	
}

void Tessellation::buildGaussBonnetPath(int i, vector<Vector> &atoms, vector<double> &radii, SASAsForMolecule &sasas, bool split){
	int k,j;
	CircularInterfacesPerAtom circles;
	CircularInterfacesPerAtom circlesFrontHemisphere;
	CircularInterfacesPerAtom circlesBackHemisphere;
	
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
	


	makeCircularInterfaces(i,origin, radius, atoms, radii, circles);
	for(j=0;j<circles.size();j++){
		determineProjection(origin, radius, circles[j]);
	}
	
	//there is room for optimisation here...
	
	if(split){
	
		circlesFrontHemisphere = coverHemisphere(frontTessellationOrigin, radius, circles, SPLITTER);
		circlesBackHemisphere = coverHemisphere(backTessellationOrigin, radius, circles, SPLITTER);
		
		
		filterCircularInterfaces(frontTessellationOrigin,radius, circlesFrontHemisphere);
		filterCircularInterfaces(backTessellationOrigin,radius, circlesBackHemisphere);
		
		determinePsiRotations(frontTessellationOrigin, circlesFrontHemisphere);
		determinePsiRotations(backTessellationOrigin, circlesBackHemisphere);
		
		reindexCircularInterfaces(circlesFrontHemisphere);
		reindexCircularInterfaces(circlesBackHemisphere);
		
		determineCircularIntersections(circlesFrontHemisphere);
		determineCircularIntersections(circlesBackHemisphere);
		
		printf("CIRC: %d, FRONT: %d, BACK: %d\n",circles.size(), circlesFrontHemisphere.size(), circlesBackHemisphere.size());
		
		buildIntersectionGraph(radius, frontTessellationOrigin, circlesFrontHemisphere, *newSasas, FRONTHEMISPHERE, string("gbonnet0.csv"));
		buildIntersectionGraph(radius, backTessellationOrigin, circlesBackHemisphere, *newSasas, BACKHEMISPHERE, string("gbonnet1.csv"));
	}
	else{
		printf("not supported\n");
		exit(-1);
		
		filterCircularInterfaces(frontTessellationOrigin, radius, circles);
		reindexCircularInterfaces(circles);
		
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




void Tessellation::determineProjection(Vector &origin, double radius, CircularInterface &circle){
	double d_i;
	double r_l, r_i;
	double g, g_normalised;
	Vector mu;
	double a;
	Matrix Identity(3,3);
	Identity.eye();
	double dlambda_dg;
	double dg_dd;
	Matrix dmu_dx(3,3);
	
	
	
	
	d_i = norm(circle.vector,2);
	mu = circle.normal;
	r_l = radius;
	r_i = circle.sphereRadius;
	
	
	g = (d_i * d_i + r_l * r_l - r_i * r_i ) / (2 * d_i);
	
	if(g<0){
		circle.form=CONCAVE;
		g = abs(g);
		circle.normal=-circle.normal;
	}
	else circle.form=CONVEX;

	circle.g = g;
	g_normalised = g/radius;
	circle.g_normalised = g_normalised;
	circle.lambda.rotation = acos(g_normalised);
	circle.d=d_i;
	circle.a= sqrt(r_l * r_l - g * g);
	circle.valid=true;

	
	//derivatives
	dlambda_dg = -1/sqrt(1-g_normalised*g_normalised);
	dg_dd = -g_normalised/d_i + 1/r_l;
	
	circle.lambda.drotation_dxi = dlambda_dg * dg_dd * circle.normal;
	circle.lambda.drotation_dxj = Vector(3).zeros();
	circle.lambda.drotation_dxl = -circle.lambda.drotation_dxi;
	
	
	
	dmu_dx = (1.0/circle.d)*(Identity - kron(circle.normal,circle.normal.t()));
	circle.dmu_dx = dmu_dx;
	


	
	

	//printf("DET: [%f %f %f] %f %f %f %f\n",circle.normal(0),circle.normal(1),circle.normal(2), r_i, r_k, g, circle.lambda);
	//printf("%f %f %f %f %f\n",d_k, r_i, r_k, g, circle.lambda);
}


void Tessellation::determinePsiRotations(Vector &tessellationOrigin, CircularInterfacesPerAtom &circles){
	int i;
	for(i=0; i < circles.size(); ++i)
		circles[i].psi=calculatePsi(tessellationOrigin, circles[i]);
}


IntersectionPair Tessellation::determineIntersectionPoints(double radius, CircularInterface &K, CircularInterface &J){
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



void Tessellation::makeCircularInterfaces(int i,Vector &origin, double radius, vector<vec> &atoms, vector<double> &radii, vector<CircularInterface> &circles){
	int j;
	Vector res,normal,v;
	CircularInterface circle;
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



void Tessellation::depleteCircularInterfaces(Vector tessellationOrigin, double radius, vector<CircularInterface> &circles){
	circles.clear();
	Vector t;
	t = -tessellationOrigin;
	circles = coverHemisphere(t, radius, circles, CONCAVE);
}

int Tessellation::filterCircularInterfaces(Vector tessellationOrigin, double radius, vector<CircularInterface> &circles){
	vector<CircularInterface>::iterator it;
	double angle;
	int i;
	bool erased;
	Vector n0,n1;
	it = circles.begin();
	printf("----------------\n");
	while(it != circles.end()){
		erased=false;
	
		//if(it->form==CONVEX) n0 = it->normal;
		//else n0 = -it->normal;
		
		n0 = it->normal;
		
		for(i=0;i<circles.size();i++){
			if(it->id != circles[i].id){

				//in case the interface is not convex, we have to reverse the normal (since we already reversed it one time when we calculated the interface)
				//if(circles[i].form == CONVEX) n1 = circles[i].normal;
				//else n1 = -circles[i].normal;
				
				n1 = circles[i].normal;
				
				angle = getAngleBetweenNormals(n0,n1);
				
				printf("FILTER: %d %d: %f\n",it->id, circles[i].id, angle);
				
				if(it->form==CONVEX){
					//convex circle IT is inside of convex circle i
					if(circles[i].form == CONVEX){
						
						if(it->lambda.rotation + angle < circles[i].lambda.rotation){
							it = circles.erase(it);
							erased=true;
							break;
						}
					}
					//convex circle is outside of concave circle i
					else{
						if(angle-it->lambda.rotation >  circles[i].lambda.rotation){
							it = circles.erase(it);
							erased=true;
							break;
						}
					}
					
				}
				else{
					//concave circle IT has a free area. This area is covered by convex circle i
					if(circles[i].form == CONVEX){
						if(it->lambda.rotation + angle < circles[i].lambda.rotation){
							depleteCircularInterfaces(tessellationOrigin, radius, circles);
							return -1;
						}
						
					}
					else{
						//concave circle IT is completely inside the occlusion of concave circle i
						if(it->lambda.rotation > circles[i].lambda.rotation + angle){
							it = circles.erase(it);
							erased=true;
							break;
						}
						//concave circle IT has a free area. This area is covered by concave circle i
						else if(angle > it->lambda.rotation + circles[i].lambda.rotation){
							depleteCircularInterfaces(tessellationOrigin, radius, circles);
							return -1;
							
						}
						
					}

				}	
				
						
					
			}
		}
		if(!erased) ++it;
	}
	
	
	//if just the splitter is left, the hemisphere is completely free. In that case, the splitter is unnecessary
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



void Tessellation::reindexCircularInterfaces(CircularInterfacesPerAtom &circles){
	int i,j;
	for(i=0;i<circles.size();i++){
		circles[i].id = i+1;
	}
}





void Tessellation::insertArtificialIntersectionPoints(CircularInterface &I, IntersectionGraph &intersectionGraph, Vector &tessellationOrigin){
	Vector v,v2,o,o2;
	int k;
	IntersectionNode a,b;
	IntersectionAddress addressA,addressB;
	Vector x0,x1;
	Rotation f,f_reverse;
	
	
	
	
	//this will create two points on the border of the circular region on opposite sides.
	f.rotation = M_PI/2;
	f_reverse.rotation = -M_PI/2;
	
	a.id0 = addressA.id0 = I.id;
	a.id1 = addressA.id1 = -I.id;
	a.pointsTo0 = a.prev0 = -I.id;
	a.pointsTo1 = a.prev1 = I.id;
	a.vector = x0;
	a.rotation1=f;
	a.rotation0=f_reverse;
	a.visited=false;

	b.id0 = addressB.id0 = -I.id;
	b.id1 = addressB.id1 = I.id;
	b.pointsTo0 = a.prev0 = I.id;
	b.pointsTo1 = a.prev1 = -I.id;
	b.vector = x1;
	b.rotation1=f_reverse;
	b.rotation0=f;
	b.visited=false;
	
	intersectionGraph[addressA]=a;
	intersectionGraph[addressB]=b;
	
			
}


int Tessellation::sgn(double d){
	if(d>=0) return 1;
	else return -1;
}

void Tessellation::determineCircularIntersections(CircularInterfacesPerAtom &circles){
	int i,k,j;
	double angle;
	CircularIntersection c;
	
	
	//this needs to be optimised
	for(k=0;k<circles.size();k++)
		for(j=k+1;j<circles.size();j++){
			if(k!=j){
				//determine if there will be intersections
				angle = getAngleBetweenNormals(circles[k].normal,circles[j].normal);
				if(angle < circles[k].lambda.rotation + circles[j].lambda.rotation)
					if(angle + circles[k].lambda.rotation > circles[j].lambda.rotation && angle + circles[j].lambda.rotation > circles[k].lambda.rotation){
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

void Tessellation::measurementPoints(Vector &p0, Vector &p1, Vector &tessellationOrigin, CircularInterface &I){
	Vector v(3);
	Vector o(3);
	Vector vx(3);
	Vector s(3), s2(3);
	
	o = tessellationOrigin;
	v = I.normal * I.g;
	
	printf("O: %f %f %f\n",o(0),o(1),o(2));
	printf("V: %f %f %f\n",v(0),v(1),v(2));
	printf("INORMAL: %f %f %f\n",I.normal(0),I.normal(1),I.normal(2));
	printf("DOT: %f\n",dot(o,I.normal));
	
	
	if(isWithinNumericalLimits(dot(o,I.normal),1.0) ){
		p1 = Vector(3);
		p1(0) = 0;
		p1(1) = -1;
		p1(2) = 0;
		
		p0 = cross(p1,-o);
		printf("FIRST\n");
		
	}
	else 	if(isWithinNumericalLimits(dot(o,I.normal),-1.0) ){
		p1 = Vector(3);
		p1(0) = 0;
		p1(1) = 1;
		p1(2) = 0;
		
		p0 = cross(p1,-o);

		printf("SECOND\n");
		
	}
	else{
		

		
		s = cross(v,o);
		s = I.a*s / norm(s,2);
		p0 = v+s;
		
		s2 = -cross(v,s);
		s2 = I.a*s2 / norm(s2,2);
		p1 = v+s2;
		
		if(dot(tessellationOrigin,I.normal)<0){
			//p0 = v-s;
			//p1 = v-s2;
		}
		
		printf("THIRD\n");
		
	}
	
}

Interfaces Tessellation::angularInterfaces(Vector &x0, Vector &x1, Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J){
	Vector p0(3), p1(3);
	Interfaces res;
	Vector v(3);
	Vector t;
	int q;
	
	/*
	//measurement points will be written into p0 and p1
	if(I.form != CONVEX){
		t = -tessellationOrigin;
		measurementPoints(p0,p1,t,I);
	}
	else 
	*/
	measurementPoints(p0,p1,tessellationOrigin,I);
	v = I.normal * I.g;
	
	printf("I: %d J: %d\n",I.id,J.id);
	printf("V: %f %f %f\n",v(0),v(1),v(2));
	printf("P0: %f %f %f\n",p0(0),p0(1),p0(2));
	printf("P1: %f %f %f\n",p1(0),p1(1),p1(2));
	
	
	
	q = 1;
	if(I.form != CONVEX) q*=-1;
	if(J.form != CONVEX) q*=-1;
	
	
	if(q>0){
		res.out = angularInterface(x0,v,p0,p1);
		res.vectorOut = x0;
		res.in = angularInterface(x1,v,p0,p1);
		res.vectorIn = x1;
	}
	else{
		res.in = angularInterface(x0,v,p0,p1);
		res.vectorIn = x0;
		res.out = angularInterface(x1,v,p0,p1);
		res.vectorOut = x1;
	}
	
	
	
	return res;
	
}






PHIContainer Tessellation::retrieveInterfaces(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J, double dij, double radius){
	double gij;
	double eta;
	double mu;
	double entryPoint, exitPoint;
	IntersectionPair ip;
	Vector x0(3), x1(3);
	Interfaces intf;
	Interfaces intf1;
	double baseAngleIJ;
	double rotationalPartIJ;
	double baseAngleJI;
	double rotationalPartJI;
	PHIContainer p;
	Vector Null(3);
	Null.zeros();
	
	
		ip = determineIntersectionPoints(radius, I, J);
		x0 = ip.k_j;
		x1 = ip.j_k;
		
		
		intf1 = angularInterfaces(x0,x1, tessellationOrigin, I, J);
		
		

	p.in.rotation = intf1.in;
	p.in.drotation_dxi = Null;
	p.in.drotation_dxj = Null;
	p.in.drotation_dxl = Null;
	p.in.vector=intf1.vectorIn;
	p.out.rotation = intf1.out;
	p.out.drotation_dxi = Null;
	p.out.drotation_dxj = Null;
	p.out.drotation_dxl = Null;
	p.out.vector=intf1.vectorOut;
	
	printf("INTERFACES: %f %f\n",p.in.rotation, p.out.rotation);
	printf("X0: %f %f %f\n",x0(0),x0(1),x0(2));
	printf("X1: %f %f %f\n",x1(0),x1(1),x1(2));

		
	return p;
	
	
}




double Tessellation::rotationalAngle(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J){
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
	double s0,s1;

	
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
		printf("NEGATIVE RANGE -1 %f\n",d);
	}
	else n0 = cross(tessellationOrigin,I.normal);
	
	n2 = cross(I.normal,n0);
	
	n1 = cross(I.normal, J.normal);
	
	s0 = sgn(norm_dot(n1,n2));
	a = asin(norm_dot(n0,n1));
	s1 = sgn(a);
	
	//if(sn<0 && a>=0) a = M_PI-a;
	//else if(sn<0 && a<0) a = -M_PI-a;
	
	a = -s0*(s1*(1-s0)*M_PI/2 - a);


	
	
	return a;
	
}



Matrix Tessellation::matrixCross(Matrix &m, Vector &v){
	Vector c(3);
	Vector r(3);
	Matrix res(3,3);
	int i,j;
	//produce a matrix in which each column is the crossproduct with the corresponding column in m and v
	
	for(i=0;i<3;++i){
		for(j=0;j<3;++j){
			c(j) = m(j,i);			
		}
		r = cross(c,v);
		for(j=0;j<3;++j){
			res(j,i) = r(j);			
		}
	}
	
	return res;
}





Rotation Tessellation::calculateOmega(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J){
	Vector ex(3);
	Vector ez(3);
	Vector ni(3);
	Vector nij(3);
	Vector nii(3);
	double sn;
	double varpi;
	double rho;
	Matrix Identity(3,3),reverseIdentity(3,3);
	Identity.eye();
	reverseIdentity.eye();
	reverseIdentity = reverseIdentity * -1;
	Rotation omega;
	double dot_ni_nij;
	Matrix dmui_dxi, dmui_dxl, dmuj_dxj, dmuj_dxl;
	double domega_dvarpi;
	Vector dvarpi_dnij, dvarpi_dni;
	Matrix dni_dmui, dnij_dmui, dnij_dmuj;
	Vector domega_dxi, domega_dxj, domega_dxl;
	
	ez(0)=0;
	ez(1)=0;
	ez(2)=1;
	double d;
	double s0,s1;

	
	d = dot(tessellationOrigin,I.normal);
	
	if(isWithinNumericalLimits(d,1.0)){
		ni = Vector(3);
		ni(0) = 0;
		ni(1) = 0;
		ni(2) = 1;
		printf("POSITIVE RANGE +1 %f\n",d);
	}
	else if(isWithinNumericalLimits(d,-1.0)){
		ni = Vector(3);
		ni(0) = 0;
		ni(1) = 0;
		ni(2) = -1;
		printf("POSITIVE RANGE -1 %f\n",d);
	}
	else ni = cross(tessellationOrigin,I.normal);
	
	nii = cross(I.normal,ni);
	
	nij = cross(I.normal, J.normal);
	dot_ni_nij = norm_dot(nij,nii);
	
	s0 = sgn(dot_ni_nij);
	varpi = asin(norm_dot(ni,nij));
	s1 = sgn(varpi);
	
	omega.rotation = -s0*(s1*(1-s0)*M_PI/2 - varpi);
	
	if(I.form != CONVEX){
		if(omega.rotation>=0)
			omega.rotation = -M_PI + omega.rotation;
	}
	
	
	dmui_dxi = I.dmu_dx;
	dmui_dxl = -I.dmu_dx;
	dmuj_dxj = J.dmu_dx;
	dmuj_dxl = -J.dmu_dx;
	
	domega_dvarpi = s0;
	dvarpi_dni = nij / sqrt(1-dot_ni_nij*dot_ni_nij);
	dvarpi_dnij = ni / sqrt(1-dot_ni_nij*dot_ni_nij);
	
	dni_dmui=Matrix(3,3);
	dni_dmui.zeros();
	dni_dmui(1,2)=1;
	dni_dmui(2,1)=-1;
	
	
	dnij_dmui = matrixCross(Identity,J.normal);
	dnij_dmuj = matrixCross(reverseIdentity,I.normal);
	
	domega_dxi = domega_dvarpi * ((dvarpi_dni.t() * dni_dmui * dmui_dxi).t() + (dvarpi_dnij.t() * dnij_dmui * dmui_dxi).t() );
	domega_dxj = domega_dvarpi * (dvarpi_dnij.t() * dnij_dmuj * dmuj_dxj).t();

	//domega_dxdelta = -domega_dxi-domega_dxj;

	domega_dxl = domega_dvarpi * ( (dvarpi_dni.t() * dni_dmui * dmui_dxi).t() + (dvarpi_dnij.t() * (dnij_dmui * dmui_dxl + dnij_dmuj * dmuj_dxl) ).t() );
	
	omega.drotation_dxi = domega_dxi;
	omega.drotation_dxj = domega_dxj;
	omega.drotation_dxl = domega_dxl;
	
	
	return omega;
	
}




Rotation Tessellation::calculateEta(CircularInterface &I, CircularInterface &J){
	double lambda_i, lambda_j, lambda_k, rho;
	double dot_IJ;
	Vector drho_mui(3), drho_muj(3);
	Vector drho_dxi(3), drho_dxj(3);
	Vector drho_dxl(3);
	double sig0,sig1,sig2;
	double deta_dlambdai, deta_dlambdaj, deta_drhoij;
	Vector deta_dxi(3), deta_dxj(3), deta_dxl(3);
	Rotation eta;
	Matrix dmui_dxi, dmuj_dxj;
	
	
	lambda_i = I.lambda.rotation;
	lambda_j = J.lambda.rotation;
	dot_IJ = norm_dot(I.normal,J.normal);
	rho = acos(dot_IJ);
	
	
	sig0 = cot(lambda_i) * cot(rho);
	sig1 = cos(lambda_j)*csc(lambda_i)*csc(rho);
	sig2 = sig0-sig1;
	
	eta.rotation = -asin(sig2);
	
	
	
	drho_mui = -J.normal/sqrt(1-dot_IJ*dot_IJ);
	drho_muj = -I.normal/sqrt(1-dot_IJ*dot_IJ);
	
	dmui_dxi = I.dmu_dx;
	dmuj_dxj = J.dmu_dx;
	
	printf("MATRIX FORM %d\n",J.form);
	dmuj_dxj.print("dmuj_dxj");
	
	drho_dxi = (drho_mui.t() * dmui_dxi).t();
	drho_dxj = (drho_muj.t() * dmuj_dxj).t();
	drho_dxl = -drho_dxi - drho_dxj;
	
	deta_dlambdai = (csc(lambda_i) * sig0/cos(lambda_i) - sig1*cos(lambda_i)) / sqrt(1-sig2*sig2);
	deta_dlambdaj = - sig1*tan(lambda_j) / sqrt(1-sig2*sig2);
	deta_drhoij = (csc(rho) * sig0/cos(rho) - sig1*cos(rho)) / sqrt(1-sig2*sig2);
	
	
	deta_dxi = deta_dlambdai * I.lambda.drotation_dxi + deta_drhoij * drho_dxi;
	deta_dxj = deta_dlambdaj * J.lambda.drotation_dxj + deta_drhoij * drho_dxj;
	deta_dxl = -deta_dxi-deta_dxj;
	
	eta.drotation_dxi = deta_dxi;
	eta.drotation_dxj = deta_dxj;
	eta.drotation_dxl = deta_dxl;
	
	
	return eta;
}



PHIContainer Tessellation::calculatePHI(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J, double dij, double radius){
	PHIContainer p,p2;
	
	Rotation eta,omega;
	
	
	//eta = M_PI/2 - acos(-(1/tan(lambda_j))*(1/tan(rho))+cos(lambda_k)*(1/sin(lambda_j))*(1/sin(rho)));
	
	eta = calculateEta(I,J);
	omega = calculateOmega(tessellationOrigin,I,J);
	

	
	p.out.rotation = eta.rotation + omega.rotation;
	if(p.out.rotation>=M_PI) p.out.rotation = (eta.rotation + omega.rotation) - 2*M_PI;
	else if(p.out.rotation<=-M_PI) p.out.rotation = (eta.rotation + omega.rotation) + 2*M_PI;
	p.out.drotation_dxi = eta.drotation_dxi + omega.drotation_dxi;
	p.out.drotation_dxj = eta.drotation_dxj + omega.drotation_dxj;
	p.out.drotation_dxl = eta.drotation_dxl + omega.drotation_dxl;
	
	
	
	p.in.rotation = (-M_PI-eta.rotation) + omega.rotation;
	if(p.in.rotation<=-M_PI) p.in.rotation = (M_PI-eta.rotation) + omega.rotation;
	p.in.drotation_dxi = -eta.drotation_dxi + omega.drotation_dxi;
	p.in.drotation_dxj = -eta.drotation_dxj + omega.drotation_dxj;
	p.in.drotation_dxl = -eta.drotation_dxl + omega.drotation_dxl;
		
	
	p2 = retrieveInterfaces(tessellationOrigin, I, J, dij, radius);
	return p2;
	
	if(I.form!=SPLITTER)
		if(!isWithinNumericalLimits(p.in.rotation,p2.in.rotation) || !isWithinNumericalLimits(p.out.rotation,p2.out.rotation)){
			printf("ROTATIONAL ERROR: in p:%f, p2:%f  out: p:%f p2:%f\n",p.in.rotation, p2.in.rotation, p.out.rotation, p2.out.rotation);
			exit(-1);
		}

	
	return p;
	
	
	
}
















Rotation Tessellation::calculatePsi(Vector tessellationOrigin, CircularInterface &circle){
	Vector dpsi_dmui;
	Rotation r;
	Vector dpsi_dxi, dpsi_dxl;
	
	
	if(circle.form != SPLITTER){
		dpsi_dmui = -tessellationOrigin/sqrt(1-circle.normal(0)*circle.normal(0));
		dpsi_dxi = (dpsi_dmui.t() * circle.dmu_dx).t();
		dpsi_dxl = -dpsi_dxi;
		
		r.rotation = getAngleBetweenNormals(tessellationOrigin,circle.normal);
		r.drotation_dxi = dpsi_dxi;
		r.drotation_dxj = Vector(3).zeros();
		r.drotation_dxl = dpsi_dxl;
	}
	else{
		//the splitter is immovable
		r.rotation = getAngleBetweenNormals(tessellationOrigin,circle.normal);
		r.drotation_dxi = Vector(3).zeros();
		r.drotation_dxj = Vector(3).zeros();
		r.drotation_dxl = Vector(3).zeros();
	}
	
	return r;
	
}






IntersectionBranches::iterator Tessellation::increaseBranchInterator(IntersectionBranches::iterator it, CircularInterfacesPerAtom &circles){
	IntersectionBranches* p;
	p = it->second.body;
	printf("FORM: %d id: %d\n",circles[it->second.id-1].form,it->second.id);
	if(circles[it->second.id-1].form == CONVEX){
		++it;
		if(it == p->end()) it=p->begin();
		
	}
	else{
		if(it == p->begin()) it=p->end();
		--it;
	}
	
	return it;
	
}

IntersectionBranches::iterator Tessellation::decreaseBranchInterator(IntersectionBranches::iterator it, CircularInterfacesPerAtom &circles){
	IntersectionBranches* p;
	p = it->second.body;
	printf("FORM: %d id: %d\n",circles[it->second.id-1].form,it->second.id);
	
	if(circles[it->second.id-1].form == CONVEX){
		if(it == p->begin()) it=p->end();
		--it;
		
	}
	else{
		++it;
		if(it == p->end()) it=p->begin();
	}
	
	return it;
}



IntersectionBranches::iterator Tessellation::increaseBranchInterator(multimap<double, IntersectionBranch>::iterator it, int ignore, CircularInterfacesPerAtom &circles){
	IntersectionBranches* p;
	p = it->second.body;
	printf("FORM: %d id: %d\n",circles[it->second.id-1].form,it->second.id);
	
	if(circles[it->second.id-1].form == CONVEX){
		++it;
		if(it == p->end()) it=p->begin();
		
	}
	else{
		if(it == p->begin()) it=p->end();
		--it;
	}
	
	if(it->second.it->second.id == ignore) it = increaseBranchInterator(it, ignore, circles);
	
	return it;
	
}

IntersectionBranches::iterator Tessellation::decreaseBranchInterator(IntersectionBranches::iterator it, int ignore, CircularInterfacesPerAtom &circles){
	IntersectionBranches* p;
	p = it->second.body;
	printf("FORM: %d id: %d\n",circles[it->second.id-1].form,it->second.id);
	
	if(circles[it->second.id-1].form == CONVEX){
		if(it == p->begin()) it=p->end();
		--it;
		
	}
	else{
		++it;
		if(it == p->end()) it=p->begin();
	}
	
	//printf("chk %d against %d\n",it->second.it->second.id, ignore);
	if(it->second.it->second.id == ignore) it = decreaseBranchInterator(it, ignore, circles);
	
	
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

void Tessellation::deleteIntersectionPoint(IntersectionBranches::iterator &it,IntersectionGraph &intersectionGraph, CircularInterfacesPerAtom &circles){
	//delete point from intersectionGraph
	IntersectionAddress address;
	IntersectionBranches::iterator it_mirror;
	int id,id_mirror;
	IntersectionBranches *body,*body_mirror;
	
	
	
	
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
	
	it_mirror = it->second.it;
	
	id = it->second.id;
	id_mirror = it_mirror->second.id;
	body = it->second.body;
	body_mirror = it_mirror->second.body;
	
	
	
	//first destroy the mirror
	body_mirror->erase(it_mirror);
	
	//then destroy the point
	body->erase(it);
	//from now on "it" is invalid...
	
	
	//if all branches have been deleted, the circular interface is completely interior and can be safely disregarded
	if(body->size()==0){
		circles[id-1].valid = false;
		printf("INVALIDATING %d\n",id);		
	}
	if(body_mirror->size()==0){
		circles[id_mirror-1].valid = false;
		printf("INVALIDATING %d\n",id_mirror);		
	}
	
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





void Tessellation::createIntersectionBranch(IntersectionAddress &address, PHIContainer &PHII, PHIContainer &PHIJ, CircularInterface &I, CircularInterface &J, IntersectionGraph &intersectionGraph){
	
	pair<double, IntersectionBranch> x;
	IntersectionBranches::iterator it0;
	IntersectionBranches::iterator it1;
	
	x.second.visited = -1;
	
	x.first = PHIJ.in.rotation;
	x.second.node=&intersectionGraph[address];
	x.second.direction = IN;
	x.second.it = I.intersectionBranches.end();
	x.second.body = &J.intersectionBranches;
	x.second.id = address.id1;
	x.second.flagged = false;
	it0 = J.intersectionBranches.insert(x);
	
	//printf("address %d-%d side %d IN\n",address.id0, address.id1, address.id1);
	
	x.first = PHII.out.rotation;
	x.second.node=&intersectionGraph[address];
	x.second.direction = OUT;
	x.second.it = it0;
	x.second.body = &I.intersectionBranches;
	x.second.id = address.id0;
	x.second.flagged = false;
	it1 = I.intersectionBranches.insert(x);
	it0->second.it=it1;

	//printf("address %d-%d side OUT\n",address.id0, address.id1, address.id0);
	
	
	intersectionGraph[address].rotation0=PHII.out;
	intersectionGraph[address].rotation1=PHIJ.in;
	
	
	
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



bool Tessellation::addToEraseList(map<IntersectionBranches::iterator, bool, IteratorComparator> &masterEraseList, map<IntersectionBranches::iterator,bool,IteratorComparator> &eraseList, IntersectionBranches::iterator &it, int limit){
	IntersectionBranches::iterator it_mirror;
	it_mirror = it->second.it;
	
	printf("[%d-%d] id!=limit: %d mirror_id!=limit: %d flagged: %d mirror_flagged: %d masterfind: %d eraselistfind: %d\n",it->second.node->id0,it->second.node->id1,it->second.id != limit,  it_mirror->second.id != limit , !it->second.flagged ,  !it_mirror->second.flagged , masterEraseList.find(it_mirror)==masterEraseList.end()  , eraseList.find(it_mirror)==eraseList.end());
	
	if(it->second.id != limit && it_mirror->second.id != limit && !it->second.flagged &&  !it_mirror->second.flagged && masterEraseList.find(it_mirror)==masterEraseList.end()  && eraseList.find(it_mirror)==eraseList.end()){
			eraseList[it]=true;
			it->second.flagged=true;
			it_mirror->second.flagged=true;
			printf("adding %d-%d\n",it->second.node->id0,it->second.node->id1);
	}
	else return false;
	
	return true;
}


bool Tessellation::addToEraseListCascade(map<IntersectionBranches::iterator, bool, IteratorComparator> &masterEraseList, map<IntersectionBranches::iterator, bool, IteratorComparator> &eraseList, IntersectionBranches::iterator &it, int limit, CircularInterfacesPerAtom &circles){
	IntersectionBranches::iterator it2, it_mirror;
	bool res;
	res = false;
	
	it_mirror = it->second.it;
	
	printf("TRYING TO ADD %d %d\n",it->second.node->id0,it->second.node->id1);
	if(it->second.id != limit && it_mirror->second.id != limit){

		printf("spreading %d-%d limit(%d)\n",it->second.node->id0,it->second.node->id1,limit);
	
		it2 = increaseBranchInterator(it, circles);
		res = addToEraseList(masterEraseList, eraseList, it2, limit) || res;

		it2 = decreaseBranchInterator(it, circles);
		res = addToEraseList(masterEraseList, eraseList, it2, limit) || res;
		

		it2 = increaseBranchInterator(it_mirror, circles);
		res = addToEraseList(masterEraseList, eraseList, it2, limit) || res;

		it2 = decreaseBranchInterator(it_mirror, circles);
		res = addToEraseList(masterEraseList, eraseList, it2, limit) || res;
	}
	else{
		printf("trying to spread %d-%d limit(%d)\n",it->second.node->id0,it->second.node->id1,limit);
		
		if(it->second.direction == IN){
			printf("IN\n");
			it2 = decreaseBranchInterator(it, circles);
			res = addToEraseList(masterEraseList, eraseList, it2, limit) || res;
		}
		else{
			printf("OUT\n");
			it2 = increaseBranchInterator(it, circles);
			res = addToEraseList(masterEraseList, eraseList, it2, limit) || res;
		}
		
		if(it_mirror->second.direction == IN){
			printf("MIRROR IN\n");
			it2 = decreaseBranchInterator(it_mirror, circles);
			res = addToEraseList(masterEraseList, eraseList, it2, limit) || res;
		}
		else{
			printf("MIRROR OUT\n");
			it2 = increaseBranchInterator(it_mirror, circles);
			res = addToEraseList(masterEraseList, eraseList, it2, limit) || res;
		}
		
	}
	
	return res;
		
}

void Tessellation::buildIntersectionGraph(double radius, Vector &tessellationOrigin, CircularInterfacesPerAtom &circles, SASAs &sasas, Hemisphere hemisphere, string filename){
	map<int, bool> processed;
	map<int, bool>::iterator it_p;
	IntersectionGraph intersectionGraph;
	IntersectionGraph::iterator it_g, it_x;
	int cid0, cid1;
	int q;
	
	SASA potentialSasa;
	SASANode sasaNode;
	
	int i,j;
	CircularInterface *I, *J;
	double tau0, tau1;
	map<int,CircularIntersection>::iterator it_j;
	
	bool empty;
	IntersectionAddress start, x, t;
	bool valid;
	map<IntersectionBranches::iterator,bool,IteratorComparator> eraseList, eraseList2, eraseList3;
	map<IntersectionBranches::iterator,bool,IteratorComparator>::iterator it_e;
	bool addedToEraseList;
	IntersectionPair ip;
	
	
	Interfaces interfacesJ, interfacesI;
	IntersectionBranches::iterator it_main;
	IntersectionBranches::iterator it0, it1;
	IntersectionBranches::iterator it, it_mirror, it_next, it_prev, it_mirror_next, it_mirror_prev, it_mirror_next_ignore, it_mirror_prev_ignore;
	IntersectionAddress addressIJ, addressJI;
	PHIContainer PHIJ, PHII;
	
	
	
	
	//iterate through all circles and add them to the intersectiongraph, one by one
	for(i=0; i < circles.size(); ++i){
		
		I = &circles[i];
		printf("PROCESSING %d id: %d\n",i,I->id);
		processed[i]=true;
		
		if(I->form!=CONVEX) printf("I FORM CONCAVE\n");
		
		
		
		//we have to sort this circle into the intersectiongraph
		//go through all intersecting circles, check whether they have been processed, and if so, calculate the intersections
		
		if(I->circularIntersections.size()==0 && I->form!=SPLITTER){
			insertArtificialIntersectionPoints(*I,intersectionGraph,tessellationOrigin);			
		}
		
		
		for(it_j=I->circularIntersections.begin(); it_j != I->circularIntersections.end(); ++it_j){
			j=(*it_j).first;
			J = &circles[j];
			printf("VALIDNESS: [%d]: %d\n",J->id,J->valid);
		}
		
		for(it_j=I->circularIntersections.begin(); it_j != I->circularIntersections.end(); ++it_j){
			j=(*it_j).first;
			J = &circles[j];
			printf("CHECKING CIRCLE %d id: %d\n",j,J->id);
			if(J->valid){
				
				
				
				
				
				it_p = processed.find(j);
				if(it_p != processed.end()){
					//printf("processing intersection with %d\n",j);
					
					if(J->form!=CONVEX) printf(" J FORM CONCAVE\n");
					
					addressIJ.id0 = I->id;
					addressIJ.id1 = J->id;
					addressJI.id0 = J->id;
					addressJI.id1 = I->id;
					
					//push the intersection points
					createIntersectionNode(addressIJ,intersectionGraph);
					createIntersectionNode(addressJI,intersectionGraph);
					
					
					//retrieve external and internal interfaces (respectively)
					PHIJ = calculatePHI(tessellationOrigin, *J, *I, it_j->second.d, radius);
					PHII = calculatePHI(tessellationOrigin, *I, *J, it_j->second.d, radius);
					
					createIntersectionBranch(addressIJ, PHII, PHIJ, *I, *J, intersectionGraph);
					createIntersectionBranch(addressJI, PHIJ, PHII, *J, *I, intersectionGraph);
					
				}
			}
		}
		
		for(int m=0; m<circles.size(); ++m){
			printf("circle %d\n",m);
			for(it_main = circles[m].intersectionBranches.begin(); it_main != circles[m].intersectionBranches.end(); ++it_main){
				if(it_main->second.direction==OUT)
					printf("[%d] (%d,%d) d: %f - (%f,%f) (OUT)\n",it_main->second.id, it_main->second.node->id0, it_main->second.node->id1, it_main->first, it_main->second.node->rotation0.rotation, it_main->second.node->rotation1.rotation);
				else
					printf("[%d] (%d,%d) d: %f - (%f,%f) (IN)\n",it_main->second.id, it_main->second.node->id0, it_main->second.node->id1, it_main->first, it_main->second.node->rotation0.rotation, it_main->second.node->rotation1.rotation);
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
				
				printBranch("it",it);
				printBranch("it_mirror",it_mirror);
				
				it_next = increaseBranchInterator(it, circles);
				printBranch("it_next",it_next);
				
				it_prev = decreaseBranchInterator(it, circles);
				printBranch("it_prev",it_prev);
				
				it_mirror_next = increaseBranchInterator(it_mirror, circles);
				printBranch("it_mirror_next",it_mirror_next);
				it_mirror_prev = decreaseBranchInterator(it_mirror, circles);
				printBranch("it_mirror_prev",it_mirror_prev);
				
				
				//if the mirror was previously empty, then we just added 2 branches. Hence, if it has two branches, it's considered empty.
				if(it_mirror->second.body->size()<=2) empty=true;
				else empty=false;
				
				
				
				if(!empty){
					it_mirror_next_ignore = increaseBranchInterator(it_mirror, it->second.id, circles);
				printBranch("it_mirror_next_ignore",it_mirror_next_ignore);
					it_mirror_prev_ignore = decreaseBranchInterator(it_mirror, it->second.id, circles);
				printBranch("it_mirror_prev_ignore",it_mirror_prev_ignore);
				}
				else{
					//printf("EMPTY\n");
					it_mirror_next_ignore = increaseBranchInterator(it_mirror, circles);
				printBranch("it_mirror_next_ignore",it_mirror_next_ignore);
					it_mirror_prev_ignore = decreaseBranchInterator(it_mirror, circles);
				printBranch("it_mirror_prev_ignore",it_mirror_prev_ignore);
				}
				
				
				
				
				
				
				printf("it_mirror_prev ig forward %d\n",it_mirror_prev_ignore->second.forward);
				printf("it_mirror_prev ig backward %d\n",it_mirror_prev_ignore->second.backward);
				printf("it_mirror_next ig forward %d\n",it_mirror_next_ignore->second.forward);
				printf("it_mirror_next ig backward %d\n",it_mirror_next_ignore->second.backward);
				
				if(empty){
					printf("EMPTY\n");
					q = 1;
					if(it->second.direction == IN) q*=-1;
					printf("Q: %d\n",q);
					//if(circles[it->second.id-1].form != CONVEX) q*=-1;
					//printf("Q: %d\n",q);
					//if(circles[it_mirror->second.id-1].form != CONVEX) q*=-1;
					//printf("Q: %d\n",q);
					
					
					
					
					if(q>0){
						printf("POSITIVE Q\n");
						it->second.forward = INTERNAL;
						it->second.backward = EXTERNAL;
						it_mirror->second.forward = EXTERNAL;
						it_mirror->second.backward = INTERNAL;
					}
					else{
						printf("NEGATIVE Q\n");
						it->second.forward = EXTERNAL;
						it->second.backward = INTERNAL;
						it_mirror->second.forward = INTERNAL;
						it_mirror->second.backward = EXTERNAL;
					}
					connectIntersectionPoints(*it_mirror->second.node, *it_next->second.node, intersectionGraph);
				}
				else{
					if(it->second.direction == IN){
						printf("IN\n");
						
						if(it_mirror_prev_ignore->second.forward == EXTERNAL && it_mirror_next_ignore->second.backward == EXTERNAL){
							printf("EXTERNAL\n");
							connectIntersectionPoints(*it_mirror_prev->second.node, *it_mirror->second.node, intersectionGraph);
							connectIntersectionPoints(*it_mirror->second.node, *it_next->second.node, intersectionGraph);
							it->second.backward=INTERNAL;
							it->second.forward=EXTERNAL;
							it_mirror->second.backward=EXTERNAL;
							it_mirror->second.forward=INTERNAL;
							
							if(it_mirror_next->second.it->second.id != it->second.id){
								printf("ERASING [0] %d %d\n",it_mirror_next->second.node->id0,it_mirror_next->second.node->id1);
								eraseList[it_mirror_next]=true;
							}
							if(it_prev->second.it->second.id != it_mirror->second.id){
								printf("ERASING [1] %d %d\n",it_prev->second.node->id0,it_prev->second.node->id1);
								eraseList[it_prev]=true;
							}
							
						}
						else{
							printf("INTERNAL\n");
							printf("ERASING [2] %d %d\n",it->second.node->id0,it->second.node->id1);
							eraseList[it]=true;
						}
					}
					
					if(it->second.direction == OUT){
						printf("OUT\n");
						
						if(it_mirror_next_ignore->second.backward == EXTERNAL && it_mirror_prev_ignore->second.forward == EXTERNAL){
							printf("EXTERNAL\n");
							connectIntersectionPoints(*it_mirror->second.node, *it_mirror_next->second.node, intersectionGraph);
							connectIntersectionPoints(*it_prev->second.node, *it->second.node, intersectionGraph);
							it->second.backward=EXTERNAL;
							it->second.forward=INTERNAL;
							it_mirror->second.backward=INTERNAL;
							it_mirror->second.forward=EXTERNAL;
							
							if(it_mirror_prev->second.it->second.id != it->second.id){
								printf("ERASING [0] %d %d\n",it_mirror_prev->second.node->id0,it_mirror_prev->second.node->id1);
								eraseList[it_mirror_prev]=true;
							}
							if(it_next->second.it->second.id != it_mirror->second.id){
								printf("ERASING [1] %d %d\n",it_next->second.node->id0,it_next->second.node->id1);
								eraseList[it_next]=true;
							}
							
						}
						else{
							printf("INTERNAL\n");
							printf("ERASING [2] %d %d\n",it->second.node->id0,it->second.node->id1);
							eraseList[it]=true;
						}
					}
					
					
					/*
					
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
					*/
					
					
				}
				
				printf("-------------------------------------------------\n");
				printIntersectionGraph(intersectionGraph);
				printf("-------------------------------------------------\n");
				
			}
			
		}
		
		printf("CASCADE\n");
		
		
		for(it_e=eraseList.begin(); it_e != eraseList.end(); ++it_e){
			it = it_e->first;
			eraseList2[it]=true;
		}
		
		do{
			printf("-------MASTER LIST--------\n");
			for(it_e=eraseList.begin(); it_e != eraseList.end(); ++it_e){
				it = it_e->first;
				printf("%d-%d\n",it->second.node->id0,it->second.node->id1);
			}
			printf("-------------------------------");
			printf("-------SEARCH LIST--------\n");
			for(it_e=eraseList2.begin(); it_e != eraseList2.end(); ++it_e){
				it = it_e->first;
				printf("%d-%d\n",it->second.node->id0,it->second.node->id1);
			}
			printf("-------------------------------");
			
			printf("ERASING %d\n",eraseList2.size());
			addedToEraseList = false;
			eraseList3.clear();
			
			int b2=0;
			for(it_e=eraseList2.begin(); it_e != eraseList2.end(); ++it_e){
				it = it_e->first;
				printf("B: %d\n",b2);
				++b2;
				addedToEraseList = addToEraseListCascade(eraseList, eraseList3, it, I->id, circles) || addedToEraseList; //order of the or operands is important
				
			}
			
			for(it_e=eraseList3.begin(); it_e != eraseList3.end(); ++it_e){
				it = it_e->first;
				eraseList[it]=true;
			}
			eraseList2 = eraseList3;
			
			
		}while(addedToEraseList);
		

		
		for(it_e=eraseList.begin(); it_e != eraseList.end(); ++it_e){
			it = it_e->first;
			printf("REMOVING IP: %d\n",&*it);
			printf("--: %d-%d\n",it->second.node->id0,it->second.node->id1);
			deleteIntersectionPoint(it,intersectionGraph,circles);
		}
		
		
		printf("-------------------------------------------------\n");
		printIntersectionGraph(intersectionGraph);
		printf("-------------------------------------------------\n");
		
		
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
				cid1=abs(x.id1)-1;
				
				sasaNode.id0 = x.id0;
				sasaNode.id1 = x.id1;
				sasaNode.index0 = circles[cid0].index;
				sasaNode.index1 = circles[cid1].index;;
				sasaNode.rotation0 = intersectionGraph[x].rotation0;
				sasaNode.rotation1 = intersectionGraph[x].rotation1;
				sasaNode.lambda = circles[cid0].lambda;
				sasaNode.psi = circles[cid0].psi;
				sasaNode.vector = intersectionGraph[x].rotation1.vector;
				sasaNode.normalForCircularInterface = circles[cid0].normal;
				sasaNode.form = circles[cid0].form;
				
				double ank = getAngle(circles[cid0].normal, circles[cid1].normal);
				//printf("VERIFICATION.. g0: %f g1: %f cosrho: %f PHI0: %f PHI1: %f rotational part: %f n0(%f, %f, %f) n1(%f, %f, %f)\n",circles[cid0].g/radius,circles[cid1].g/radius, cos(ank), sasaNode.angle0,  sasaNode.angle1, rotationalAngle(tessellationOrigin, circles[cid0],circles[cid1]),circles[cid0].normal(0),circles[cid0].normal(1),circles[cid0].normal(2),circles[cid1].normal(0),circles[cid1].normal(1),circles[cid1].normal(2)); 
				
				
				
				
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












void Tessellation::outputGaussBonnetData(string filename, double radius, CircularInterfacesPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph){
	FILE* file;
	
	file = fopen (filename.c_str(),"w");
	
	SASA sasa;
	SASANodeList::iterator it;
	CircularInterface circle;
	double area=0;
	int k;
	Vector origin(3), v(3);
	
	
	fprintf(file, "atom radius %f\n", radius);
	
	
	for(int i=0; i< circles.size(); ++i){
		circle = circles[i];
		if(circle.form==SPLITTER){
			fprintf(file, "circularregion %d radius %f vector %f %f %f form %d\n", circle.id, circle.a, 0.0001, 0.0, 0.0, 2);
		}
		else{
			fprintf(file, "circularregion %d radius %f vector %f %f %f form %d\n", circle.id, circle.a, circle.normal(0)*circle.g, circle.normal(1)*circle.g, circle.normal(2)*circle.g, circle.form);
			fprintf(file, "intersector %d radius %f vector %f %f %f\n", i, circle.sphereRadius, circle.vector(0), circle.vector(1), circle.vector(2));
		}
	} 
	/*
	origin = atoms[13];
	
	for(int i=0; i< atoms.size(); ++i){
		if(i!=13){
			v = atoms[i] - origin;
			fprintf(file, "intersector %d radius %f vector %f %f %f\n", i, radii[i], v(0), v(1), v(2));
		}
	}
	*/
	
	
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
