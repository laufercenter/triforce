#include "tessellation.h"


Tessellation::Tessellation(Molecule &m){
	molecule = m;
}



void Tessellation::build(bool split){
	CircularInterfacesPerAtom circlesPerAtom;
	vector<int> neighbourlist;
	sasasForMolecule.clear();
	
	//molecule.update();
	atoms = molecule.fetchCoordinates();
	radii = molecule.fetchRadii();
	
	
	
	//iterate over all atoms and build the tessellation for each of them
	for(unsigned int i=0; i<atoms.size(); ++i){
		//int i=12;{
		neighbourlist = molecule.getNeighborListFor(i);

		
		//printf("NEIGHBOURS RECEIVED: %d\n",neighbourlist.size());
		/*vector<int>::iterator it;
		for(it=neighbourlist.begin(); it!=neighbourlist.end(); ++it){
			printf(" %d",*it);
		}
		printf("\n");
		*/
		
		//printf("NEIGHBOURS RECEIVED: %d\n",neighbourlist.size());
		buildGaussBonnetPath(i, atoms, radii, sasasForMolecule, split, neighbourlist);
		/*
		if(i==10){
		fprintf(stderr,"--GBATOM-- %d\n",i);
			for(int j=0; j<sasasForMolecule[i].sasas.size(); ++j){
				fprintf(stderr,"SASA %d (%d)\n",j,sasasForMolecule[i].sasas[j].hemisphere);
				outputGaussBonnetPath(sasasForMolecule[i].sasas[j]);
			}
		}
		*/

	}
	
	
	
}



CircularInterfacesPerAtom Tessellation::coverHemisphere(Vector tessellationAxis, double radius, CircularInterfacesPerAtom circles, CircularInterfaceForm form){
	CircularInterface C;
	Vector v(3);
	Matrix NullMatrix(3,3);
	NullMatrix.zeros();
	
	v = tessellationAxis;
	
	C.lambdaRotation.rotation = M_PI * 0.5;
	C.lambda.rotation = M_PI * 0.5;
	C.psi.rotation=0;
	C.g=0;
	C.normal = v;
	C.form = form;
	C.id = circles.size()+1;
	C.valid = true;
	C.index = -1;
	C.hasDerivatives=true;
	
	C.dmu_dx = NullMatrix;
	C.lambda.drotation_dxi = Vector(3).zeros();
	C.lambda.drotation_dxj = Vector(3).zeros();
	C.lambda.drotation_dxl = Vector(3).zeros();

	C.psi.drotation_dxi = Vector(3).zeros();
	C.psi.drotation_dxj = Vector(3).zeros();
	C.psi.drotation_dxl = Vector(3).zeros();
	
	
	circles.push_back(C);
	
	return circles;
	
}

void Tessellation::buildGaussBonnetPath(int i, vector<Vector> &atoms, vector<double> &radii, SASAs &sasas, bool split, vector<int> &neighbourlist){
	CircularInterfacesPerAtom circles;
	CircularInterfacesPerAtom circlesFrontHemisphere;
	CircularInterfacesPerAtom circlesBackHemisphere;
	
	
	Vector frontTessellationAxis(3);
	Vector backTessellationAxis(3);
	Vector origin;
	double radius;
	bool isInsideAnotherSphere;
	Vector v;
	
	
	origin = atoms[i];
	torigin = origin;
	radius = radii[i];
	tradius = radius;
	ti = i;
	
	
	
	frontTessellationAxis(0) = 1.0;
	frontTessellationAxis(1) = 0.0;
	frontTessellationAxis(2) = 0.0;
	
	/*
	v = randu<vec>(3);
	v(0) = v(0)-0.5;
	v(1) = v(1)-0.5;
	v(2) = v(2)-0.5;
	frontTessellationAxis=v;
	*/
	
	
	frontTessellationAxis = frontTessellationAxis/norm(frontTessellationAxis,2);
	
	
	


	backTessellationAxis = -frontTessellationAxis;
	SASASegmentList *newSasas;
	
	srand(2);
	
	SASASegmentList sasa;
	
	
	sasas.push_back(sasa);
	


	isInsideAnotherSphere = makeCircularInterfaces(i,origin, radius, atoms, radii, circles, neighbourlist);
	
	if(isInsideAnotherSphere) return;
	
	for(unsigned int j=0;j<circles.size();j++){
		determineProjection(origin, radius, circles[j]);
	}
	
	//there is room for optimisation here...
	
	if(split){
	
		circlesFrontHemisphere = coverHemisphere(frontTessellationAxis, radius, circles, SPLITTER);
		circlesBackHemisphere = coverHemisphere(backTessellationAxis, radius, circles, SPLITTER);
		
		
		filterCircularInterfaces(frontTessellationAxis,radius, circlesFrontHemisphere);
		filterCircularInterfaces(backTessellationAxis,radius, circlesBackHemisphere);
		
		//for(unsigned int j=0;j<circlesFrontHemisphere.size();j++) if(circlesFrontHemisphere[j].form!=SPLITTER) calculateProjectionDerivatives(circlesFrontHemisphere[j]);
		//for(unsigned int j=0;j<circlesBackHemisphere.size();j++) if(circlesBackHemisphere[j].form!=SPLITTER) calculateProjectionDerivatives(circlesBackHemisphere[j]);
		
		//determinePsiRotations(frontTessellationAxis, circlesFrontHemisphere);
		//determinePsiRotations(backTessellationAxis, circlesBackHemisphere);
		
		reindexCircularInterfaces(circlesFrontHemisphere);
		reindexCircularInterfaces(circlesBackHemisphere);
		
		determineCircularIntersections(circlesFrontHemisphere);
		determineCircularIntersections(circlesBackHemisphere);
		
		
		buildIntersectionGraph(radius, frontTessellationAxis, circlesFrontHemisphere, sasas[sasas.size()-1], FRONTHEMISPHERE, string("gbonnet0.csv"));
		buildIntersectionGraph(radius, backTessellationAxis, circlesBackHemisphere, sasas[sasas.size()-1], BACKHEMISPHERE, string("gbonnet1.csv"));
	}
	else{
		printf("not supported\n");
		exit(-1);
		
		filterCircularInterfaces(frontTessellationAxis, radius, circles);
		reindexCircularInterfaces(circles);
		
		determineCircularIntersections(circles);
		
		buildIntersectionGraph(radius, frontTessellationAxis, circles, *newSasas, FRONTHEMISPHERE, string("gbonnet0.csv"));
	}
	
	

	
}






SASAs &Tessellation::sasas(){
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
	if(eps-(v+THRESHOLD_NUMERICAL) <= 0) return true;
	else return false;
}

bool Tessellation::isWithinNumericalLimits(double x, double l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) return true;
	else return false;
}


bool Tessellation::isWithinStrongNumericalLimits(double x, double l){
	if(abs(x-l)<=THRESHOLD_STRONG_NUMERICAL) return true;
	else return false;
}


LambdaRotation Tessellation::calculateLambda(double d_i, double r_l, double r_i, Vector &mu_i){
	CircularInterfaceForm form;
	return calculateLambda(0,d_i, r_l, r_i, mu_i, form);
}


LambdaRotation Tessellation::calculateLambda(int index_i, double d_i, double r_l, double r_i, Vector &mu_i, CircularInterfaceForm &form){
	double g;
	double dg_dd;
	double dlambda_dg;
	LambdaRotation r;
	double q;
	Vector mu_i_original;
	Vector dlambda_dxi(3), dlambda_dxl(3);
	Vector dd_dxi;
	
	mu_i_original = mu_i;
	g = (d_i * d_i + r_l * r_l - r_i * r_i ) / (2 * d_i * r_l);

	if(g<0){
		form=CONCAVE;
		g = abs(g);
		mu_i=-mu_i;
	}
	else form=CONVEX;
	
	r.rotation = acos(g);
	r.d_i = d_i;
	r.r_l = r_l;
	r.g = g;
	
	return r;
	
}



Rotation Tessellation::calculateLambdaDerivatives(LambdaRotation &r, CircularInterface &circle){
	double g;
	double dg_dd;
	double dlambda_dg;
	double q;
	Vector mu_i_original;
	Vector dlambda_dxi(3), dlambda_dxl(3);
	Vector dd_dxi;
	Rotation lambda;
	double d_i;
	double r_l;
	double r_i;
	Vector mu_i;
	CircularInterfaceForm form;
	
	g = r.g;
	d_i = r.d_i;
	r_l = r.r_l;
	mu_i = circle.normal;
	form = circle.form;

	if(g > 1.0-MINISCULE) g = 1.0-MINISCULE;
	if(g < -1.0+MINISCULE) g = -1.0+MINISCULE;
	

	dlambda_dg = -1/sqrt(1-g*g);
	
	
	
	
	if(form==CONVEX){
		dg_dd = -g/d_i + 1/r_l;
		q = 1;
	}
	else{
		dg_dd = -(g/d_i + 1/r_l);
		q = -1;
	}
	

	
	dd_dxi = q*mu_i;
	
	dlambda_dxi = dlambda_dg * dg_dd * dd_dxi;
	//smoothing
	//if(norm(dlambda_dxi,2)>1.0) dlambda_dxi = normalise(dlambda_dxi);
	dlambda_dxl = -dlambda_dxi;
	
	lambda.rotation = r.rotation;
	lambda.drotation_dxi = dlambda_dxi;
	lambda.drotation_dxj = Vector(3).zeros();
	lambda.drotation_dxl = dlambda_dxl;
	
	
	
	
	
	return lambda;
	
	
}



void Tessellation::determineProjection(Vector &origin, double r_l, CircularInterface &circle){
	circle.lambdaRotation = calculateLambda(circle.index, circle.d, r_l, circle.sphereRadius, circle.normal, circle.form);
	circle.lambda.rotation = circle.lambdaRotation.rotation;
	circle.valid=true;
}



void Tessellation::calculateProjectionAndDerivatives(Vector &tessellationAxis, CircularInterface &circle){
	Matrix Identity(3,3);
	Matrix dmu_dx(3,3);
	Matrix fd_dmu_dx(3,3);
	Identity.eye();
	
	
	
	

	circle.lambda = calculateLambdaDerivatives(circle.lambdaRotation, circle);
	
	//derivatives
	if(circle.form==CONVEX)
		dmu_dx = (1.0/circle.d)*(Identity - kron(circle.normal,circle.normal.t())).t();
	else
		dmu_dx = -(1.0/circle.d)*(Identity - kron(circle.normal,circle.normal.t())).t();
	
	circle.dmu_dx = dmu_dx;
	
	
	
	
	circle.psi = calculatePsi(tessellationAxis, circle);
	
	circle.hasDerivatives=true;
	
	
	
}



void Tessellation::determinePsiRotations(Vector &tessellationAxis, CircularInterfacesPerAtom &circles){
	for(unsigned i=0; i < circles.size(); ++i)
		circles[i].psi=calculatePsi(tessellationAxis, circles[i]);
}


IntersectionPair Tessellation::determineIntersectionPoints(double radius, CircularInterface &K, CircularInterface &J){
	IntersectionPair res;
	double g_k, g_j;
	double r_i;
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


Vector Tessellation::calculateInterfaceNormal(const Vector &v_l, const Vector &v_i){
	return calculateInterfaceNormal(v_l,v_i);
}


Vector Tessellation::calculateInterfaceNormal(const Vector &v_l, const Vector &v_i, double &d){
	Vector v;
	v = v_i - v_l;
	d = norm(v,2);
	v = v/d;
	return v;
	
}


bool Tessellation::makeCircularInterfaces(int l,Vector &origin, double r_l, vector<vec> &atoms, vector<double> &radii, vector<CircularInterface> &circles, vector<int> &neighbourlist){
	CircularInterface circle;
	double r_i;
	double d_i;
	Vector mu_i;
	vector<int>::iterator it;
	int j;
	

	for(it=neighbourlist.begin(); it!=neighbourlist.end(); ++it){
		j=*it;
		if(l != j){
			r_i = radii[j];
			
			mu_i = calculateInterfaceNormal(origin, atoms[j], d_i);
			
			//reject, if no intersection
			if(d_i < r_l + r_i && d_i+r_i > r_l && d_i+r_l > r_i){
				circle.id = circles.size()+1;
				circle.normal = mu_i;
				circle.d = d_i;
				circle.sphereRadius = r_i;
				circle.intersect = false;
				circle.index=j;
				circle.hasDerivatives=false;
				circles.push_back(circle);
			}
			//in this case, the atom is completely inside another atom and has no SASA
			else if(d_i+r_l <= r_i){
				circles.clear();
				return true;
			}
			
		}
	}
	return false;
}



void Tessellation::depleteCircularInterfaces(Vector tessellationAxis, double radius, vector<CircularInterface> &circles){
	circles.clear();
	Vector t;
	t = -tessellationAxis;
	circles = coverHemisphere(t, radius, circles, CONCAVE);
}

int Tessellation::filterCircularInterfaces(Vector tessellationAxis, double radius, vector<CircularInterface> &circles){
	vector<CircularInterface>::iterator it;
	double angle;
	bool erased;
	Vector n0,n1;
	it = circles.begin();
	while(it != circles.end()){
		erased=false;
	
		
		n0 = it->normal;
		
		for(unsigned i=0;i<circles.size();i++){
			if(it->id != circles[i].id){
				
				n1 = circles[i].normal;
				
				angle = getAngle(n0,n1);
				
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
							circles.clear();
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
							circles.clear();
							return -1;
							
						}
						
					}

				}	
				
						
					
			}
		}
		if(!erased) ++it;
	}
	
	
	//if just the splitter is left, the hemisphere is completely free. In that case, the splitter is unnecessary
	//if(circles.size()==1 && circles[0].form==SPLITTER)
	//	circles.clear();
	
	return 0;

}






void Tessellation::outputGaussBonnetPath(SASAs &points){
	/*
	SASASegmentList::iterator it;
	int i;
	for(it=points.sasa.begin(), i=0; it!=points.sasa.end(); ++it, ++i)
		fprintf(stderr,"GBPATH[%d] %d - %d indexes: %d - %d form: %d %d\n", i, it->id0, it->id1, it->index0, it->index1, it->form0, it->form1);
	*/
}



void Tessellation::reindexCircularInterfaces(CircularInterfacesPerAtom &circles){
	for(unsigned i=0;i<circles.size();i++){
		circles[i].id = i+1;
	}
}





void Tessellation::insertArtificialIntersectionPoints(CircularInterface &I, Vector &tessellationAxis, Hemisphere hemisphere, SASASegmentList &sasa){
	Rotation f,f_reverse;
	SASASegment sasaSegment;
	
	
	
	
	//this will create two points on the border of the circular region on opposite sides.

	
	f.rotation = -M_PI/2;
	f.drotation_dxi=Vector(3).zeros();
	f.drotation_dxj=Vector(3).zeros();
	f.drotation_dxl=Vector(3).zeros();

	f_reverse.rotation = M_PI/2;
	f_reverse.drotation_dxi=Vector(3).zeros();
	f_reverse.drotation_dxj=Vector(3).zeros();
	f_reverse.drotation_dxl=Vector(3).zeros();
	
	
	
	
	calculateProjectionAndDerivatives(tessellationAxis, I);
	
	
	sasaSegment.id0 = I.id;
	sasaSegment.id1 = I.id;
	sasaSegment.id2 = I.id;
	
	sasaSegment.index0 = I.index;
	sasaSegment.index1 = I.index;
	sasaSegment.index2 = I.index;
	sasaSegment.lambda = I.lambda;
	
	
	sasaSegment.psi = I.psi;
	sasaSegment.normalForCircularInterface = I.normal;
	sasaSegment.form0 = I.form;
	sasaSegment.form1 = I.form;
	sasaSegment.form2 = I.form;
	
	sasaSegment.tessellationAxis = tessellationAxis;
	sasaSegment.hemisphere = hemisphere;
	
	
	sasaSegment.rotation0 = f;
	sasaSegment.rotation1 = f_reverse;
	
	sasa.push_back(sasaSegment);
	
	sasaSegment.rotation1 = f;
	sasaSegment.rotation0 = f_reverse;
	
	sasa.push_back(sasaSegment);
	
	
			
}




/*


void Tessellation::insertArtificialIntersectionPoints(CircularInterface &I, IntersectionGraph &intersectionGraph, Vector &tessellationAxis){
	Vector v,v2,o,o2;
	int k;
	IntersectionNode a,b;
	IntersectionAddress addressA,addressB;
	Vector x0,x1;
	Rotation f,f_reverse;
	int q=1;
	
	
	
	//this will create two points on the border of the circular region on opposite sides.
	//if(tessellationAxis(0)<0) q=-q;
	//if(I.form==CONVEX) q=-q;

	
	if(q<0){
		f.rotation = M_PI/2;
		f_reverse.rotation = -M_PI/2;
	}
	else{
		f.rotation = -M_PI/2;
		f_reverse.rotation = M_PI/2;
	}
	
	a.id0 = addressA.id0 = I.id;
	a.id1 = addressA.id1 = -I.id;
	a.pointsTo0 = a.prev0 = -I.id;
	a.pointsTo1 = a.prev1 = I.id;
	a.vector = x0;
	a.rotation0=f;
	a.rotation1=f;
	a.visited=false;
	a.rotation0.drotation_dxi=Vector(3).zeros();
	a.rotation0.drotation_dxj=Vector(3).zeros();
	a.rotation0.drotation_dxl=Vector(3).zeros();
	a.rotation1.drotation_dxi=Vector(3).zeros();
	a.rotation1.drotation_dxj=Vector(3).zeros();
	a.rotation1.drotation_dxl=Vector(3).zeros();

	b.id0 = addressB.id0 = -I.id;
	b.id1 = addressB.id1 = I.id;
	b.pointsTo0 = b.prev0 = I.id;
	b.pointsTo1 = b.prev1 = -I.id;
	b.vector = x1;
	b.rotation0=f_reverse;
	b.rotation1=f_reverse;
	b.visited=false;
	b.rotation0.drotation_dxi=Vector(3).zeros();
	b.rotation0.drotation_dxj=Vector(3).zeros();
	b.rotation0.drotation_dxl=Vector(3).zeros();
	b.rotation1.drotation_dxi=Vector(3).zeros();
	b.rotation1.drotation_dxj=Vector(3).zeros();
	b.rotation1.drotation_dxl=Vector(3).zeros();
	
	
	intersectionGraph[addressA]=a;
	intersectionGraph[addressB]=b;
	
			
}
*/





int Tessellation::sgn(double d){
	if(d>=0) return 1;
	else return -1;
}

void Tessellation::determineCircularIntersections(CircularInterfacesPerAtom &circles){
	double angle;
	RhoContainer c;
	
	//this doesn't really need to be optimised. The reason is, that the number of excluded interface-interface intersections is much smaller and at most equal to the number of included intersections
	for(unsigned k=0;k<circles.size();k++)
		for(unsigned j=k+1;j<circles.size();j++){
			if(k!=j){
				//determine if there will be intersections
				angle = getAngle(circles[k].normal,circles[j].normal);
				if(angle < circles[k].lambda.rotation + circles[j].lambda.rotation){
					if(angle + circles[k].lambda.rotation > circles[j].lambda.rotation && angle + circles[j].lambda.rotation > circles[k].lambda.rotation){
						c.rho=angle;
						c.id = k;
						circles[j].circularIntersections.push_back(c);
						c.id = j;
						circles[k].circularIntersections.push_back(c);
					}
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

/*
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

void Tessellation::measurementPoints(Vector &p0, Vector &p1, Vector &tessellationAxis, CircularInterface &I){
	Vector v(3);
	Vector o(3);
	Vector vx(3);
	Vector s(3), s2(3);
	
	o = tessellationAxis;
	v = I.normal * I.g;
	
	
	
	if(isWithinNumericalLimits(dot(o,I.normal),1.0) ){
		p1 = Vector(3);
		p1(0) = 0;
		p1(1) = -1;
		p1(2) = 0;
		
		p0 = cross(p1,-o);
		
	}
	else 	if(isWithinNumericalLimits(dot(o,I.normal),-1.0) ){
		p1 = Vector(3);
		p1(0) = 0;
		p1(1) = 1;
		p1(2) = 0;
		
		p0 = cross(p1,-o);

		
	}
	else{
		

		
		s = cross(v,o);
		s = I.a*s / norm(s,2);
		p0 = v+s;
		
		s2 = -cross(v,s);
		s2 = I.a*s2 / norm(s2,2);
		p1 = v+s2;
		
		if(dot(tessellationAxis,I.normal)<0){
			//p0 = v-s;
			//p1 = v-s2;
		}
		
		
	}
	
}

Interfaces Tessellation::angularInterfaces(Vector &x0, Vector &x1, Vector &tessellationAxis, CircularInterface &I, CircularInterface &J){
	Vector p0(3), p1(3);
	Interfaces res;
	Vector v(3);
	Vector t;
	int q;
	
	measurementPoints(p0,p1,tessellationAxis,I);
	v = I.normal * I.g;
	
	
	
	
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






PHIContainer Tessellation::retrieveInterfaces(Vector &tessellationAxis, CircularInterface &I, CircularInterface &J, double dij, double radius){
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
		
		
		intf1 = angularInterfaces(x0,x1, tessellationAxis, I, J);
		
		

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
	

		
	return p;
	
	
}




double Tessellation::rotationalAngle(Vector &tessellationAxis, CircularInterface &I, CircularInterface &J){
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

	
	d = dot(tessellationAxis,I.normal);
	
	if(isWithinNumericalLimits(d,1.0)){
		n0 = Vector(3);
		n0(0) = 0;
		n0(1) = 0;
		n0(2) = 1;
	}
	else if(isWithinNumericalLimits(d,-1.0)){
		n0 = Vector(3);
		n0(0) = 0;
		n0(1) = 0;
		n0(2) = -1;
	}
	else n0 = cross(tessellationAxis,I.normal);
	
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

*/

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

Vector Tessellation::normalise(Vector x){
	double l;
	Vector v(3);
	l = norm(x,2);
	v = x/l;
	return v;
	
}

Vector Tessellation::normalise(Vector x, double &l){
	Vector v(3);
	l = norm(x,2);
	v = x/l;
	return v;
	
}





double Tessellation::sacos(Vector &a, Vector &b){
	double la = norm(a,2);
	double lb = norm(b,2);
	return acos(dot(a,b)/(la*lb));
	
}


double Tessellation::sdot(Vector &a, Vector &b){
	double la = norm(a,2);
	double lb = norm(b,2);
	return dot(a,b)/(la*lb);
	
}

double Tessellation::l(Vector &a){
	double s = 0;
	for(int i=0; i<3; ++i)
		s += a(i)*a(i);
	return sqrt(s);
}


OmegaRotation Tessellation::calculateOmega(Vector &tessellationAxis, CircularInterface &I, CircularInterface &J, RhoContainer &rhoContainer){
	return calculateOmega(tessellationAxis, I.normal, J.normal, I.id-1, J.id-1, I.form, J.form, rhoContainer);

}


OmegaRotation Tessellation::calculateOmega(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, RhoContainer &rhoContainer){
	return calculateOmega(tessellationAxis, mu_i, mu_j, 0, 0, CONVEX, CONVEX, rhoContainer);
}



OmegaRotation Tessellation::calculateOmega(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer){
	Vector ex(3);
	Vector ez(3);
	Vector ey(3);
	Vector ni(3),ni2(3);
	Vector nij(3);
	Vector nii(3);
	double varpi;
	Matrix Identity(3,3),reverseIdentity(3,3);
	Identity.eye();
	reverseIdentity.eye();
	reverseIdentity = reverseIdentity * -1;
	OmegaRotation omega;
	double dot_ni_nij;
	bool problem=false;
	
	ez(0)=0;
	ez(1)=0;
	ez(2)=1;
	ey(0)=0;
	ey(1)=1;
	ey(2)=0;
	double dotChi_mui, dotmui_muj;
	double s0,s2,s3;
	Vector pp(3),p(3),p2(3);
	double di;
	double dij;
	
	dotChi_mui = dot(tessellationAxis,mu_i);
	
	s2 = 1;
	if(isWithinStrongNumericalLimits(dotChi_mui,1.0)){
		/*pp(0) = 0;
		pp(1) = 0;
		pp(2) = 1;
		
		p = normalise(cross(pp,mu_i));
		*/
		
		p(0) = 0;
		p(1) = 1;
		p(2) = 0;
		
		ni=normalise(cross(tessellationAxis,p),di);
		di=1;
		problem=true;
		
		if(dot(ni,mu_j)<0)
			s2 = -1;
		
	}
	else if(isWithinNumericalLimits(dotChi_mui,-1.0)){
		p(0) = 0;
		p(1) = -1;
		p(2) = 0;
		
		ni=normalise(cross(tessellationAxis,p),di);
		di=1;
		if(dot(ni,mu_j)<0)
			s2 = -1;
	}
	else{
		p = mu_i;
		ni = normalise(cross(tessellationAxis,p),di);
		
		
	}
		
	
	
	
	dotmui_muj = rhoContainer.dot_mui_muj;
	if(isWithinNumericalLimits(dotmui_muj,1.0)){
		p2(0) = 0;
		p2(1) = -1;
		p2(2) = 0;
		nij = normalise(cross(tessellationAxis,p2),dij);
		dij=1;
		//exit(-2);
	}
	else if(isWithinNumericalLimits(dotmui_muj,-1.0)){
		p2(0) = 0;
		p2(1) = 1;
		p2(2) = 0;
		nij = normalise(cross(tessellationAxis,p2),dij);
		dij=1;
		//exit(-3);
	}
	else{
		nij = normalise(cross(mu_i, mu_j),dij);
	}
	
	
	
	dot_ni_nij = dot(ni,nij);
	

	
	
	
	if(dot_ni_nij > 1.0) dot_ni_nij=1.0;
	if(dot_ni_nij < -1.0) dot_ni_nij=-1.0;
	
	varpi = acos(dot_ni_nij);
	
	s0 = -sgn(dot(nij,tessellationAxis));
	
	if(problem){
		s0=-1;
	}
	
	
	
	omega.rotation = s2 * s0 * varpi;
	
	
	if(isnan(omega.rotation)){
		fprintf(stderr,"omega evaluated to nan\n");
		exit(-2);
	}
	
	
	omega.di = di;
	omega.dij = dij;
	omega.ni = ni;
	omega.nij = nij;
	omega.dot_ni_nij = dot_ni_nij;
	omega.id_i = id_i;
	omega.id_j = id_j;
	omega.s0 = s0;
	omega.s2 = s2;
	
	return omega;
	
}
	
	
	
	
Rotation Tessellation::calculateOmegaDerivatives(OmegaRotation &r, CircularInterfacesPerAtom &circles, Vector &tessellationAxis){
	Matrix dmui_dxl, dmuj_dxl;
	double domega_dvarpi;
	Vector dvarpi_dnij, dvarpi_dni;
	Matrix dvi_dmui, dvij_dmui, dvij_dmuj;
	Matrix dni_dvi, dnij_dvij;
	Vector domega_dxi, domega_dxj, domega_dxl;

	double di,dij;
	Vector ni(3);
	Vector nij(3);
	Matrix Identity(3,3),reverseIdentity(3,3);
	Identity.eye();
	reverseIdentity.eye();
	reverseIdentity = reverseIdentity * -1;
	Rotation omega;
	double dot_ni_nij;
	Vector mu_i;
	Vector mu_j;
	Matrix dmui_dxi;
	Matrix dmuj_dxj;
	int id_i;
	int id_j;
	double s0,s2;
	double denominator;
	
	omega.rotation=r.rotation;
	
	di = r.di;
	dij = r.dij;
	ni = r.ni;
	nij = r.nij;
	dot_ni_nij = r.dot_ni_nij;
	
	id_i = r.id_i;
	id_j = r.id_j;
	
	s0=r.s0;
	s2=r.s2;
	
	
	mu_i = circles[id_i].normal;
	mu_j = circles[id_j].normal;
	dmui_dxi = circles[id_i].dmu_dx;
	dmuj_dxj = circles[id_j].dmu_dx;
	
	
	
	
	
	dmui_dxl = -dmui_dxi;
	dmuj_dxl = -dmuj_dxj;
	
	if(s2>0)
		dvi_dmui=-matrixCross(Identity,tessellationAxis);
	else
		dvi_dmui = Matrix(3,3).zeros();
	
	
	
	dvij_dmui = matrixCross(Identity,mu_j);
	dvij_dmuj = matrixCross(reverseIdentity,mu_i);
	
	dni_dvi = (1.0/(di))*(Identity - kron(ni,ni.t())).t();
	dnij_dvij = (1.0/(dij))*(Identity - kron(nij,nij.t())).t();
	
	
	domega_dvarpi = s0;
	
	
	
	if(isWithinStrongNumericalLimits(abs(dot_ni_nij),1.0)){
		
		dvarpi_dni = Vector(3).zeros();
		dvarpi_dnij = Vector(3).zeros();
		domega_dxi = Vector(3).zeros();
		domega_dxj = Vector(3).zeros();
		domega_dxl = Vector(3).zeros();
		
	
		
	}
	else{

		if(dot_ni_nij > 1.0-MINISCULE) dot_ni_nij = 1.0-MINISCULE;
		if(dot_ni_nij < -1.0+MINISCULE) dot_ni_nij = -1.0+MINISCULE;

	
		denominator = (sqrt(1-dot_ni_nij*dot_ni_nij));
		
	
		dvarpi_dni = -(nij) / denominator;
		dvarpi_dnij = -(ni) / denominator;
		
		
		//if(norm(dvarpi_dni,2) > 1.0) dvarpi_dni = normalise(dvarpi_dni);
		//if(norm(dvarpi_dnij,2) > 1.0) dvarpi_dnij = normalise(dvarpi_dnij);

		
		domega_dxi = s2* domega_dvarpi * ((dvarpi_dni.t() * dni_dvi * dvi_dmui * dmui_dxi).t() + (dvarpi_dnij.t() * dnij_dvij * dvij_dmui * dmui_dxi).t() );
		
		domega_dxj = s2* domega_dvarpi * (dvarpi_dnij.t() * dnij_dvij * dvij_dmuj * dmuj_dxj).t();

		domega_dxl = s2* domega_dvarpi * ( (dvarpi_dni.t() * dni_dvi * dvi_dmui * dmui_dxl).t() +  (dvarpi_dnij.t() * (dnij_dvij * dvij_dmui * dmui_dxl + dnij_dvij * dvij_dmuj * dmuj_dxl) ).t() );
	}
	
	omega.drotation_dxi = domega_dxi;
	omega.drotation_dxj = domega_dxj;
	omega.drotation_dxl = domega_dxl;
	


	
	return omega;
	
}




double Tessellation::calculateRho(Vector &mu_i, Vector &mu_j){
	return calculateRho(mu_i, mu_j, false);
}


double Tessellation::calculateRho(Vector &mu_i, Vector &mu_j, bool derivatives){
	double dot_IJ;
	double rho;
	dot_IJ = norm_dot(mu_i,mu_j);
	rho = acos(dot_IJ);
	return rho;
}


EtaRotation Tessellation::calculateEta(Vector &tessellationAxis, CircularInterface &I, CircularInterface &J, RhoContainer &rhoContainer){
	return calculateEta(tessellationAxis, I.normal, J.normal, I.lambda, J.lambda, I.id-1, J.id-1, I.form, J.form, rhoContainer);
}


EtaRotation Tessellation::calculateEta(Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, RhoContainer &rhoContainer){
	Vector tessellationAxis(3);
	
	tessellationAxis.zeros();
	return calculateEta(tessellationAxis, mu_i, mu_j, lambda_i, lambda_j, 0, 0, CONVEX, CONVEX, rhoContainer);
}


EtaRotation Tessellation::calculateEta(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer){
	double rho;
	double dot_IJ;
	double sig0,sig1,sig2;
	EtaRotation eta;
	
	
	
	
	dot_IJ = rhoContainer.dot_mui_muj;
	
	rho = rhoContainer.rho;
	
	
	
	if(isWithinNumericalLimits(rho,0)){
		eta.rotation = 0;
	}
	else{
	
		sig0 = cot(lambda_i.rotation) * cot(rho);
		sig1 = cos(lambda_j.rotation)*csc(lambda_i.rotation)*csc(rho);
		sig2 = sig0-sig1;
		
		
		eta.rotation = acos(sig2);
		
		eta.sig0 = sig0;
		eta.sig1 = sig1;
		eta.sig2 = sig2;
		eta.dot_IJ = dot_IJ;
		eta.rho = rho;

		
		
		
	}
	
	if(isnan(eta.rotation)){
		fprintf(stderr,"eta evaluated to nan\n");
		fprintf(stderr,"sig: %f %f %f\n",sig0,sig1,sig2);
		fprintf(stderr,"angles: %f %f %f\n",lambda_i.rotation, lambda_j.rotation, rho);
		exit(-2);
	}
	
	
	eta.id_i = id_i;
	eta.id_j = id_j;
	
	
	return eta;
}
	
	
	
Rotation Tessellation::calculateEtaDerivatives(EtaRotation &r, CircularInterfacesPerAtom &circles){
					       
	Vector drho_dmui(3), drho_dmuj(3);
	Vector drho_dxi(3), drho_dxj(3);
	Vector drho_dxl(3);
	double sig0,sig1,sig2;
	double deta_dlambdai, deta_dlambdaj, deta_drhoij;
	Vector deta_dxi(3), deta_dxj(3), deta_dxl(3);
	Rotation eta;
	Matrix dmui_dxl, dmuj_dxl;
	Vector mu_i;
	Vector mu_j;
	Rotation lambda_i;
	Rotation lambda_j;
	Matrix dmui_dxi;
	Matrix dmuj_dxj;
	int id_i;
	int id_j;
	double dot_IJ;
	double rho;
	
	sig0 = r.sig0;
	sig1 = r.sig1;
	sig2 = r.sig2;
	dot_IJ = r.dot_IJ;
	rho = r.rho;
	id_i = r.id_i;
	id_j = r.id_j;
	
	mu_i = circles[id_i].normal;
	mu_j = circles[id_j].normal;
	lambda_i = circles[id_i].lambda;
	lambda_j = circles[id_j].lambda;
	dmui_dxi = circles[id_i].dmu_dx;
	dmuj_dxj = circles[id_j].dmu_dx;

	eta.rotation = r.rotation;

	if(isWithinNumericalLimits(rho,0)){
		deta_dxi = Vector(3).zeros();
		deta_dxj = Vector(3).zeros();
		deta_dxl = Vector(3).zeros();
	}
	else{
		
		if(dot_IJ > 1.0-MINISCULE) dot_IJ = 1.0-MINISCULE;
		if(dot_IJ < -1.0+MINISCULE) dot_IJ = -1.0+MINISCULE;
		
	
		drho_dmui = -mu_j/(sqrt(1-dot_IJ*dot_IJ));
		drho_dmuj = -mu_i/(sqrt(1-dot_IJ*dot_IJ));
		
		dmui_dxl = -dmui_dxi;
		dmuj_dxl = -dmuj_dxj;
		
		
		drho_dxi = (drho_dmui.t() * dmui_dxi).t();
		drho_dxj = (drho_dmuj.t() * dmuj_dxj).t();
		
		drho_dxl = ((drho_dmui.t() * dmui_dxl).t() + (drho_dmuj.t() * dmuj_dxl).t());
		
		if(sig2 > 1.0-MINISCULE) sig2 = 1.0-MINISCULE;
		if(sig2 < -1.0+MINISCULE) sig2 = -1.0+MINISCULE;
		
		
		deta_dlambdai = csc(lambda_i.rotation) * (sig0/cos(lambda_i.rotation) - sig1*cos(lambda_i.rotation)) / sqrt(1-sig2*sig2);
		deta_dlambdaj =  -sin(lambda_j.rotation)*csc(lambda_i.rotation)*csc(rho) / sqrt(1-sig2*sig2);
		deta_drhoij = csc(rho) * (sig0/cos(rho) - sig1*cos(rho)) / sqrt(1-sig2*sig2);
		
		
		
		deta_dxi = deta_dlambdai * lambda_i.drotation_dxi + deta_drhoij * drho_dxi;
		deta_dxj = deta_dlambdaj * lambda_j.drotation_dxi + deta_drhoij * drho_dxj;
		deta_dxl = deta_dlambdai * lambda_i.drotation_dxl + deta_dlambdaj * lambda_j.drotation_dxl + deta_drhoij * drho_dxl;
		
		eta.drotation_dxi = deta_dxi;
		eta.drotation_dxj = deta_dxj;
		eta.drotation_dxl = deta_dxl;
	}
		


	
	
	return eta;
}

PHIContainer Tessellation::calculatePHI(Vector &tessellationAxis, CircularInterface &I, CircularInterface &J, double radius, RhoContainer &rhoContainer){
	return calculatePHI(tessellationAxis, I.normal, J.normal, I.lambda, J.lambda, I.id-1, J.id-1, I.form, J.form, rhoContainer);
	
}


PHIContainer Tessellation::calculatePHI(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, RhoContainer &rhoContainer){
	return calculatePHI(tessellationAxis, mu_i, mu_j, lambda_i, lambda_j, 0, 0, CONVEX, CONVEX, rhoContainer);
}


PHIContainer Tessellation::calculatePHI(Vector &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer){
	PHIContainer p,p2,p3;
	
	EtaRotation eta;
	OmegaRotation omega;
	int q;
	

	eta = calculateEta(tessellationAxis, mu_i, mu_j, lambda_i, lambda_j, id_i, id_j, form_i, form_j, rhoContainer);
	omega = calculateOmega(tessellationAxis, mu_i, mu_j, id_i, id_j, form_i, form_j, rhoContainer);
	
	
	
	p.out.rotation = omega.rotation + eta.rotation;
	if(p.out.rotation > M_PI)
		p.out.rotation = -2*M_PI + p.out.rotation;
	
	p.in.rotation = omega.rotation - eta.rotation;
	if(p.in.rotation < -M_PI)
		p.in.rotation = 2*M_PI + p.in.rotation;

	
	
	p.in.omega = omega;
	p.in.eta = eta;
	p.out.omega = omega;
	p.out.eta = eta;

	p.out.s = 1;
	p.in.s = -1;
	
	
	
	q = 1;
	if(form_i != CONVEX) q*=-1;
	if(form_j != CONVEX) q*=-1;
	
	
	
	
	if(q<0){
		p3=p;
		p.in=p3.out;
		p.out=p3.in;
	}

	
	


	
	return p;
	
	
	 
}


Rotation Tessellation::calculatePHIDerivatives(PHIRotation &r, CircularInterfacesPerAtom &circles, Vector &tessellationAxis){
	Rotation omega,eta;
	Rotation PHI;
	int s;
	
	omega = calculateOmegaDerivatives(r.omega, circles, tessellationAxis);
	eta = calculateEtaDerivatives(r.eta, circles);
	
	PHI.rotation = r.rotation;
	s=r.s;
		
	PHI.drotation_dxi = s*eta.drotation_dxi + omega.drotation_dxi;
	PHI.drotation_dxj = s*eta.drotation_dxj + omega.drotation_dxj;
	PHI.drotation_dxl = s*eta.drotation_dxl + omega.drotation_dxl;
	
	
	
	return PHI;
	
			
}









Rotation Tessellation::calculatePsi(Vector &tessellationAxis, CircularInterface &circle){
	return calculatePsi(tessellationAxis, circle.normal, circle.dmu_dx, circle.form, circle.id-1);
}

Rotation Tessellation::calculatePsi(Vector &tessellationAxis, Vector &mu_i){
	Matrix dmu_dx(3,3);
	CircularInterfaceForm form;
	int id;
	
	dmu_dx.zeros();
	form = CONVEX;
	id=0;
	
	return calculatePsi(tessellationAxis, mu_i, dmu_dx, form, id);
}

Rotation Tessellation::calculatePsi(Vector &tessellationAxis, Vector &mu_i, Matrix &dmu_dx, CircularInterfaceForm form, int index){
	Vector dpsi_dmui;
	Rotation r;
	Vector dpsi_dxi, dpsi_dxl;
	
	Vector n_i, nn_i;
	Vector ey(3);
	Matrix I(3,3);
	I.eye();
	
	Matrix dni_dmui;
	Matrix dnni_dni;
	Vector dpsi_dnni;

	Matrix dni_dxi;
	Matrix dnni_dxi;
	double dot_mui_chi;
	
	
	r.drotation_dxi = Vector(3).zeros();
	r.drotation_dxj = Vector(3).zeros();
	r.drotation_dxl = Vector(3).zeros();
	
	r.rotation = getAngle(tessellationAxis,mu_i);
	

	if(form != SPLITTER){
		
		
		dot_mui_chi = dot(mu_i,tessellationAxis);

			
		if(isWithinNumericalLimits(dot_mui_chi,1.0) || isWithinNumericalLimits(dot_mui_chi,-1.0)){
			dpsi_dmui = Vector(3).zeros();
			dpsi_dxi = Vector(3).zeros();
			dpsi_dxl = Vector(3).zeros();
		}
		else{
			
			//if(mu_i(0) > 1.0-MINISCULE) mu_i(0) = 1.0-MINISCULE;
			//if(mu_i(0) < -1.0+MINISCULE) mu_i(0) = -1.0+MINISCULE;
			
			if(dot_mui_chi > 1.0-MINISCULE) dot_mui_chi = 1.0-MINISCULE;
			if(dot_mui_chi < -1.0+MINISCULE) dot_mui_chi = -1.0+MINISCULE;
			
			dpsi_dmui = -tessellationAxis/(sqrt(1-dot_mui_chi*dot_mui_chi));
			
			
			dpsi_dxi = (dpsi_dmui.t() * dmu_dx).t();
			//smoothing
			//if(norm(dpsi_dxi,2)>1.0) dpsi_dxi = normalise(dpsi_dxi);
			
			dpsi_dxl = -dpsi_dxi;
		}
		
		
		r.drotation_dxi = dpsi_dxi;
		r.drotation_dxj = Vector(3).zeros();
		r.drotation_dxl = dpsi_dxl;
	}
	else{
		//the splitter is immovable
		//r.drotation_dxi = Vector(3).zeros();
		//r.drotation_dxj = Vector(3).zeros();
		//r.drotation_dxl = Vector(3).zeros();
	}
	

	
	
	
	return r;
	
}




IntersectionBranches::iterator Tessellation::increaseBranchInterator(IntersectionBranches::iterator it, CircularInterface &circle){
	IntersectionBranches* p;
	p = it->second.body;
	if(circle.form == CONVEX){
		++it;
		if(it == p->end()) it=p->begin();
		
	}
	else{
		if(it == p->begin()) it=p->end();
		--it;
	}
	
	return it;
	
}

IntersectionBranches::iterator Tessellation::decreaseBranchInterator(IntersectionBranches::iterator it, CircularInterface &circle){
	IntersectionBranches* p;
	p = it->second.body;
	
	if(circle.form == CONVEX){
		if(it == p->begin()) it=p->end();
		--it;
		
	}
	else{
		++it;
		if(it == p->end()) it=p->begin();
	}
	
	return it;
}






void Tessellation::createIntersectionBranch(PHIContainer &PHII, PHIContainer &PHIJ, CircularInterface &I, CircularInterface &J){
	
	pair<double, IntersectionBranch> x;
	IntersectionBranches::iterator it0;
	IntersectionBranches::iterator it1;
	
	x.second.visited = -1;
	
	x.first = PHII.in.rotation;
	x.second.direction = IN;
	x.second.body = &I.intersectionBranches;
	x.second.id = J.id;
	x.second.flagged = false;
	x.second.PHI = PHII.in;
	it0 = I.intersectionBranches.insert(x);
	
	
	x.first = PHII.out.rotation;
	x.second.direction = OUT;
	x.second.body = &I.intersectionBranches;
	x.second.id = J.id;
	x.second.flagged = false;
	x.second.PHI = PHII.out;
	it1 = I.intersectionBranches.insert(x);

	
		
	
	
}








void Tessellation::printBranch(const char* s, multimap<double, IntersectionBranch>::iterator &it){
	if(it->second.direction == OUT)
		printf("iterator %s: %d-%d (OUT)\n",s,it->second.node->id0, it->second.node->id1);
	else
		printf("iterator %s: %d-%d (IN)\n",s,it->second.node->id0, it->second.node->id1);
}

void Tessellation::printIntersectionGraph(IntersectionGraph &g,CircularInterfacesPerAtom &circles){
	IntersectionGraph::iterator it;
	
	for(it = g.begin(); it != g.end(); ++it){
		//fprintf(stderr,"NODE[%d,%d]: (%d,%d) -> (%d,%d)\n", circles[it->first.id0-1].index, circles[it->first.id1-1].index, circles[it->second.id0-1].index, circles[it->second.id1-1].index, circles[it->second.pointsTo0-1].index, circles[it->second.pointsTo1-1].index);
		fprintf(stderr,"NODE[%d,%d]: (%d,%d) -> (%d,%d)\n", it->first.id0, it->first.id1, it->second.id0, it->second.id1, it->second.pointsTo0, it->second.pointsTo1);
	}
}









void Tessellation::buildIntersectionGraph(double radius, Vector &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, string filename){
	SASASegment sasaSegment;
	
	unsigned int j;
	CircularInterface *I, *J;
	vector<RhoContainer>::iterator it_j;
	
	
	Interfaces interfacesJ, interfacesI;
	IntersectionBranches::iterator it, it2;
	PHIContainer PHIJ, PHII;
	
	int searchMode;
	bool erased;
	bool done;
	RhoContainer rhoContainer;
	
	
	//iterate through all circles and add them to the intersectiongraph, one by one
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		
		
		if(I->circularIntersections.size()==0){
			insertArtificialIntersectionPoints(*I,tessellationAxis,hemisphere,sasa);			
		}
		
		
		
		for(it_j=I->circularIntersections.begin(); it_j != I->circularIntersections.end(); ++it_j){
			j=it_j->id;
			rhoContainer = *it_j;
			J = &circles[j];
			
			//retrieve external and internal interfaces (respectively)
			PHIJ = calculatePHI(tessellationAxis, *J, *I, radius, rhoContainer);
			PHII = calculatePHI(tessellationAxis, *I, *J, radius, rhoContainer);
			
			createIntersectionBranch(PHII, PHIJ, *I, *J);
			
			

			
			
			
			for(it = I->intersectionBranches.begin(); it != I->intersectionBranches.end(); ++it){
				if(it->second.id == J->id) searchMode=0;
				else searchMode=1;
				
				if(it->second.direction==OUT){
				
				
					done=false;
					it2 = it;
					while(!done){
						it2=increaseBranchInterator(it2, *I);
						
						
						
						if(	(it2->second.direction==IN && searchMode==0 && it2->second.id == J->id) ||
							(it2->second.direction==IN && searchMode==1 && it2->second.id != J->id)){
								done=true;
							}
							else{
								it2->second.flagged=true;
							}
							
							
					}
				}
			}
			
			it = I->intersectionBranches.begin();
			while(it != I->intersectionBranches.end()){
				erased=false;
				if(it->second.flagged){
					it2=it;
					++it2;
					I->intersectionBranches.erase(it);
					it=it2;
					erased=true;
				}
				
				if(!erased) ++it;
					
			}
			
			
			if(I->intersectionBranches.size()==0)
				break;
			
			
			
					
		}
	}
	
	
	
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		
		for(it = I->intersectionBranches.begin(); it != I->intersectionBranches.end(); ++it){
			if(it->second.direction==IN){
				it2=increaseBranchInterator(it,*I);
				
				
				sasaSegment.id0 = it->second.id;
				sasaSegment.id1 = I->id;
				sasaSegment.id2 = it2->second.id;
				
				sasaSegment.index0 = circles[it->second.id-1].index;
				sasaSegment.index1 = circles[I->id-1].index;
				sasaSegment.index2 = circles[it2->second.id-1].index;
				
				if(!circles[it->second.id-1].hasDerivatives) calculateProjectionAndDerivatives(tessellationAxis, circles[it->second.id-1]);
				if(!circles[I->id-1].hasDerivatives) calculateProjectionAndDerivatives(tessellationAxis, circles[I->id-1]);
				if(!circles[it2->second.id-1].hasDerivatives) calculateProjectionAndDerivatives(tessellationAxis, circles[it2->second.id-1]);
				
				sasaSegment.lambda = circles[I->id-1].lambda;
				sasaSegment.psi = circles[I->id-1].psi;
				
				
				
				sasaSegment.rotation0 = calculatePHIDerivatives(it->second.PHI, circles, tessellationAxis);
				sasaSegment.rotation1 = calculatePHIDerivatives(it2->second.PHI, circles, tessellationAxis);
				sasaSegment.normalForCircularInterface = circles[I->id-1].normal;
				sasaSegment.form0 = circles[it->second.id-1].form;
				sasaSegment.form1 = circles[I->id-1].form;
				sasaSegment.form2 = circles[it2->second.id-1].form;
				
				sasaSegment.tessellationAxis = tessellationAxis;
				sasaSegment.hemisphere = hemisphere;
				
			
				
				
				sasa.push_back(sasaSegment);
				
			}
				
				
				
				
				
				
		}
		
		
	}
	





	
	
}
















void Tessellation::outputGaussBonnetData(string filename, double radius, CircularInterfacesPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph){
/*	FILE* file;
	
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
			//fprintf(file, "intersector %d radius %f vector %f %f %f\n", i, circle.sphereRadius, circle.vector(0), circle.vector(1), circle.vector(2));
		}
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
	
	*/	
	

}
