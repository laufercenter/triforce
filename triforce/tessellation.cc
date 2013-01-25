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
		//int i=12;{
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
	C.normal = v;
	C.form = form;
	C.id = circles.size()+1;
	C.a=radius;
	C.valid = true;
	C.index = -1;
	
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
	torigin = origin;
	radius = radii[i];
	tradius = radius;
	
	
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


bool Tessellation::isWithinStrongNumericalLimits(double x, double l){
	if(abs(x-l)<=THRESHOLD_STRONG_NUMERICAL) return true;
	else return false;
}


Rotation Tessellation::calculateLambda(double d_i, double r_l, double r_i, Vector &mu_i){
	double g2;
	CircularInterfaceForm form;
	return calculateLambda(0,d_i, r_l, r_i, mu_i, form, g2, false);
}


Rotation Tessellation::calculateLambda(int index_i, double d_i, double r_l, double r_i, Vector &mu_i, CircularInterfaceForm &form, double &g2, bool derivatives){
	double g;
	double lambda;
	double dg_dd;
	double dlambda_dg;
	Rotation r;
	double q;
	Vector mu_i_original;
	Vector dlambda_dxi(3), dlambda_dxl(3);
	Vector dd_dxi;
	
	mu_i_original = mu_i;
	g = (d_i * d_i + r_l * r_l - r_i * r_i ) / (2 * d_i * r_l);

	g2 = (d_i * d_i + r_l * r_l - r_i * r_i ) / (2 * d_i);

	if(g<0){
		form=CONCAVE;
		g = abs(g);
		mu_i=-mu_i;
	}
	else form=CONVEX;
	
	r.rotation = acos(g);
	
	
	if(derivatives){
	
	
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
		dlambda_dxl = -dlambda_dxi;
		
		r.drotation_dxi = dlambda_dxi;
		r.drotation_dxj = Vector(3).zeros();
		r.drotation_dxl = dlambda_dxl;
		
		
		
		
		
		//finite differences
		double fd_dlambda_dg;
		double fd_dg_dd;
		Vector fd_dmu(3);
		Vector x(3), xp(3), xn(3);
		Vector np(3), nn(3);
		Vector vp(3), vn(3);
		double lenvn, lenvp;
		double dp, dn;
		Rotation lambdap, lambdan;
		Vector fd_dlambda_dxi(3);
		Vector fd_dd_dxi(3);
		
		//dlambda_dg
		fd_dlambda_dg = (acos(g+FD)-acos(g-FD))/(2*FD);
		printf("dlambda_dg: %f FD: %f\n",dlambda_dg, fd_dlambda_dg);
		if(abs(dlambda_dg - fd_dlambda_dg) > FDT) exit(-1);
		
		//dg_dd
		if(form==CONVEX){
			fd_dg_dd = (((d_i+FD) * (d_i+FD) + r_l * r_l - r_i * r_i ) / (2 * (d_i+FD) * r_l) - ((d_i-FD) * (d_i-FD) + r_l * r_l - r_i * r_i ) / (2 * (d_i-FD) * r_l))/(2*FD);
		}
		else{
			fd_dg_dd = (-((d_i+FD) * (d_i+FD) + r_l * r_l - r_i * r_i ) / (2 * (d_i+FD) * r_l) - -((d_i-FD) * (d_i-FD) + r_l * r_l - r_i * r_i ) / (2 * (d_i-FD) * r_l))/(2*FD);
		}
		printf("dg_dd: %f FD: %f\n",dg_dd,fd_dg_dd);
		if(abs(dg_dd - fd_dg_dd) > FDT) exit(-1);
		
		//dd_dxi
		x = atoms[index_i];
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			np = calculateInterfaceNormal(torigin, xp, dp);
			nn = calculateInterfaceNormal(torigin, xn, dn);
			
			fd_dd_dxi(j) = (dp - dn) / (2*FD);
		}
		printf("dd_dxi: %f %f %f FD: %f %f %f\n",dd_dxi(0),dd_dxi(1),dd_dxi(2), fd_dd_dxi(0), fd_dd_dxi(1), fd_dd_dxi(2));
		for(int i=0; i<3; ++i){
			if(abs(dd_dxi(i) - fd_dd_dxi(i)) > FDT) exit(-1);
		}
			
		
		
		
		//dlambda_dxi
		x = atoms[index_i];
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;

			np = calculateInterfaceNormal(torigin, xp, dp);
			nn = calculateInterfaceNormal(torigin, xn, dn);
			lambdap = calculateLambda(dp, tradius, r_i, np);
			lambdan = calculateLambda(dn, tradius, r_i, nn);
			
			fd_dlambda_dxi(j) = (lambdap.rotation - lambdan.rotation) / (2*FD);
		}
		printf("dlambda_dxi: %f %f %f FD: %f %f %f\n",dlambda_dxi(0),dlambda_dxi(1),dlambda_dxi(2), fd_dlambda_dxi(0), fd_dlambda_dxi(1), fd_dlambda_dxi(2));
		for(int i=0; i<3; ++i){
			if(abs(dlambda_dxi(i) - fd_dlambda_dxi(i)) > FDT) exit(-1);
		}
		
			
		
		
		
	}
	
	
	
	
	return r;
	
}


void Tessellation::determineProjection(Vector &origin, double r_l, CircularInterface &circle){
	double d_i;
	double r_i;
	double g, g_normalised;
	double a;
	Matrix Identity(3,3);
	double dlambda_dg;
	double dg_dd;
	Matrix dmu_dx(3,3);
	Matrix fd_dmu_dx(3,3);
	
	
	Identity.eye();
	
	
	
	

	circle.lambda = calculateLambda(circle.index, circle.d, r_l, circle.sphereRadius, circle.normal, circle.form, circle.g, true);
	circle.a= sqrt(r_l * r_l - circle.g * circle.g);
	circle.valid=true;

	
	//derivatives
	if(circle.form==CONVEX)
		dmu_dx = (1.0/circle.d)*(Identity - kron(circle.normal,circle.normal.t())).t();
	else
		dmu_dx = -(1.0/circle.d)*(Identity - kron(circle.normal,circle.normal.t())).t();
	
	circle.dmu_dx = dmu_dx;
	
	
	
	
	
	
	
	
	Vector fd_dmu(3);
	Vector x(3), xp(3), xn(3);
	Vector np(3), nn(3);
	Vector vp(3), vn(3);
	double lenvn, lenvp;
	int i,j;
	
	
	
	
	//dmu_dx
	x = atoms[circle.index];
	
	for(j=0; j<3; ++j){
		xp=x;
		xp(j) = x(j) + FD;
		xn=x;
		xn(j) = x(j) - FD;
		
		
		vp= xp - origin;
		lenvp = norm(vp,2);
		np=vp/lenvp;
		if(circle.form!=CONVEX)
			np = -np;

		vn= xn - origin;
		lenvn = norm(vn,2);
		nn=vn/lenvn;
		if(circle.form!=CONVEX)
			nn = -nn;
		
		fd_dmu = (np - nn)/(2*FD);
		
		for(i=0; i<3; ++i){
			fd_dmu_dx(i,j) = fd_dmu(i);
		}
	}
	
	
	for(j=0; j<3; ++j){
		printf("dmu_dx: %f %f %f\n",dmu_dx(j,0),dmu_dx(j,1),dmu_dx(j,2));
	}
	for(j=0; j<3; ++j){
		printf("fd_dmu_dx: %f %f %f\n",fd_dmu_dx(j,0),fd_dmu_dx(j,1),fd_dmu_dx(j,2));
	}
	
	for(i=0; i<3; ++i)
		for(j=0; j<3; ++j)
			if(abs(dmu_dx(i,j)-fd_dmu_dx(i,j)) > FDT) exit(-1);
	
	
	


	
	

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


Vector Tessellation::calculateInterfaceNormal(const Vector &v_l, const Vector &v_i){
	double d;
	return calculateInterfaceNormal(v_l,v_i);
}


Vector Tessellation::calculateInterfaceNormal(const Vector &v_l, const Vector &v_i, double &d){
	Vector v;
	v = v_i - v_l;
	d = norm(v,2);
	v = v/d;
	return v;
	
}


void Tessellation::makeCircularInterfaces(int l,Vector &origin, double r_l, vector<vec> &atoms, vector<double> &radii, vector<CircularInterface> &circles){
	CircularInterface circle;
	double r_i;
	double d_i;
	Vector mu_i;
	
	printf("MAKECIRCLE %d %f\n",l, r_l);

	for(int j=0;j<atoms.size();j++){
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
				circles.push_back(circle);
				
				printf("ADDCIRCLE %d %f %f\n",j,r_i,d_i);
				
				
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
				
				angle = getAngle(n0,n1);
				printf("FILTER: %d %d: %f %f %f\n",it->index, circles[i].index, it->lambda.rotation, circles[i].lambda.rotation, angle);
				
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



Rotation Tessellation::calculateOmega(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J){
	return calculateOmega(tessellationOrigin, I.normal, J.normal, I.dmu_dx, J.dmu_dx, I.index, J.index, I.form, J.form, true);

}


Rotation Tessellation::calculateOmega(Vector &tessellationOrigin, Vector &mu_i, Vector &mu_j){
	Matrix dmu_dx(3,3);
	dmu_dx.zeros();
	
	return calculateOmega(tessellationOrigin, mu_i, mu_j, dmu_dx, dmu_dx, 0, 0, CONVEX, CONVEX, false);
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


Rotation Tessellation::calculateOmega(Vector &tessellationOrigin, Vector &mu_i, Vector &mu_j, Matrix &dmui_dxi, Matrix &dmuj_dxj, int index_i, int index_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, bool derivatives){
	Vector ex(3);
	Vector ez(3);
	Vector ey(3);
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
	double dot_ni_nij, dot_ni_nij2;
	Matrix dmui_dxl, dmuj_dxl;
	double domega_dvarpi;
	Vector dvarpi_dnij, dvarpi_dni;
	Matrix dvi_dmui, dvij_dmui, dvij_dmuj;
	Matrix dni_dvi, dnij_dvij;
	Vector domega_dxi, domega_dxj, domega_dxl;
	
	ez(0)=0;
	ez(1)=0;
	ez(2)=1;
	ey(0)=0;
	ey(1)=1;
	ey(2)=0;
	double dotChi_mui, dotmui_muj, d;
	double s0,s1,s2;
	Vector p(3),p2(3);
	double di;
	double dij;
	Vector vi,vij;
	
	dotChi_mui = dot(tessellationOrigin,mu_i);
	s2 = 1;
	if(isWithinNumericalLimits(dotChi_mui,1.0)){
		//printf("POSITIVE RANGE +1 %f\n",dotChi_mui);
		p(0) = 0;
		p(1) = 1;
		p(2) = 0;
		
		ni=normalise(cross(tessellationOrigin,p),di);
		vi=cross(tessellationOrigin,p);
		
		if(dot(ni,mu_j)<0)
			s2 = -1;
	}
	else if(isWithinNumericalLimits(dotChi_mui,-1.0)){
		//printf("POSITIVE RANGE -1 %f\n",dotChi_mui);
		p(0) = 0;
		p(1) = -1;
		p(2) = 0;
		
		ni=normalise(cross(tessellationOrigin,p),di);
		vi=cross(tessellationOrigin,p);
		if(dot(ni,mu_j)<0)
			s2 = -1;
	}
	else{
		p = mu_i;
		ni = normalise(cross(tessellationOrigin,p),di);
		vi=cross(tessellationOrigin,p);
	}
		
	
	
	//printf("CROSS: %f,%f,%f\n",ni(0),ni(1),ni(2));
	//printf("P: %f %f %f (%f)\n",p(0),p(1),p(2),dotChi_mui);
	//printf("TESSELLATION: %f,%f,%f\n",tessellationOrigin(0),tessellationOrigin(1),tessellationOrigin(2));
	
	//nii = cross(mu_i,ni);
	
	dotmui_muj = dot(mu_i, mu_j);
	if(isWithinNumericalLimits(dotmui_muj,1.0)){
		p2(0) = 0;
		p2(1) = -1;
		p2(2) = 0;
		nij = normalise(cross(tessellationOrigin,p2),dij);
		vij=cross(tessellationOrigin,p2);
		//exit(-2);
	}
	else if(isWithinNumericalLimits(dotmui_muj,-1.0)){
		p2(0) = 0;
		p2(1) = 1;
		p2(2) = 0;
		nij = normalise(cross(tessellationOrigin,p2),dij);
		vij=cross(tessellationOrigin,p2);
		//exit(-3);
	}
	else{
		nij = normalise(cross(mu_i, mu_j),dij);
		vij=cross(mu_i, mu_j);
	}
	
	
	
	dot_ni_nij = dot((ni),(nij));
	double dot_ni_nijtest = dot(ni,nij);
	
	dot_ni_nij2 = sdot(ni,nij);

	
	if(!isWithinNumericalLimits(dot_ni_nij,dot_ni_nijtest)){
		printf("DOTFAIL: %f %f\n",dot_ni_nij, dot_ni_nijtest);
	}
	
	
	double dninij=norm_dot((ni),(nij));
//	if(dninij > 1.0) dninij=2.0 - dninij;
//	if(dninij < -1.0) dninij=-2.0 - dninij;
	if(dninij > 1.0) dninij=1.0;
	if(dninij < -1.0) dninij=-1.0;
	
		varpi = acos(dninij);
	
	s0 = -sgn(dot(nij,tessellationOrigin));
	
	//omega.rotation = -s0*(s1*(1-s0)*M_PI/2 - varpi);
	
	omega.rotation = s0*varpi;
	/*
	if(omega.rotation < -M_PI){
		printf("omega.rotation < -M_PI\n");
		omega.rotation = 2*M_PI + omega.rotation;
	}

	if(omega.rotation > M_PI){
		printf("omega.rotation > M_PI\n");
		omega.rotation = -2*M_PI + omega.rotation;
	}
	*/
	omega.rotation = omega.rotation * s2;
	
	
	//printf("s0: %f, s2: %f\n",s0,s2);
	
	//printf("LENGTH mui: %f muj: %f, ni: %f, nij: %f\n",norm(mu_i,2), norm(mu_j,2), norm(ni,2), norm(nij,2));
	//printf("LENGTH2 mui: %f muj: %f, ni: %f, nij: %f\n",l(mu_i), l(mu_j), l(ni), l(nij));
	
	if(isnan(omega.rotation)) exit(-2);
	
	
	if(derivatives){
	
		//if(dot_ni_nij2 > 1.0 - MINISCULE) dot_ni_nij2 = 1.0-MINISCULE;
		//if(dot_ni_nij2 < -1.0 + MINISCULE) dot_ni_nij2 = -1.0+MINISCULE;
		
		
	printf("dotninij: %f %f\n",dot_ni_nij, dot_ni_nij2);
	
	printf("ni: %f %f %f nij: %f %f %f, dot: %f\n",ni(0),ni(1),ni(2),nij(0),nij(1),nij(2),dot_ni_nij);
		

		
		dmui_dxl = -dmui_dxi;
		dmuj_dxl = -dmuj_dxj;
		
		if(s2>0)
			dvi_dmui=-matrixCross(Identity,tessellationOrigin);
		else
			dvi_dmui = Matrix(3,3).zeros();
		
		
		
		dvij_dmui = matrixCross(Identity,mu_j);
		dvij_dmuj = matrixCross(reverseIdentity,mu_i);
		
		dni_dvi = (1.0/(di))*(Identity - kron(ni,ni.t())).t();
		dnij_dvij = (1.0/(dij))*(Identity - kron(nij,nij.t())).t();
		
		
		domega_dvarpi = s0;
		
		
		
		if(isWithinStrongNumericalLimits(abs(dot_ni_nij),1.0)){
			printf("ZERO LEVEL\n");
			
			/*
			dvarpi_dni = -(nij);
			dvarpi_dnij = -(ni);
			
			double s3;
			s3 = sgn(dot_ni_nij2);
			
		
			Vector a2 = (ni.t() * dni_dmui * dmui_dxi).t();
			Vector b2 = (nij.t() * dnij_dmui * dmui_dxi).t();
			
			printf("A2: (%f %f %f) B2: (%f %f %f)\n",a2(0),a2(1),a2(2),b2(0),b2(1),b2(2));
			
			
			domega_dxi = s2* domega_dvarpi * ((ni.t() * dni_dmui * dmui_dxi).t() + (nij.t() * dnij_dmui * dmui_dxi).t() );
			
			if(dot(a2,b2)>0){
				domega_dxi *=-s3*M_PI/2.0;
			}
			else{
				domega_dxi *=-s3*M_PI;
			}
			
			domega_dxj = s2* domega_dvarpi * (nij.t() * dnij_dmuj * dmuj_dxj).t();
			domega_dxj *=-s3*M_PI/2.0;
			
			
		Vector a = (ni.t() * dni_dmui * dmui_dxl).t();
		Vector b = (nij.t() * (dnij_dmui * dmui_dxl + dnij_dmuj * dmuj_dxl) ).t();
		
		printf("A: (%f %f %f) B: (%f %f %f)\n",a(0),a(1),a(2),b(0),b(1),b(2));
			
			
			domega_dxl = s2* domega_dvarpi * (  (ni.t() * dni_dmui * dmui_dxl).t() +  (nij.t() * (dnij_dmui * dmui_dxl + dnij_dmuj * dmuj_dxl) ).t() );
			domega_dxl *= -s3*M_PI/2.0;
			
			*/
			
			dvarpi_dni = Vector(3).zeros();
			dvarpi_dnij = Vector(3).zeros();
			domega_dxi = Vector(3).zeros();
			domega_dxj = Vector(3).zeros();
			domega_dxl = Vector(3).zeros();
			
		
			
		}
		else{
		
		
		
			dvarpi_dni = -(nij) / (sqrt(1-dot_ni_nij*dot_ni_nij));
			dvarpi_dnij = -(ni) / (sqrt(1-dot_ni_nij*dot_ni_nij));
			
			

			
			domega_dxi = s2* domega_dvarpi * ((dvarpi_dni.t() * dni_dvi * dvi_dmui * dmui_dxi).t() + (dvarpi_dnij.t() * dnij_dvij * dvij_dmui * dmui_dxi).t() );

			
			domega_dxj = s2* domega_dvarpi * (dvarpi_dnij.t() * dnij_dvij * dvij_dmuj * dmuj_dxj).t();

			domega_dxl = s2* domega_dvarpi * ( (dvarpi_dni.t() * dni_dvi * dvi_dmui * dmui_dxl).t() +  (dvarpi_dnij.t() * (dnij_dvij * dvij_dmui * dmui_dxl + dnij_dvij * dvij_dmuj * dmuj_dxl) ).t() );
		}
		
		omega.drotation_dxi = domega_dxi;
		omega.drotation_dxj = domega_dxj;
		omega.drotation_dxl = domega_dxl;
		
		
		
		//finite differences
		Vector fd_dvarpi_dni(3);
		Vector fd_dvarpi_dnij(3);
		Vector fd_dvi_dmui_row(3);
		Vector fd_dvij_dmui_row(3);
		Vector fd_dvij_dmuj_row(3);
		Matrix fd_dvi_dmui(3,3);
		Matrix fd_dvij_dmui(3,3);
		Matrix fd_dvij_dmuj(3,3);
		Vector fd_domega_dxi(3);
		Vector fd_domega_dxj(3);
		Vector fd_domega_dxl(3);
		Vector fd_dni_dvi_row(3);
		Vector fd_dnij_dvij_row(3);
		Matrix fd_dni_dvi(3,3);
		Matrix fd_dnij_dvij(3,3);
		
		
		Vector x(3), xp(3), xn(3);
		Vector np(3), nn(3);
		Vector np2(3), nn2(3);
		Vector np3(3), nn3(3);
		Vector vp(3), vn(3);
		double lenvn, lenvp;
		Vector nip,nin, nijp, nijn;
		double omp,omn;
		
		//dvarpi_dni
		x = normalise(ni);
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			if(abs(dot(xp, normalise(nij))) > 1.0 || abs(dot(xn,normalise(nij))) > 1.0){
				xp=x;
				xp(j) = x(j) + 2 * FD;
				xn=x;
			}
			if(abs(dot(xp, normalise(nij))) > 1.0 || abs(dot(xn,normalise(nij))) > 1.0){
				xp=x;
				xn=x;
				xn(j) = x(j) - 2 * FD;
			}
			
			

			fd_dvarpi_dni(j) = (acos(dot((xp),normalise(nij))) - acos(dot((xn),normalise(nij))))/(2*FD);
			
		}		
		printf("dvarpi_dni: (%f %f %f) fd_dvarpi_dni: (%f %f %f)\n",dvarpi_dni(0),dvarpi_dni(1),dvarpi_dni(2),fd_dvarpi_dni(0),fd_dvarpi_dni(1),fd_dvarpi_dni(2));
		for(int j=0; j<3; ++j){
			//if(abs(dvarpi_dni(j) - fd_dvarpi_dni(j)) > FDT || isnan(dvarpi_dni(j)) || isnan(fd_dvarpi_dni(j))) exit(-1);
			//if(isnan(dvarpi_dni(j)) ) exit(-1);
		}
		
		
		//dvarpi_dnij
		x = normalise(nij);
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			if(abs(dot(normalise(ni),xp)) > 1.0 || abs(dot(normalise(ni),xn)) > 1.0){
				xp=x;
				xp(j) = x(j) + 2 * FD;
				xn=x;
			}
			if(abs(dot(normalise(ni),xp)) > 1.0 || abs(dot(normalise(ni),xn)) > 1.0){
				xp=x;
				xn=x;
				xn(j) = x(j) - 2 * FD;
			}
			
			

			fd_dvarpi_dnij(j) = (acos(dot(normalise(ni),xp)) - acos(dot(normalise(ni),xn)))/(2*FD);
		}		
		printf("dvarpi_dnij: (%f %f %f) fd_dvarpi_dnij: (%f %f %f)\n",dvarpi_dnij(0),dvarpi_dnij(1),dvarpi_dnij(2),fd_dvarpi_dnij(0),fd_dvarpi_dnij(1),fd_dvarpi_dnij(2));
		for(int j=0; j<3; ++j){
			//if(abs(dvarpi_dnij(j) - fd_dvarpi_dnij(j)) > FDT || isnan(dvarpi_dnij(j)) || isnan(fd_dvarpi_dnij(j))) exit(-1);
			//if(isnan(dvarpi_dnij(j)) ) exit(-1);
		}
		
		
		//dvi_dmui
		if(s2>0){
			x = normalise(mu_i);
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				

				fd_dvi_dmui_row = (cross(tessellationOrigin,xp) - cross(tessellationOrigin,xn))/(2*FD);
				for(int i=0; i<3; ++i){
					fd_dvi_dmui(i,j) = fd_dvi_dmui_row(i);
				}
				
			}		
			for(int j=0; j<3; ++j){
				printf("dvi_dmui: %f %f %f \t fd_dvi_dmui: %f %f %f\n",dvi_dmui(j,0),dvi_dmui(j,1),dvi_dmui(j,2),fd_dvi_dmui(j,0),fd_dvi_dmui(j,1),fd_dvi_dmui(j,2));
			}
			
			for(int i=0; i<3; ++i)
				for(int j=0; j<3; ++j){
					if(abs(dvi_dmui(i,j)-fd_dvi_dmui(i,j)) > FDT) exit(-1);
				}
		}
		
		
		
		//dvij_dmui
		x = normalise(mu_i);
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;

			fd_dvij_dmui_row = ((cross(xp, mu_j)) - (cross(xn, mu_j)))/(2*FD);
			for(int i=0; i<3; ++i){
				fd_dvij_dmui(i,j) = fd_dvij_dmui_row(i);
			}
			
		}		
		for(int j=0; j<3; ++j){
			printf("dvij_dmui: %f %f %f \t fd_dvij_dmui: %f %f %f\n",dvij_dmui(j,0),dvij_dmui(j,1),dvij_dmui(j,2),fd_dvij_dmui(j,0),fd_dvij_dmui(j,1),fd_dvij_dmui(j,2));
		}
		
		for(int i=0; i<3; ++i)
			for(int j=0; j<3; ++j){
				if(abs(dvij_dmui(i,j)-fd_dvij_dmui(i,j)) > FDT) exit(-1);
			}
		
			
		//dvij_dmuj
		x = mu_j;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;

			fd_dvij_dmuj_row = (cross(mu_i, xp) - cross(mu_i, xn))/(2*FD);
			for(int i=0; i<3; ++i){
				fd_dvij_dmuj(i,j) = fd_dvij_dmuj_row(i);
			}
			
		}		
		for(int j=0; j<3; ++j){
			printf("dvij_dmuj: %f %f %f \t fd_dvij_dmuj: %f %f %f\n",dvij_dmuj(j,0),dvij_dmuj(j,1),dvij_dmuj(j,2),fd_dvij_dmuj(j,0),fd_dvij_dmuj(j,1),fd_dvij_dmuj(j,2));
		}
		
		for(int i=0; i<3; ++i)
			for(int j=0; j<3; ++j){
				if(abs(dvij_dmuj(i,j)-fd_dvij_dmuj(i,j)) > FDT) exit(-1);
			}
			
			
			
			
		//dni_dvi
		x = vi;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			

			fd_dni_dvi_row = (normalise(xp) - normalise(xn))/(2*FD);
			
			for(int i=0; i<3; ++i){
				fd_dni_dvi(i,j) = fd_dni_dvi_row(i);
			}
			
		}		
		for(int j=0; j<3; ++j){
			printf("dni_dvi: %f %f %f \t fd_dni_dvi: %f %f %f\n",dni_dvi(j,0),dni_dvi(j,1),dni_dvi(j,2),fd_dni_dvi(j,0),fd_dni_dvi(j,1),fd_dni_dvi(j,2));
		}
		
		for(int i=0; i<3; ++i)
			for(int j=0; j<3; ++j){
				if(abs(dni_dvi(i,j)-fd_dni_dvi(i,j)) > FDT) exit(-1);
			}
		
				
	
		//dnij_dvij
		x = vij;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			

			fd_dnij_dvij_row = (normalise(xp) - normalise(xn))/(2*FD);
			
			for(int i=0; i<3; ++i){
				fd_dnij_dvij(i,j) = fd_dnij_dvij_row(i);
			}
			
		}		
		for(int j=0; j<3; ++j){
			printf("dnij_dvij: %f %f %f \t fd_dnij_dvij: %f %f %f\n",dnij_dvij(j,0),dnij_dvij(j,1),dnij_dvij(j,2),fd_dnij_dvij(j,0),fd_dnij_dvij(j,1),fd_dnij_dvij(j,2));
		}
		
		for(int i=0; i<3; ++i)
			for(int j=0; j<3; ++j){
				if(abs(dnij_dvij(i,j)-fd_dnij_dvij(i,j)) > FDT) exit(-1);
			}
	
			
			
		
		
		
		//domega_dxi
		if(form_i!=SPLITTER){
			x = atoms[index_i];
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + 2*FD;
				xn=x;
				xn(j) = x(j) - 2*FD;
				
				
				vp= xp - torigin;
				lenvp = norm(vp,2);
				np=vp/lenvp;
				if(form_i!=CONVEX)
					np=-np;

				vn= xn - torigin;
				lenvn = norm(vn,2);
				nn=vn/lenvn;
				if(form_i!=CONVEX)
					nn=-nn;
				
				
				omp = (calculateOmega(tessellationOrigin, np, mu_j).rotation);
				omn = (calculateOmega(tessellationOrigin, nn, mu_j).rotation);
				
				printf("OM: %f %f\n",omp,omn);
				

				fd_domega_dxi(j) = sgn(omp-omn)*abs(abs(omp)-abs(omn))/(2*FD);
				
			}		
		}
		else fd_domega_dxi=Vector(3).zeros();
		
		printf("domega_dxi: (%f %f %f) fd_domega_dxi: (%f %f %f)\n",domega_dxi(0),domega_dxi(1),domega_dxi(2),fd_domega_dxi(0),fd_domega_dxi(1),fd_domega_dxi(2));
		for(int j=0; j<3; ++j){
			//if(abs(domega_dxi(j) - fd_domega_dxi(j)) > FDT || isnan(domega_dxi(j)) || isnan(fd_domega_dxi(j))) exit(-1);
		}
		
		
		//domega_dxj
		if(form_j!=SPLITTER){
			x = atoms[index_j];
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				
				vp= xp - torigin;
				lenvp = norm(vp,2);
				np=vp/lenvp;
				if(form_j!=CONVEX)
					np=-np;

				vn= xn - torigin;
				lenvn = norm(vn,2);
				nn=vn/lenvn;
				if(form_j!=CONVEX)
					nn=-nn;
				
				
				omp = (calculateOmega(tessellationOrigin, mu_i, np).rotation);
				omn = (calculateOmega(tessellationOrigin, mu_i, nn).rotation);
				
				
				printf("O: %f %f\n",omp,omn);
				
				fd_domega_dxj(j) = sgn(omp-omn)*abs((abs(omp)-abs(omn)))/(2*FD);
			}		
		}
		else fd_domega_dxj=Vector(3).zeros();
		
		printf("domega_dxj: (%f %f %f) fd_domega_dxj: (%f %f %f)\n",domega_dxj(0),domega_dxj(1),domega_dxj(2),fd_domega_dxj(0),fd_domega_dxj(1),fd_domega_dxj(2));
		for(int j=0; j<3; ++j){
			//if(abs(domega_dxj(j) - fd_domega_dxj(j)) > FDT || isnan(domega_dxj(j)) || isnan(fd_domega_dxj(j))) exit(-1);
		}
		
			
			
			
		
		//domega_dxl
		x = torigin;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
			
			if(form_i!=SPLITTER){
				vp= atoms[index_i] - xp;
				lenvp = norm(vp,2);
				np=vp/lenvp;
				if(form_i!=CONVEX)
					np=-np;
				np3 = np;
			}
			else{
				np = tessellationOrigin;
				np3 = p;
			}

			if(form_i!=SPLITTER){
				vn= atoms[index_i] - xn;
				lenvn = norm(vn,2);
				nn=vn/lenvn;
				if(form_i!=CONVEX)
					nn=-nn;
				nn3 = nn;
			}
			else{
				nn = tessellationOrigin;
				nn3 = p;
			}
			
			
			if(form_j!=SPLITTER){
				vp= atoms[index_j] - xp;
				lenvp = norm(vp,2);
				np2=vp/lenvp;
				if(form_j!=CONVEX)
					np2=-np2;
			}else np2 = tessellationOrigin;

			if(form_j!=SPLITTER){
				vn= atoms[index_j] - xn;
				lenvn = norm(vn,2);
				nn2=vn/lenvn;
				if(form_j!=CONVEX)
					nn2=-nn2;
			}else nn2 = tessellationOrigin;
			
			
			omp = (calculateOmega(tessellationOrigin, np, np2).rotation);
			omn = (calculateOmega(tessellationOrigin, nn, nn2).rotation);
				
			fd_domega_dxl(j) = sgn(omp-omn)*abs((abs(omp)-abs(omn)))/(2*FD);
		}		
		printf("domega_dxl: (%f %f %f) fd_domega_dxl: (%f %f %f)\n",domega_dxl(0),domega_dxl(1),domega_dxl(2),fd_domega_dxl(0),fd_domega_dxl(1),fd_domega_dxl(2));
		for(int j=0; j<3; ++j){
			//if(abs(domega_dxl(j) - fd_domega_dxl(j)) > FDT || isnan(domega_dxl(j)) || isnan(fd_domega_dxl(j))) exit(-1);
		}
	}

	
	return omega;
	
}



Rotation Tessellation::calculateEta(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J){
	printf("CALCULATING ETA FOR: %d %d\n",I.index, J.index);
	return calculateEta(tessellationOrigin, I.normal, J.normal, I.lambda, J.lambda, I.dmu_dx, J.dmu_dx, I.index, J.index, I.form, J.form, true);
}


Rotation Tessellation::calculateEta(Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j){
	Vector tessellationOrigin(3);
	Matrix dmu_dx(3,3);
	
	dmu_dx.zeros();
	
	tessellationOrigin.zeros();
	return calculateEta(tessellationOrigin, mu_i, mu_j, lambda_i, lambda_j, dmu_dx ,dmu_dx, 0, 0, CONVEX, CONVEX, false);
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

Rotation Tessellation::calculateEta(Vector &tessellationOrigin, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, Matrix &dmui_dxi, Matrix &dmuj_dxj, int index_i, int index_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, bool derivatives){
	double rho;
	double dot_IJ;
	Vector drho_dmui(3), drho_dmuj(3);
	Vector drho_dxi(3), drho_dxj(3);
	Vector drho_dxl(3);
	double sig0,sig1,sig2;
	double deta_dlambdai, deta_dlambdaj, deta_drhoij;
	Vector deta_dxi(3), deta_dxj(3), deta_dxl(3);
	Rotation eta;
	Matrix dmui_dxl, dmuj_dxl;
	
	
	
	
	if(derivatives)
		dot_IJ = norm_dot(mu_i,mu_j);
	else
		dot_IJ = dot(mu_i,mu_j);
	rho = acos(dot_IJ);
	
	printf("etav: %f %f %f\n",lambda_i.rotation, lambda_j.rotation, rho);
	
	
	if(isWithinNumericalLimits(rho,0)){
		eta.rotation = 0;
	}
	else{
	
		sig0 = cot(lambda_i.rotation) * cot(rho);
		sig1 = cos(lambda_j.rotation)*csc(lambda_i.rotation)*csc(rho);
		sig2 = sig0-sig1;
		
		printf("SIG: %f %f %f\n",sig0, sig1, sig2);
		
		eta.rotation = acos(sig2);
	}
	
	printf("ETA: %f\n",eta.rotation);
	if(isnan(eta.rotation)) exit(-2);
	
	if(derivatives){
	

		if(isWithinNumericalLimits(rho,0)){
			deta_dxi = Vector(3).zeros();
			deta_dxj = Vector(3).zeros();
			deta_dxl = Vector(3).zeros();
		}
		else{
		
			drho_dmui = -mu_j/(sqrt(1-dot_IJ*dot_IJ));
			drho_dmuj = -mu_i/(sqrt(1-dot_IJ*dot_IJ));
			
			dmui_dxl = -dmui_dxi;
			dmuj_dxl = -dmuj_dxj;
			
			printf("MATRIX FORM %d\n",form_j);
			dmuj_dxj.print("dmuj_dxj");
			
			drho_dxi = (drho_dmui.t() * dmui_dxi).t();
			drho_dxj = (drho_dmuj.t() * dmuj_dxj).t();
			printf("drho_dmui: %f %f %f drho_dmuj: %f %f %f\n",drho_dmui(0),drho_dmui(1),drho_dmui(2),drho_dmuj(0),drho_dmuj(1),drho_dmuj(2));
			for(int i=0; i<3; ++i)
				printf("dmui_dxl: %f %f %f \t\t dmuj_dxl: %f %f %f\n",dmui_dxl(i,0),dmui_dxl(i,1),dmui_dxl(i,2),dmuj_dxl(i,0),dmuj_dxl(i,1),dmuj_dxl(i,2));
			
			drho_dxl = ((drho_dmui.t() * dmui_dxl).t() + (drho_dmuj.t() * dmuj_dxl).t());
			
			deta_dlambdai = csc(lambda_i.rotation) * (sig0/cos(lambda_i.rotation) - sig1*cos(lambda_i.rotation)) / sqrt(1-sig2*sig2);
			deta_dlambdaj =  -sin(lambda_j.rotation)*csc(lambda_i.rotation)*csc(rho) / sqrt(1-sig2*sig2);
			deta_drhoij = csc(rho) * (sig0/cos(rho) - sig1*cos(rho)) / sqrt(1-sig2*sig2);
			
			
			printf("DATA: %f, (%f %f %f), %f, (%f %f %f)\n",deta_dlambdai, lambda_i.drotation_dxi(0), lambda_i.drotation_dxi(1), lambda_i.drotation_dxi(2), deta_drhoij, drho_dxi(0), drho_dxi(1), drho_dxi(2));
			
			deta_dxi = deta_dlambdai * lambda_i.drotation_dxi + deta_drhoij * drho_dxi;
			deta_dxj = deta_dlambdaj * lambda_j.drotation_dxi + deta_drhoij * drho_dxj;
			deta_dxl = deta_dlambdai * lambda_i.drotation_dxl + deta_dlambdaj * lambda_j.drotation_dxl + deta_drhoij * drho_dxl;
			
			eta.drotation_dxi = deta_dxi;
			eta.drotation_dxj = deta_dxj;
			eta.drotation_dxl = deta_dxl;
			
			
			
			
			//finite differences
			Vector fd_drho_dmui(3);
			Vector fd_drho_dmuj(3);
			Vector fd_drho_dxi(3);
			Vector fd_drho_dxj(3);
			Vector fd_drho_dxl(3);
			double fd_deta_dlambdai;
			double fd_deta_dlambdaj;
			double fd_deta_drhoij;
			Vector fd_deta_dxi(3);
			Vector fd_deta_dxj(3);
			Vector fd_deta_dxl(3);
			double sig0p,sig1p,sig2p;
			double sig0n,sig1n,sig2n;
			Vector x(3), xp(3), xn(3);
			Vector np(3), nn(3);
			Vector np2(3), nn2(3);
			Vector vp(3), vn(3);
			double lenvn, lenvp;
			
			
			//drho_dmui
			x = mu_i;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				fd_drho_dmui(j) = (getAngleBetweenNormals(xp, mu_j) -  getAngleBetweenNormals(xn, mu_j))/(2*FD);
			}		
			printf("drho_dmui: (%f %f %f) fd_drho_mui: (%f %f %f)\n",drho_dmui(0),drho_dmui(1),drho_dmui(2),fd_drho_dmui(0),fd_drho_dmui(1),fd_drho_dmui(2));
			for(int j=0; j<3; ++j){
				if(abs(drho_dmui(j) - fd_drho_dmui(j)) > FDT) exit(-1);
			}
			
			
			//drho_dmuj
			x = mu_j;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				fd_drho_dmuj(j) = (getAngleBetweenNormals(mu_i, xp) -  getAngleBetweenNormals(mu_i, xn))/(2*FD);
			}		
			printf("drho_dmuj: (%f %f %f) fd_drho_muj: (%f %f %f)\n",drho_dmuj(0),drho_dmuj(1),drho_dmuj(2),fd_drho_dmuj(0),fd_drho_dmuj(1),fd_drho_dmuj(2));
			for(int j=0; j<3; ++j){
				if(abs(drho_dmuj(j) - fd_drho_dmuj(j)) > FDT) exit(-1);
			}
			

			
			//drho_dxi
			if(form_i!=SPLITTER){
				x = atoms[index_i];
				for(int j=0; j<3; ++j){
					xp=x;
					xp(j) = x(j) + FD;
					xn=x;
					xn(j) = x(j) - FD;
					
					
					vp= xp - torigin;
					lenvp = norm(vp,2);
					np=vp/lenvp;
					if(form_i!=CONVEX)
						np=-np;

					vn= xn - torigin;
					lenvn = norm(vn,2);
					nn=vn/lenvn;
					if(form_i!=CONVEX)
						nn=-nn;
					
					fd_drho_dxi(j) = (getAngleBetweenNormals(np, mu_j) -  getAngleBetweenNormals(nn, mu_j))/(2*FD);
				}		
			}
			else{
				for(int j=0; j<3; ++j)
					fd_drho_dxi(j) = 0;
			}
			printf("drho_dxi: (%f %f %f) fd_drho_dxi: (%f %f %f)\n",drho_dxi(0),drho_dxi(1),drho_dxi(2),fd_drho_dxi(0),fd_drho_dxi(1),fd_drho_dxi(2));
			for(int j=0; j<3; ++j){
				if(abs(drho_dxi(j) - fd_drho_dxi(j)) > FDT) exit(-1);
			}
			
			//drho_dxj
			if(form_j!=SPLITTER){
				x = atoms[index_j];
				for(int j=0; j<3; ++j){
					xp=x;
					xp(j) = x(j) + FD;
					xn=x;
					xn(j) = x(j) - FD;
					
					
					vp= xp - torigin;
					lenvp = norm(vp,2);
					np=vp/lenvp;
					if(form_j!=CONVEX)
						np=-np;

					vn= xn - torigin;
					lenvn = norm(vn,2);
					nn=vn/lenvn;
					if(form_j!=CONVEX)
						nn=-nn;
					
					fd_drho_dxj(j) = (getAngleBetweenNormals(mu_i, np) -  getAngleBetweenNormals(mu_i, nn))/(2*FD);
				}
			}
			else{
				for(int j=0; j<3; ++j)
					fd_drho_dxj(j) = 0;
			}
			
			printf("drho_dxj: (%f %f %f) fd_drho_dxj: (%f %f %f)\n",drho_dxj(0),drho_dxj(1),drho_dxj(2),fd_drho_dxj(0),fd_drho_dxj(1),fd_drho_dxj(2));
			for(int j=0; j<3; ++j){
				if(abs(drho_dxj(j) - fd_drho_dxj(j)) > FDT) exit(-1);
			}

			
			//drho_dxl
			x = torigin;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				
				if(form_i!=SPLITTER){
					vp= atoms[index_i] - xp;
					lenvp = norm(vp,2);
					np=vp/lenvp;
					if(form_i!=CONVEX)
						np=-np;
					
				}
				else np = tessellationOrigin;


				if(form_i!=SPLITTER){
					vn= atoms[index_i] - xn;
					lenvn = norm(vn,2);
					nn=vn/lenvn;
					if(form_i!=CONVEX)
						nn=-nn;
					
				}
				else nn = tessellationOrigin;
				
				
				if(form_j!=SPLITTER){
					vp= atoms[index_j] - xp;
					lenvp = norm(vp,2);
					np2=vp/lenvp;
					if(form_j!=CONVEX)
						np2=-np2;
					
				}
				else np2 = tessellationOrigin;

				if(form_j!=SPLITTER){
					vn= atoms[index_j] - xn;
					lenvn = norm(vn,2);
					nn2=vn/lenvn;
					if(form_j!=CONVEX)
						nn2=-nn2;
				}
				else nn2 = tessellationOrigin;
				
				
				fd_drho_dxl(j) = (getAngleBetweenNormals(np, np2) -  getAngleBetweenNormals(nn, nn2))/(2*FD);
			}
			printf("drho_dxl: (%f %f %f) fd_drho_dxl: (%f %f %f)\n",drho_dxl(0),drho_dxl(1),drho_dxl(2),fd_drho_dxl(0),fd_drho_dxl(1),fd_drho_dxl(2));
			for(int j=0; j<3; ++j){
				if(abs(drho_dxl(j) - fd_drho_dxl(j)) > FDT) exit(-1);
			}
			
			

			//deta_dlambdai
			sig0p = cot(lambda_i.rotation+FD) * cot(rho);
			sig1p = cos(lambda_j.rotation)*csc(lambda_i.rotation+FD)*csc(rho);
			sig2p = sig0p-sig1p;

			sig0n = cot(lambda_i.rotation-FD) * cot(rho);
			sig1n = cos(lambda_j.rotation)*csc(lambda_i.rotation-FD)*csc(rho);
			sig2n = sig0n-sig1n;
			

			fd_deta_dlambdai = (acos(sig2p) - acos(sig2n)) / (2*FD);
			printf("deta_dlambdai: %f FD: %f\n",deta_dlambdai, fd_deta_dlambdai);
			if(abs(deta_dlambdai - fd_deta_dlambdai) > FDT) exit(-1);
			
			//deta_dlambdaj
			sig0p = cot(lambda_i.rotation) * cot(rho);
			sig1p = cos(lambda_j.rotation+FD)*csc(lambda_i.rotation)*csc(rho);
			sig2p = sig0p-sig1p;

			sig0n = cot(lambda_i.rotation) * cot(rho);
			sig1n = cos(lambda_j.rotation-FD)*csc(lambda_i.rotation)*csc(rho);
			sig2n = sig0n-sig1n;
			
			fd_deta_dlambdaj = (acos(sig2p) - acos(sig2n)) / (2*FD);
			printf("deta_dlambdaj: %f FD: %f\n",deta_dlambdaj, fd_deta_dlambdaj);
			if(abs(deta_dlambdaj - fd_deta_dlambdaj) > FDT) exit(-1);
			
			
			//deta_drhoij
			sig0p = cot(lambda_i.rotation) * cot(rho+FD);
			sig1p = cos(lambda_j.rotation)*csc(lambda_i.rotation)*csc(rho+FD);
			sig2p = sig0p-sig1p;

			sig0n = cot(lambda_i.rotation) * cot(rho-FD);
			sig1n = cos(lambda_j.rotation)*csc(lambda_i.rotation)*csc(rho-FD);
			sig2n = sig0n-sig1n;
			
			fd_deta_drhoij = (acos(sig2p) - acos(sig2n)) / (2*FD);
			printf("deta_drhoij: %f FD: %f\n",deta_drhoij, fd_deta_drhoij);
			if(abs(deta_drhoij - fd_deta_drhoij) > FDT) exit(-1);
			
			
			
			

			Vector ai(3), aj(3), al(3);
			Vector nip,nin,nj;
			Vector njp,njn,ni;
			int indexi, indexl, indexj;
			double ri, rl, rj;
			Rotation lambdaip;
			Rotation lambdain;
			Rotation lambdaj;
			double dip,din,dj;
			Vector fd_dPHIij_dxi_in(3);
			Vector fd_dPHIij_dxi_out(3);
			PHIContainer PHIp, PHIn;
			Rotation lambdajp;
			Rotation lambdajn;
			Rotation lambdai;
			double djp,djn,di;
			Vector fd_dPHIij_dxj_in(3);
			Vector fd_dPHIij_dxj_out(3);
			double rhop, rhon;
			


			
			//deta_dxi
			indexi = index_i;
			indexj = index_j;
			if(indexi >= 0 && indexj >= 0){
				ai = atoms[indexi];
				aj = atoms[indexj];
				al = torigin;
				ri = radii[indexi];
				rj = radii[indexj];
				rl = tradius;
				
				printf("ai: %f %f %f aj: %f %f %f\n",ai(0),ai(1),ai(2), aj(0),aj(1),aj(2));
				
				
				x=ai;
				for(int j=0; j<3; ++j){
					xp=x;
					xp(j) = x(j) + FD;
					xn=x;
					xn(j) = x(j) - FD;
					
					nip = calculateInterfaceNormal(al, xp, dip);
					nin = calculateInterfaceNormal(al, xn, din);
					nj = calculateInterfaceNormal(al, aj, dj);
					
					lambdaip = calculateLambda(dip, rl, ri, nip);
					lambdain = calculateLambda(din, rl, ri, nin);
					lambdaj = calculateLambda(dj, rl, rj, nj);
					
					rhop = calculateRho(nip, nj);
					rhon = calculateRho(nin, nj);
					
					sig0p = cot(lambdaip.rotation) * cot(rhop);
					sig1p = cos(lambdaj.rotation)*csc(lambdaip.rotation)*csc(rhop);
					sig2p = sig0p-sig1p;

					sig0n = cot(lambdain.rotation) * cot(rhon);
					sig1n = cos(lambdaj.rotation)*csc(lambdain.rotation)*csc(rhon);
					sig2n = sig0n-sig1n;
				
					fd_deta_dxi(j) = (acos(sig2p) - acos(sig2n)) / (2*FD);
					
				}
				
				printf("deta_dxi: %f %f %f FD: %f %f %f\n",deta_dxi(0),deta_dxi(1),deta_dxi(2), fd_deta_dxi(0), fd_deta_dxi(1), fd_deta_dxi(2));
				for(int i=0; i<3; ++i){
					if(abs(deta_dxi(i) - fd_deta_dxi(i)) > FDT) exit(-1);
				}
				
			}
		
			
			//deta_dxj
			indexi = index_i;
			indexj = index_j;
			if(indexi >= 0 && indexj >= 0){
				ai = atoms[indexi];
				aj = atoms[indexj];
				al = torigin;
				ri = radii[indexi];
				rj = radii[indexj];
				rl = tradius;
				
				printf("ai: %f %f %f aj: %f %f %f\n",ai(0),ai(1),ai(2), aj(0),aj(1),aj(2));
				
				x=aj;
				for(int j=0; j<3; ++j){
					xp=x;
					xp(j) = x(j) + FD;
					xn=x;
					xn(j) = x(j) - FD;
					
					njp = calculateInterfaceNormal(al, xp, djp);
					njn = calculateInterfaceNormal(al, xn, djn);
					ni = calculateInterfaceNormal(al, ai, di);
					
					lambdajp = calculateLambda(djp, rl, rj, njp);
					lambdajn = calculateLambda(djn, rl, rj, njn);
					lambdai = calculateLambda(di, rl, ri, ni);
					
					rhop = calculateRho(ni, njp);
					rhon = calculateRho(ni, njn);
					
					sig0p = cot(lambdai.rotation) * cot(rhop);
					sig1p = cos(lambdajp.rotation)*csc(lambdai.rotation)*csc(rhop);
					sig2p = sig0p-sig1p;

					sig0n = cot(lambdai.rotation) * cot(rhon);
					sig1n = cos(lambdajn.rotation)*csc(lambdai.rotation)*csc(rhon);
					sig2n = sig0n-sig1n;
				
					fd_deta_dxj(j) = (acos(sig2p) - acos(sig2n)) / (2*FD);
					
				}
				
				printf("deta_dxj: %f %f %f FD: %f %f %f\n",deta_dxj(0),deta_dxj(1),deta_dxj(2), fd_deta_dxj(0), fd_deta_dxj(1), fd_deta_dxj(2));
				for(int i=0; i<3; ++i){
					if(abs(deta_dxj(i) - fd_deta_dxj(i)) > FDT) exit(-1);
				}
				
			}
		
			
			//deta_dxl
			indexi = index_i;
			indexj = index_j;
			if(indexi >= 0 && indexj >= 0){
				ai = atoms[indexi];
				aj = atoms[indexj];
				al = torigin;
				ri = radii[indexi];
				rj = radii[indexj];
				rl = tradius;
				
				x=al;
				for(int j=0; j<3; ++j){
					xp=x;
					xp(j) = x(j) + FD;
					xn=x;
					xn(j) = x(j) - FD;
					
					nip = calculateInterfaceNormal(xp, ai, dip);
					nin = calculateInterfaceNormal(xn, ai, din);
					njp = calculateInterfaceNormal(xp, aj, djp);
					njn = calculateInterfaceNormal(xn, aj, djn);
					
					lambdaip = calculateLambda(dip, rl, ri, nip);
					lambdain = calculateLambda(din, rl, ri, nin);
					lambdajp = calculateLambda(djp, rl, rj, njp);
					lambdajn = calculateLambda(djn, rl, rj, njn);
					
					rhop = calculateRho(nip, njp);
					rhon = calculateRho(nin, njn);
					
					sig0p = cot(lambdaip.rotation) * cot(rhop);
					sig1p = cos(lambdajp.rotation)*csc(lambdaip.rotation)*csc(rhop);
					sig2p = sig0p-sig1p;

					sig0n = cot(lambdain.rotation) * cot(rhon);
					sig1n = cos(lambdajn.rotation)*csc(lambdain.rotation)*csc(rhon);
					sig2n = sig0n-sig1n;
				
					fd_deta_dxl(j) = (acos(sig2p) - acos(sig2n)) / (2*FD);
					
				}
				
				printf("deta_dxl: %f %f %f FD: %f %f %f\n",deta_dxl(0),deta_dxl(1),deta_dxl(2), fd_deta_dxl(0), fd_deta_dxl(1), fd_deta_dxl(2));
				for(int i=0; i<3; ++i){
					if(abs(deta_dxl(i) - fd_deta_dxl(i)) > FDT) exit(-1);
				}
				
			}
		
		}
	}
	
	
	return eta;
}

PHIContainer Tessellation::calculatePHI(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J, double dij, double radius){
	PHIContainer p,p2;	
	p =  calculatePHI(tessellationOrigin, I.normal, J.normal, I.lambda, J.lambda, I.dmu_dx, J.dmu_dx, I.index, J.index, I.form, J.form, true);
	p2 = retrieveInterfaces(tessellationOrigin, I, J, dij, radius);
	
	p.in.vector = p2.in.vector;
	p.out.vector = p2.out.vector;
	
	printf("ROTATIONAL INFO: in p:%f, p2:%f  out: p:%f p2:%f\n",p.in.rotation, p2.in.rotation, p.out.rotation, p2.out.rotation);
	
	if(!isWithinNumericalLimits(p.in.rotation,p2.in.rotation) || !isWithinNumericalLimits(p.out.rotation,p2.out.rotation)){
	printf("I and J: c(%f,%f,%f),c(%f,%f,%f)\n",I.normal(0),I.normal(1),I.normal(2),J.normal(0),J.normal(1),J.normal(2));
		
		printf("ROTATIONAL ERROR: in p:%f, p2:%f  out: p:%f p2:%f\n",p.in.rotation, p2.in.rotation, p.out.rotation, p2.out.rotation);
		//exit(-1);
	}
	
	return p;
	
}


PHIContainer Tessellation::calculatePHI(Vector &tessellationOrigin, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j){
	Matrix dmu_dx(3,3);
	dmu_dx.zeros();
	return calculatePHI(tessellationOrigin, mu_i, mu_j, lambda_i, lambda_j, dmu_dx, dmu_dx, 0, 0, CONVEX, CONVEX, false);
}


PHIContainer Tessellation::calculatePHI(Vector &tessellationOrigin, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, Matrix &dmui_dxi, Matrix &dmuj_dxj, int index_i, int index_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, bool derivatives){
	PHIContainer p,p2,p3;
	
	Rotation eta,omega;
	int q;
	double S1;
	
	
	//eta = M_PI/2 - acos(-(1/tan(lambda_j))*(1/tan(rho))+cos(lambda_k)*(1/sin(lambda_j))*(1/sin(rho)));
	
	eta = calculateEta(tessellationOrigin, mu_i, mu_j, lambda_i, lambda_j, dmui_dxi, dmuj_dxj, index_i, index_j, form_i, form_j, derivatives);
	omega = calculateOmega(tessellationOrigin, mu_i, mu_j, dmui_dxi, dmuj_dxj, index_i, index_j, form_i, form_j, derivatives);
	
	//printf("OMEGA: %f ETA: %f\n",omega.rotation, eta.rotation);

	S1 = sgn(M_PI-p.out.rotation);

	
	p.out.rotation = omega.rotation + eta.rotation;
	if(p.out.rotation > M_PI)
		p.out.rotation = -2*M_PI + p.out.rotation;
	
	p.in.rotation = omega.rotation - eta.rotation;
	if(p.in.rotation < -M_PI)
		p.in.rotation = 2*M_PI + p.in.rotation;

	
	
	p.out.rotation2 = omega.rotation2 + eta.rotation;
	if(p.out.rotation2 > M_PI)
		p.out.rotation2 = -2*M_PI + p.out.rotation2;
	
	p.in.rotation2 = omega.rotation2 - eta.rotation;
	if(p.in.rotation2 < -M_PI)
		p.in.rotation2 = 2*M_PI + p.in.rotation2;

	
	
	if(derivatives){
		
		//if(p.out.rotation>=M_PI) p.out.rotation = (eta.rotation + omega.rotation) - 2*M_PI;
		//else if(p.out.rotation<=-M_PI) p.out.rotation = (eta.rotation + omega.rotation) + 2*M_PI;
		p.out.drotation_dxi = eta.drotation_dxi + omega.drotation_dxi;
		p.out.drotation_dxj = eta.drotation_dxj + omega.drotation_dxj;
		p.out.drotation_dxl = eta.drotation_dxl + omega.drotation_dxl;
		
		
		p.in.drotation_dxi = -eta.drotation_dxi + omega.drotation_dxi;
		p.in.drotation_dxj = -eta.drotation_dxj + omega.drotation_dxj;
		p.in.drotation_dxl = -eta.drotation_dxl + omega.drotation_dxl;
		
		if(form_i == SPLITTER){
			if(!isWithinNumericalLimits(dot(p.in.drotation_dxi, p.out.drotation_dxj),0.0)){
				printf("SPLITTER CHECK: I in (%f %f %f) out (%f %f %f)\n", p.in.drotation_dxi(0), p.in.drotation_dxi(1), p.in.drotation_dxi(2), p.out.drotation_dxi(0), p.out.drotation_dxi(1), p.out.drotation_dxi(2));
				exit(-4);
			}
		}
		
		if(form_j == SPLITTER){
			if(!isWithinNumericalLimits(dot(p.in.drotation_dxj, p.out.drotation_dxj),0.0)){
				printf("SPLITTER CHECK: J in (%f %f %f) out (%f %f %f)\n", p.in.drotation_dxj(0), p.in.drotation_dxj(1), p.in.drotation_dxj(2), p.out.drotation_dxj(0), p.out.drotation_dxj(1), p.out.drotation_dxj(2));
				exit(-4);
			}
		}
		

		
		

		Vector ai(3), aj(3), al(3), x(3), xp(3), xn(3);
		Vector nip,nin,nj;
		Vector njp,njn,ni;
		int indexi, indexl, indexj;
		double ri, rl, rj;
		Rotation lambdaip;
		Rotation lambdain;
		Rotation lambdaj;
		double dip,din,dj;
		Vector fd_dPHIij_dxi_in(3);
		Vector fd_dPHIij_dxi_out(3);
		PHIContainer PHIp, PHIn;
		Rotation lambdajp;
		Rotation lambdajn;
		Rotation lambdai;
		double djp,djn,di;
		Vector fd_dPHIij_dxj_in(3);
		Vector fd_dPHIij_dxj_out(3);
		


		
		//dphi_dxi
		indexi = index_i;
		indexj = index_j;
		if(indexi >= 0 && indexj >= 0){
			ai = atoms[indexi];
			aj = atoms[indexj];
			al = torigin;
			ri = radii[indexi];
			rj = radii[indexj];
			rl = tradius;
			
			
			nip = calculateInterfaceNormal(al, ai, dip);
			nj = calculateInterfaceNormal(al, aj, dj);
			
			lambdaip = calculateLambda(dip, rl, ri, nip);
			lambdaj = calculateLambda(dj, rl, rj, nj);
			
			PHIp = calculatePHI(tessellationOrigin, nip, nj, lambdaip, lambdaj);
			
			printf("PRECOMPARE: %f %f / %f %f\n",p.in.rotation, PHIp.in.rotation, p.out.rotation, PHIp.out.rotation);
			
			
			
			x=ai;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				nip = calculateInterfaceNormal(al, xp, dip);
				nin = calculateInterfaceNormal(al, xn, din);
				nj = calculateInterfaceNormal(al, aj, dj);
				
				lambdaip = calculateLambda(dip, rl, ri, nip);
				lambdain = calculateLambda(din, rl, ri, nin);
				lambdaj = calculateLambda(dj, rl, rj, nj);
				
				PHIp = calculatePHI(tessellationOrigin, nip, nj, lambdaip, lambdaj);
				PHIn = calculatePHI(tessellationOrigin, nin, nj, lambdain, lambdaj);
				
				fd_dPHIij_dxi_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
				fd_dPHIij_dxi_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
			}
			printf("dPHIij_dxi_in: (%f %f %f) FD: (%f, %f, %f) dPHIij_dxi_out: (%f %f %f) FD: (%f %f %f)\n",p.in.drotation_dxi(0),p.in.drotation_dxi(1),p.in.drotation_dxi(2), fd_dPHIij_dxi_in(0), fd_dPHIij_dxi_in(1), fd_dPHIij_dxi_in(2), p.out.drotation_dxi(0),p.out.drotation_dxi(1),p.out.drotation_dxi(2), fd_dPHIij_dxi_out(0), fd_dPHIij_dxi_out(1), fd_dPHIij_dxi_out(2));
			for(int i=0; i<3; ++i){
				//if(abs(p.in.drotation_dxi(i) - fd_dPHIij_dxi_in(i)) > FDT) exit(-1);
				//if(abs(p.out.drotation_dxi(i) - fd_dPHIij_dxi_out(i)) > FDT) exit(-1);
			}
			
		}
	
		
		
		
		
		//dphi_dxj
		indexi = index_i;
		indexj = index_j;
		if(indexi >= 0 && indexj >= 0){
			ai = atoms[indexi];
			aj = atoms[indexj];
			al = torigin;
			ri = radii[indexi];
			rj = radii[indexj];
			rl = tradius;
			
			
			x=aj;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				njp = calculateInterfaceNormal(al, xp, djp);
				njn = calculateInterfaceNormal(al, xn, djn);
				ni = calculateInterfaceNormal(al, ai, di);
				
				lambdajp = calculateLambda(djp, rl, rj, njp);
				lambdajn = calculateLambda(djn, rl, rj, njn);
				lambdai = calculateLambda(di, rl, ri, ni);
				
				PHIp = calculatePHI(tessellationOrigin, ni, njp, lambdai, lambdajp);
				PHIn = calculatePHI(tessellationOrigin, ni, njn, lambdai, lambdajn);
				
				fd_dPHIij_dxj_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
				fd_dPHIij_dxj_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
			}
			printf("dPHIij_dxj_in: (%f %f %f) FD: (%f, %f, %f) dPHIij_dxj_out: (%f %f %f) FD: (%f %f %f)\n",p.in.drotation_dxj(0),p.in.drotation_dxj(1),p.in.drotation_dxj(2), fd_dPHIij_dxj_in(0), fd_dPHIij_dxj_in(1), fd_dPHIij_dxj_in(2), p.out.drotation_dxj(0),p.out.drotation_dxj(1),p.out.drotation_dxj(2), fd_dPHIij_dxj_out(0), fd_dPHIij_dxj_out(1), fd_dPHIij_dxj_out(2));
			for(int i=0; i<3; ++i){
				//if(abs(p.in.drotation_dxj(i) - fd_dPHIij_dxj_in(i)) > FDT) exit(-1);
				//if(abs(p.out.drotation_dxj(i) - fd_dPHIij_dxj_out(i)) > FDT) exit(-1);
			}
		}
			
		
		
		
		
	}
	
	
	q = 1;
	if(form_i != CONVEX) q*=-1;
	if(form_j != CONVEX) q*=-1;
	
	
	
	
	if(q<0){
		printf("EXCHANGE\n");
		p3=p;
		p.in=p3.out;
		p.out=p3.in;
	}



	
	return p;
	
	
	 
}













Rotation Tessellation::calculatePsi(Vector &tessellationOrigin, CircularInterface &circle){
	return calculatePsi(tessellationOrigin, circle.normal, circle.dmu_dx, circle.form, circle.index, true);
}

Rotation Tessellation::calculatePsi(Vector &tessellationOrigin, Vector &mu_i){
	Matrix dmu_dx(3,3);
	CircularInterfaceForm form;
	int index;
	
	dmu_dx.zeros();
	form = CONVEX;
	index=0;
	
	return calculatePsi(tessellationOrigin, mu_i, dmu_dx, form, index, false);
}

Rotation Tessellation::calculatePsi(Vector &tessellationOrigin, Vector &mu_i, Matrix &dmu_dx, CircularInterfaceForm form, int index, bool derivatives){
	Vector dpsi_dmui;
	Rotation r;
	Vector dpsi_dxi, dpsi_dxl;
	
	Vector n_i, nn_i;
	double psi2;
	Vector ey(3);
	Matrix I(3,3);
	I.eye();
	
	double dotmui_nni; 
	Matrix dni_dmui;
	Matrix dnni_dni;
	Vector dpsi_dnni;

	Matrix dni_dxi;
	Matrix dnni_dxi;
	double acos_dotmui_nni;
	double q0,q1,q2,q3;
	
	
	r.drotation_dxi = Vector(3).zeros();
	r.drotation_dxj = Vector(3).zeros();
	r.drotation_dxl = Vector(3).zeros();
	
	r.rotation = getAngle(tessellationOrigin,mu_i);
	r.rotation2 = r.rotation;
	
	/*
	if((r.rotation < M_PI/4 || r.rotation > 3*M_PI/4)){
		if(isWithinNumericalLimits(r.rotation,0)){
			ey.zeros();
			ey(1) = 1;
			n_i = cross(tessellationOrigin, ey);
			nn_i = cross(n_i, tessellationOrigin);
		}
		else{
			n_i = cross(tessellationOrigin, mu_i);
			nn_i = cross(n_i, tessellationOrigin);
		
		}
		
		//dotmui_nni = norm_dot(nn_i , mu_i);
		dotmui_nni = dot(nn_i , mu_i);
		
		q0=1;
		q1=1;
		q2=0;
		q3=-1;
		
		
		acos_dotmui_nni = acos(dotmui_nni);
		
		if(acos_dotmui_nni <= M_PI/2) q0=-1;
		if(r.rotation >= M_PI/2) q1=-1;
		if(r.rotation >= M_PI/4 && acos_dotmui_nni >= M_PI/2){
			q2 = 1;
			q3 = 1;
			
		}
		
		
		
		psi2 = q3*(q2*M_PI - (q0*(q1*acos_dotmui_nni - M_PI/2)));
		
		r.rotation2 = psi2;	
		
	}
*/
	
	
	
	
	
	if(derivatives){
	
		if(form != SPLITTER){
			
			
			if(false && (r.rotation < M_PI/4 || r.rotation > 3*M_PI/4)){
				if(isWithinNumericalLimits(r.rotation,0)){
					ey.zeros();
					ey(1) = 1;
					n_i = cross(tessellationOrigin, ey);
					nn_i = cross(n_i, tessellationOrigin);
					printf("ZERO\n");
				}
				else{
					n_i = cross(tessellationOrigin, mu_i);
					nn_i = cross(n_i, tessellationOrigin);
					printf("NORMAL\n");
				
				}
				
				//dotmui_nni = norm_dot(nn_i , mu_i);
				dotmui_nni = dot(nn_i , mu_i);
				
				q0=1;
				q1=1;
				q2=0;
				q3=-1;
				
				
				acos_dotmui_nni = acos(dotmui_nni);
				
				if(acos_dotmui_nni <= M_PI/2) q0=-1;
				if(r.rotation >= M_PI/2) q1=-1;
				if(r.rotation >= M_PI/4 && acos_dotmui_nni >= M_PI/2){
					q2 = 1;
					q3 = 1;
					
				}
				
				printf("acos %f\n",acos_dotmui_nni);
				printf("q: %f %f %f %f\n",q0,q1,q2,q3);
				
				
				
				psi2 = q3*(q2*M_PI - (q0*(q1*acos_dotmui_nni - M_PI/2)));
				
				r.rotation2 = psi2;
				
				printf("psi2: %f %f\n",r.rotation, psi2);
				printf("to: %f %f %f\n",tessellationOrigin(0),tessellationOrigin(1),tessellationOrigin(2));
				//if(!isWithinNumericalLimits(r.rotation,psi2)) exit(-1);
				
				
				
				dni_dmui = -matrixCross(I,tessellationOrigin);
				dnni_dni = matrixCross(I, tessellationOrigin);
				dpsi_dnni = -2*q0*q1*q3*mu_i / sqrt(1-dotmui_nni*dotmui_nni);
				dpsi_dmui = -2*q0*q1*q3*nn_i / sqrt(1-dotmui_nni*dotmui_nni);

				dni_dxi = dni_dmui * dmu_dx;
				dnni_dxi = dnni_dni * dni_dxi;
				dpsi_dxi = 0.5*((dpsi_dnni.t() * dnni_dxi).t() + (dpsi_dmui.t() * dmu_dx).t());
				
				dpsi_dxl = -dpsi_dxi;
				
				
				
			}
			else{

				/*
				if(isWithinNumericalLimits(mu_i(0),1.0)){
					printf("ZERO LEVELII\n");
					mu_i(0) = 1.0-MINISCULE;
				}
				if(isWithinNumericalLimits(mu_i(0),-1.0)){
					printf("ZERO LEVELIII\n");
					mu_i(0) = -1.0+MINISCULE;
				}
				*/
				
				
				if(isWithinNumericalLimits(mu_i(0),1.0) || isWithinNumericalLimits(mu_i(0),-1.0)){
					printf("ZERO LEVELIII\n");
					dpsi_dmui = Vector(3).zeros();
					dpsi_dxi = Vector(3).zeros();
					dpsi_dxl = Vector(3).zeros();
				}
				else{
					
					dpsi_dmui = -tessellationOrigin/(sqrt(1-mu_i(0)*mu_i(0)));
					dpsi_dxi = (dpsi_dmui.t() * dmu_dx).t();
					dpsi_dxl = -dpsi_dxi;
				}
			}
			
			
			
			
			
			
			r.drotation_dxi = dpsi_dxi;
			r.drotation_dxj = Vector(3).zeros();
			r.drotation_dxl = dpsi_dxl;
		}
		else{
			//the splitter is immovable
			r.drotation_dxi = Vector(3).zeros();
			r.drotation_dxj = Vector(3).zeros();
			r.drotation_dxl = Vector(3).zeros();
		}
		
	

		//finite differences
		Vector fd_dpsi_dmui(3);
		Vector fd_dpsi_dxi(3);
		Vector fd_dpsi_dxl(3);
		Vector x(3), xp(3), xn(3);
		Vector np(3), nn(3);
		Vector vp(3), vn(3);
		double lenvn, lenvp;
		Vector nip, nnip, nin, nnin;
		double a0,a1;
		
		
		if(form != SPLITTER){
			//dpsi_dmui
			

				x = mu_i;
				for(int j=0; j<3; ++j){
					xp=x;
					xp(j) = x(j) + FD;
					xn=x;
					xn(j) = x(j) - FD;
					if(abs(dot(tessellationOrigin, xp)) > 1.0 || abs(dot(tessellationOrigin, xn)) > 1.0){
						xp=x;
						xp(j) = x(j) + 2 * FD;
						xn=x;
					}
					if(abs(dot(tessellationOrigin, xp)) > 1.0 || abs(dot(tessellationOrigin, xn)) > 1.0){
						xp=x;
						xn=x;
						xn(j) = x(j) - 2 * FD;
					}
					
					fd_dpsi_dmui(j) = (getAngleBetweenNormals(tessellationOrigin, xp) -  getAngleBetweenNormals(tessellationOrigin, xn))/(2*FD);
				}		
				printf("dpsi_dmui: (%f %f %f) fd_dpsi_mui: (%f %f %f)\n",dpsi_dmui(0),dpsi_dmui(1),dpsi_dmui(2),fd_dpsi_dmui(0),fd_dpsi_dmui(1),fd_dpsi_dmui(2));
				for(int j=0; j<3; ++j){
					//if(abs(dpsi_dmui(j) - fd_dpsi_dmui(j)) > FDT || isnan(dpsi_dmui(j)) || isnan(fd_dpsi_dmui(j)) ) exit(-1);
				}
		
		
		
			//dpsi_dxi
				x = atoms[index];
				for(int j=0; j<3; ++j){
					xp=x;
					xp(j) = x(j) + FD;
					xn=x;
					xn(j) = x(j) - FD;
					
					
					vp= xp - torigin;
					lenvp = norm(vp,2);
					np=vp/lenvp;
					if(form!=CONVEX)
						np=-np;

					vn= xn - torigin;
					lenvn = norm(vn,2);
					nn=vn/lenvn;
					if(form!=CONVEX)
						nn=-nn;
					
					a0 = calculatePsi(tessellationOrigin, np).rotation;
					a1 = calculatePsi(tessellationOrigin, nn).rotation;
					fd_dpsi_dxi(j) =  (a0-a1)/(2*FD);
				}		
				
				printf("dpsi_dxi: (%f %f %f) fd_dpsi_dxi: (%f %f %f)\n",dpsi_dxi(0),dpsi_dxi(1),dpsi_dxi(2),fd_dpsi_dxi(0),fd_dpsi_dxi(1),fd_dpsi_dxi(2));
				for(int j=0; j<3; ++j){
				//	if(abs(dpsi_dxi(j) - fd_dpsi_dxi(j)) > FDT) exit(-1);
				}
		
			//dpsi_dxl
				x = torigin;
				for(int j=0; j<3; ++j){
					xp=x;
					xp(j) = x(j) + FD;
					xn=x;
					xn(j) = x(j) - FD;
					
					
					vp= atoms[index] - xp;
					lenvp = norm(vp,2);
					np=vp/lenvp;
					if(form!=CONVEX)
						np=-np;
					

					vn= atoms[index] - xn;
					lenvn = norm(vn,2);
					nn=vn/lenvn;
					if(form!=CONVEX)
						nn=-nn;
					
					fd_dpsi_dxl(j) = (getAngleBetweenNormals(tessellationOrigin, np) -  getAngleBetweenNormals(tessellationOrigin, nn))/(2*FD);
				}		
				printf("dpsi_dxl: (%f %f %f) fd_dpsi_dxl: (%f %f %f)\n",dpsi_dxl(0),dpsi_dxl(1),dpsi_dxl(2),fd_dpsi_dxl(0),fd_dpsi_dxl(1),fd_dpsi_dxl(2));
				for(int j=0; j<3; ++j){
				//	if(abs(dpsi_dxl(j) - fd_dpsi_dxl(j)) > FDT) exit(-1);
				}
		}
		
	
		
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
				printf("ACTIVATING INDEX: %d AND %d (%d,%d)\n",sasaNode.index0, sasaNode.index1, x.id0, x.id1);
				sasaNode.rotation0 = intersectionGraph[x].rotation0;
				sasaNode.rotation1 = intersectionGraph[x].rotation1;
				sasaNode.lambda = circles[cid0].lambda;
				sasaNode.psi = circles[cid0].psi;
				sasaNode.vector = intersectionGraph[x].rotation1.vector;
				sasaNode.normalForCircularInterface = circles[cid0].normal;
				sasaNode.form = circles[cid0].form;
				sasaNode.form0 = circles[cid0].form;
				sasaNode.form1 = circles[cid1].form;
				
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
			//fprintf(file, "intersector %d radius %f vector %f %f %f\n", i, circle.sphereRadius, circle.vector(0), circle.vector(1), circle.vector(2));
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
