#include "tessellation.h"


Tessellation::Tessellation(Molecule &m){
	molecule = m;
}


vector<vector<double*> > &Tessellation::fetchForcePointers(){
	return forces;
}

vector<double*> &Tessellation::fetchAreaPointers(){
	return areas;
}



void Tessellation::build(bool split){
	CircularInterfacesPerAtom circlesPerAtom;
	sasasForMolecule.clear();
	
	//molecule.update();
	atoms = molecule.fetchCoordinates();
	radii = molecule.fetchRadii();
	forces = molecule.fetchForcePointers();
	areas = molecule.fetchAreaPointers();
	
	
	
	
	//atoms.size()
	//iterate over all atoms and build the tessellation for each of them
	for(int i=0; i<atoms.size(); ++i){
		buildGaussBonnetPath(i, atoms, radii, sasasForMolecule, split);
		
		
		
		//for(int i=0; i<sasasForMolecule[0].sasas.size(); ++i){
		//	printf("SASA %d\n",i);
		//	outputGaussBonnetPath(sasasForMolecule[0].sasas[i]);
		//}
		//outputCircularInterfaces(*circles[circles.size()-1]);

	}
	
	
	
}


CircularInterfacesPerAtom Tessellation::coverHemisphere(Vector tessellationOrigin, double radius, CircularInterfacesPerAtom circles, CircularInterfaceForm form){
	CircularInterface C;
	Vector v(3);
	
	v = -tessellationOrigin;
	
	C.lambda = M_PI * 0.5;
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
	CircularInterfacesPerAtom circles;
	CircularInterfacesPerAtom circlesFrontHemisphere;
	CircularInterfacesPerAtom circlesBackHemisphere;
	Vector frontTessellationOrigin(3);
	Vector backTessellationOrigin(3);
	Vector origin;
	double radius;
	SASAsForAtom sasasForAtom;
	SASAs *newSasas;
	
	
	origin = atoms[i];
	radius = radii[i];
	
	
	frontTessellationOrigin(0) = 1;
	frontTessellationOrigin(1) = 0;
	frontTessellationOrigin(2) = 0;
	backTessellationOrigin = -frontTessellationOrigin;

	//set the random number generator to some starting value
	srand(2);

	
	sasasForAtom.radius = radius;
	sasasForAtom.vector = origin;
	
	newSasas = &(sasas.insert(sasas.end(),sasasForAtom)->sasas);
	
	//in this block, basic qantities of intersecting spheres are calculated (i.e. size of the interface)
	makeCircularInterfaces(i,origin, radius, atoms, radii, circles);
	for(j=0;j<circles.size();j++){
		determineProjection(origin, radius, circles[j]);
	}
	
	//there is room for optimisation here...
	if(split){
	
		//we add a concave splitter interface, thereby cutting the sphere in two
		circlesFrontHemisphere = coverHemisphere(frontTessellationOrigin, radius, circles, SPLITTER);
		circlesBackHemisphere = coverHemisphere(backTessellationOrigin, radius, circles, SPLITTER);
		
		//filter out completely occluded interfaces
		filterCircularInterfaces(frontTessellationOrigin,radius, circlesFrontHemisphere);
		filterCircularInterfaces(backTessellationOrigin,radius, circlesBackHemisphere);
		
		//we deleted some interfaces, so we should readjust indexes
		reindexCircularInterfaces(circlesFrontHemisphere);
		reindexCircularInterfaces(circlesBackHemisphere);
		
		//now we find intersections between colliding interfaces
		determineCircularIntersections(circlesFrontHemisphere);
		determineCircularIntersections(circlesBackHemisphere);
		
		//we have all information to actually build boundary of the sasa
		buildIntersectionGraph(radius, frontTessellationOrigin, circlesFrontHemisphere, *newSasas, FRONTHEMISPHERE, string("gbonnet0.csv"));
		buildIntersectionGraph(radius, backTessellationOrigin, circlesBackHemisphere, *newSasas, BACKHEMISPHERE, string("gbonnet1.csv"));
	}
	else{
		printf("non-splitting is actually not supported at the moment...");
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




void Tessellation::determineProjection(Vector &tessellationOrigin, Vector &origin, double radius, CircularInterface &circle){
	double d_i;
	double r_l, r_i;
	double g;
	Vector mu;
	double a;
	Matrix Identity(3,3);
	Identity.eye();
	
	
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
	circle.g_normalised = g/radius;
	circle.lambda.rotation = acos(g_normalised);
	circle.d=d_i;
	
	//derivatives
	dlambda_dg = -1/sqrt(1-g_normalised*g_normalised);
	dg_dd = -g_normalised/d_i + 1/r_l;
	
	
	circle.lambda.drotation_dxi = dlambda_dg * dg_dd * circle.normal;
	circle.lambda.drotation_dxj = 0;
	circle.lambda.drotation_dxdelta = -circle.lambda.drotation_dxi;
	
	
	
	dmui_dxi = 1.0/I.d*(Identity - kron(I.normal,I.normal.t());
	
	circle.dmui_dxi = dmui_dxi;
	
	dpsi_dmui = -1.0/sqrt(1-circle.normal(0)*circle.normal(0));
	dpsi_dxi = dpsi_dmui * dmui_dxi;
	dpsi_dxdelta = -dpsi_dxi;
	
	circle.psi.rotation = getAngleBetweenNormals(tessellationOrigin,circle.normal);
	circle.psi.drotation_dxi = dpsi_dxi;
	circle.psi.drotation_dxj = 0;
	circle.psi.drotation_dxdelta = dpsi_dxdelta;
	
	
	

	//printf("DET: [%f %f %f] %f %f %f %f\n",circle.normal(0),circle.normal(1),circle.normal(2), r_i, r_k, g, circle.lambda);
	//printf("%f %f %f %f %f\n",d_k, r_i, r_k, g, circle.lambda);
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
				circle.index = j;
				circle.vector = v;
				circle.normal = normal;
				circle.sphereRadius = r_k;
				circle.intersect = false;
				circles.push_back(circle);
				
				//printf("CIRCLE[%d]: (%f, %f, %f) %f\n",circle.id, circle.vector(0), circle.vector(1), circle.vector(2), circle.sphereRadius);
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
	
		if(it->form==CONVEX) n0 = it->normal;
		else n0 = -it->normal;
		
		for(i=0;i<circles.size();i++){
			if(it->id != circles[i].id){

				//in case the interface is not convex, we have to reverse the normal (since we already reversed it one time when we calculated the interface)
				if(circles[i].form == CONVEX) n1 = circles[i].normal;
				else n1 = -circles[i].normal;
				
				
				
				angle = getAngleBetweenNormals(n0,n1);
				
				if(it->form==CONVEX){
					//convex circle IT is inside of convex circle i
					if(circles[i].form == CONVEX){
						
						if(it->lambda + angle < circles[i].lambda){
							it = circles.erase(it);
							erased=true;
							break;
						}
					}
					//convex circle is outside of concave circle i
					else{
						if(angle-it->lambda >  circles[i].lambda){
							it = circles.erase(it);
							erased=true;
							break;
						}
					}
					
				}
				else{
					//concave circle IT has a free area. This area is covered by convex circle i
					if(circles[i].form == CONVEX){
						if(it->lambda + angle < circles[i].lambda){
							depleteCircularInterfaces(tessellationOrigin, radius, circles);
							return -1;
						}
						
					}
					else{
						//concave circle IT is completely inside the occlusion of concave circle i
						if(it->lambda > circles[i].lambda + angle){
							it = circles.erase(it);
							erased=true;
							break;
						}
						//concave circle IT has a free area. This area is covered by concave circle i
						else if(angle > it->lambda + circles[i].lambda){
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
	Interfaces f;
	
	
	
	
	//this will create two points on the border of the circular region on opposite sides.
	f.in = M_PI/2;
	f.out = -M_PI/2;

	
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
				if(angle < circles[k].lambda + circles[j].lambda)
					if(angle + circles[k].lambda > circles[j].lambda && angle + circles[j].lambda > circles[k].lambda){
						c.d=angle;
						c.visited=false;
						circles[j].circularIntersections[k]=c;
						circles[k].circularIntersections[j]=c;
						
						//printf("circle intersection %d-%d\n",k,j);
					}
			}
		}
	
}




Rotation Tessellation::calculateOmega(Vector &tessellationOrigin, CircularInterface &I, CircularInterface &J){
	Vector ex(3);
	Vector ez(3);
	Vector ni(3);
	Vector nij(3);
	Vector nii(3);
	double sn;
	double omega, varpi;
	double rho;
	Matrix Identity(3,3);
	Identity.eye();
	Rotation omega;
	
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
		printf("POSITIVE RANGE -1 %f\n",d);
	}
	else ni = cross(tessellationOrigin,I.normal);
	
	nii = cross(I.normal,ni);
	
	nij = cross(I.normal, J.normal);
	dot_ni_nij = norm_dot(nij,nii)
	
	s0 = sgn(dot_ni_nij);
	varpi = asin(norm_dot(ni,nij));
	s1 = sgn(a);
	
	omega.rotation = -s0*(s1*(1-s0)*M_PI/2 - varpi);
	
	if(I.form != CONVEX){
		if(omega>=0)
			omega = -M_PI + omega;
	}
	
	
	dmui_dxi = I.dmu_dx;
	dmui_dxdelta = -I.dmu_dx;
	dmuj_dxj = J.dmu_dx;
	dmuj_dxdelta = -J.dmu_dx;
	
	domega_dvarpi = S0;
	dvarpi_dni = nij / sqrt(1-dot_ni_nij*dot_ni_nij);
	dvarpi_dnij = ni / sqrt(1-dot_ni_nij*dot_ni_nij);
	
	dni_dmui=Matrix(3,3);
	dni_dmui.zeros();
	dni_dmui(2,3)=1;
	dni_dmui(3,2)=-1;
	
	dni_muj = -dni_mui;
	
	dnij_dmui = matrixCross(Identity,J.normal);
	dnij_dmuj = matrixCross(-1*Identity,I.normal);
	
	domega_dxi = domega_dvarpi * (dvarpi_dni * dni_dmui * dmui_dxi + dvarpi_dnij * dnij_dmui * dmui_dxi);
	domega_dxj = domega_dvarpi * dvarpi_dnij * dnij_dmuj * dmui_dxj;

	//domega_dxdelta = -domega_dxi-domega_dxj;

	domega_dxdelta = domega_dvarpi * (dvarpi_dni * dni_dmui * dmui_dxi + dvarpi_dnij * (dnij_dmui * dmui_dxdelta + dnij_dmuj * dmuj_dxdelta);
	
	omega.drotation_dxi = domega_dxi;
	omega.drotation_dxj = domega_dxj;
	omega.drotation_dxdelta = domega_dxdelta;
	
	
	return omega;
	
}


Rotation Tessellation::calculateEta(CircularInterface &I, CircularInterface &J){
	double lambda_j, lambda_k, rho;
	double dot_IJ;
	double drho_mui, drho_muj;
	double drho_dxi, drho_dxj;
	double drho_dxdelta;
	double sig0,sig1,sig2;
	double deta_dlambdai, deta_dlambdaj;
	double deta_drhoij;
	double deta_dxi, deta_dxj, deta_dxdelta;
	Rotation eta;
	
	
	lambda_j = I.lambda.rotation;
	lambda_k = J.lambda.rotation;
	dot_IJ = norm_dot(I.normal,J.normal);
	rho = acos(dot_IJ);
	
	drho_mui = -J.normal/sqrt(1-dotIJ*dotIJ);
	drho_muj = -I.normal/sqrt(1-dotIJ*dotIJ);
	
	
	dmui_dxi = I.dmu_dx;
	dmui_dxj = J.dmu_dx;
	
	drho_dxi = drho_mui * dmui_dxi;
	drho_dxj = drho_muj * dmuj_dxj;
	drho_dxdelta = -drho_dxi - drho_dxj;
	
	sig0 = cot(lambda_i) * cot(rho);
	sig1 = cos(lambda_j)*csc(lambda_i)*csc(rho);
	sig2 = sig0-sig1;
	
	eta.rotation = -asin(sig2);
	
	
	deta_dlambdai = (csc(lambda_i) * sig0/cos(lambda_i) - sig1*cos(lambda_i)) / sqrt(1-sig2*sig2);
	deta_dlambdaj = - sig1*tan(lambda_j) / sqrt(1-sig2*sig2);
	deta_drhoij = (csc(rho) * sig0/cos(rho) - sig1*cos(rho)) / sqrt(1-sig2*sig2);
	
	
	deta_dxi = deta_dlambdai * I.lambda.drotation_dxi + deta_drhoij * drho_dxi;
	deta_dxj = deta_dlambdaj * J.lambda.drotation_dxj + deta_drhoij * drho_dxj;
	deta_dxdelta = -deta_dxi-deta_dxj;
	
	eta.drotation_dxi = deta_dxi;
	eta.drotation_dxj = deta_dxj;
	eta.drotation_dxdelta = deta_dxdelta;
	
	return eta;
}



PHIContainer Tessellation::calculatePHI(Vector &tessellationOrigin, CircularInterface &I, CircularInterface J, double dij, double radius){
	double eta;
	double omega;
	PhiContainer p;
	
	Rotation eta,omega;
	
	
	//eta = M_PI/2 - acos(-(1/tan(lambda_j))*(1/tan(rho))+cos(lambda_k)*(1/sin(lambda_j))*(1/sin(rho)));
	
	eta = calculateOmega(I,J);
	omega = calculateOmega(tessellationOrigin,I,J);
	

	
	p.out.rotation = eta.rotation + omega.rotation;
	if(p.out.rotation>=M_PI) p.out.rotation = (eta.rotation + omega.rotation) - 2*M_PI;
	p.out.drotation_dxi = eta.drotation_dxi + omega.drotation_dxi;
	p.out.drotation_dxj = eta.drotation_dxj + omega.drotation_dxj;
	p.out.drotation_dxdelta = eta.drotation_dxdelta + omega.drotation_dxdelta;
	
	
	
	p.in.rotation = (-M_PI-eta.rotation) + omega.rotation;
	if(p.in.rotation<=-M_PI) p.in.rotation = (M_PI-eta.rotation) + omega.rotation;
	p.in.drotation_dxi = -eta.drotation_dxi + omega.drotation_dxi;
	p.in.drotation_dxj = -eta.drotation_dxj + omega.drotation_dxj;
	p.in.drotation_dxdelta = -eta.drotation_dxdelta + omega.drotation_dxdelta;
			    
	
	return p;
	
	
	
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
	it0 = J.intersectionBranches.insert(x);
	
	//printf("address %d-%d side %d IN\n",address.id0, address.id1, address.id1);
	
	x.first = PHII.out.rotation;
	x.second.node=&intersectionGraph[address];
	x.second.direction = OUT;
	x.second.it = it0;
	x.second.body = &I.intersectionBranches;
	x.second.id = address.id0;
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


void Tessellation::buildIntersectionGraph(double radius, Vector &tessellationOrigin, CircularInterfacesPerAtom &circles, SASAs &sasas, Hemisphere hemisphere, string filename){
	map<int, bool> processed;
	map<int, bool>::iterator it_p;
	IntersectionGraph intersectionGraph;
	IntersectionGraph::iterator it_g, it_x;
	int cid0;
	
	SASA potentialSasa;
	SASANode sasaNode;
	
	int i,j;
	CircularInterface *I, *J;
	double tau0, tau1;
	map<int,CircularIntersection>::iterator it_j;
	
	bool empty;
	IntersectionAddress start, x, t;
	bool valid;
	map<IntersectionBranches::iterator,bool,IteratorComparator> eraseList;
	map<IntersectionBranches::iterator,bool,IteratorComparator>::iterator it_e;
	
	IntersectionPair ip;
	
	
	PHIContainer interfacesJ, interfacesI;
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
		
		//we don't need to calculate areas, if it's the splitter interface
		if(I->circularIntersections.size()==0 && I->form!=SPLITTER){
			insertArtificialIntersectionPoints(*I,intersectionGraph,tessellationOrigin);			
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
				
				
				//retrieve external and internal PHI values (respectively)
				PHIJ = calculatePHI(tessellationOrigin, *J, *I, it_j->second.d, radius);
				PHII = calculatePHI(tessellationOrigin, *I, *J, it_j->second.d, radius);
				
				createIntersectionBranch(addressIJ, PHII, PHIJ, *I, *J, intersectionGraph);
				createIntersectionBranch(addressJI, PHIJ, PHII, *J, *I, intersectionGraph);
				
			}
		}
		
		
		/*
		for(int m=0; m<circles.size(); ++m){
			printf("circle %d\n",m);
			for(it_main = circles[m].intersectionBranches.begin(); it_main != circles[m].intersectionBranches.end(); ++it_main){
				if(it_main->second.direction==OUT)
					printf("[%d] (%d,%d) d: %f - (%f,%f) (OUT)\n",it_main->second.id, it_main->second.node->id0, it_main->second.node->id1, it_main->first, it_main->second.node->angle0, it_main->second.node->angle1);
				else
					printf("[%d] (%d,%d) d: %f - (%f,%f) (IN)\n",it_main->second.id, it_main->second.node->id0, it_main->second.node->id1, it_main->first, it_main->second.node->angle0, it_main->second.node->angle1);
			}
		}
		*/
			
			
		//all intersectionpoints have been added, it is time to change topologies
		//start at one of I's intersectionbranches
		eraseList.clear();
		for(it_main = I->intersectionBranches.begin(); it_main != I->intersectionBranches.end(); ++it_main){
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
				
				
				//if no intersections have been added to the interface yet, do a simple interface-interface connect
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
		
		for(it_e=eraseList.begin(); it_e != eraseList.end(); ++it_e){
			it = it_e->first;
			deleteIntersectionPoint(it,intersectionGraph);
		}
		
	}
	
	int s = 0;
	for(it_g = intersectionGraph.begin(); it_g != intersectionGraph.end(); ++it_g){
		if(!it_g->second.visited){
			start.id0 = it_g->second.id0;
			start.id1 = it_g->second.id1;
			
			x.id0 = start.id0;
			x.id1 = start.id1;
			
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
				sasaNode.rotation0 = intersectionGraph[x].rotation0;
				sasaNode.rotation1 = intersectionGraph[x].rotation1;
				sasaNode.lambda = circles[cid0].lambda;
				sasaNode.psi = circles[cid0].psi;
				sasaNode.vector = intersectionGraph[x].vector;
				sasaNode.normalForCircularInterface = circles[cid0].normal;
				sasaNode.form = circles[cid0].form;
				
				
				potentialSasa.sasa.push_back(sasaNode);
				
				t.id0 = intersectionGraph[x].pointsTo0;
				t.id1 = intersectionGraph[x].pointsTo1;
				
				x = t;
				
			}
			while(!(intersectionGraph[x].id0 == start.id0 && intersectionGraph[x].id1 == start.id1));

			if(valid){
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


