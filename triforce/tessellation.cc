#include "tessellation.h"





Tesselation::Tesselation(Molecule *m){
	molecule = m;
	
	build();
}

void Tesselation::build(){
	molecule->update();
	atoms = molecule->coordinates();
	int i;
	
	//iterate over all atoms and build the tesselation for each of them
	for(i=0; i<atoms.size(); ++i){
		
		//buildGaussBonnetPath(atoms[i], radius, &atoms, &radii, &circularRegions);
		
	}
	
	

}





double Tesselation::cot(double a){
	return 1.0/tan(a);
}

double Tesselation::csc(double a){
	return 1.0/sin(a);
}




double Tesselation::getAngleBetweenNormals(Vector &a, Vector &b){
	return acos(dot(a,b));
}


double Tesselation::getAngle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}


bool Tesselation::isZero(double v){
	if(fabs(v) <= THRESHOLD_NUMERICAL) return true;
	else return false;
}



void Tesselation::determineProjection(Vector &origin, double radius, CircularRegion &circle){
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
	
	if(g<0) circle.form=CONVEX;
	else circle.form=CONCAVE;

	circle.openingAngle = acos(g / r_i);
	circle.g = g;
	circle.a = a;

	//NSWarnLog(@"DET: [%f %f %f] %f %f %f %f",circle.normal.vector[0],circle.normal.vector[1],circle.normal.vector[2], r_i, r_k, g, circle.openingAngle);


	//NSWarnLog(@"%f %f %f %f %f",d_k, r_i, r_k, g, circle.openingAngle);
}


IntersectionPair Tesselation::determineIntersectionPoints(double radius, CircularRegion &K, CircularRegion &J){
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



void Tesselation::makeCircularRegions(Vector &origin, double radius, std::vector<vec> &atoms, std::vector<int> &indizee, std::vector<double> &radii, std::vector<CircularRegion> &circles){
	int i;
	Vector res,normal,v;
	CircularRegion circle;
	double r_i = radius;
	double r_k;
	int lenv;

	for(i=0;i<atoms.size();i++){
		v=atoms[i] - origin;
		r_k = radii[i];
		
		lenv = norm(v,2);

		//reject, if no intersection
		if(lenv < r_i + r_k && lenv+r_k > r_i){
			normal = v / lenv;
			circle.id = circles.size();
			circle.index=indizee[i];
			circle.vector = v;
			circle.normal = normal;
			circle.sphereRadius = r_k;
			circle.intersect = false;
			circles.push_back(circle);
			
			//NSWarnLog(@"CIRCLE[%d]: (%f, %f, %f) %f",circle.id, circle.vector.vector[0], circle.vector.vector[1], circle.vector.vector[2], circle.sphereRadius);
		}
	}
}



void Tesselation::filterCircularRegions(double radius, std::vector<CircularRegion> &circles){
	std::vector<CircularRegion>::iterator it;
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
				if(it->openingAngle + angle < circles[i].openingAngle){
					it = circles.erase(it);
					//NSWarnLog(@"deleting region");
					erased=true;
					break;
				}
			}
		}
		if(!erased) ++it;
	}

}

void Tesselation::filterEmptyCircularRegions(std::vector<CircularRegion> &circles){
	std::vector<CircularRegion>::iterator it;
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



void Tesselation::reindexCircularRegions(std::vector<CircularRegion> &circles){
	std::list<IntersectionPoint>::iterator it;
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


void Tesselation::outputCircularRegions(std::vector<CircularRegion> &circles){
	std::list<IntersectionPoint>::iterator it;
	int i,u;
	fprintf(stdout, "......... outputting circular regions ...........");
	for(i=0;i<circles.size();i++){
		fprintf(stdout,"- circular region %d with id: %d:",i,circles[i].id);
		u=0;
		for(it = circles[i].forwardIntersections.begin(); it != circles[i].forwardIntersections.end(); ++it){
			fprintf(stdout,"intersection %d from: %d with: %d (%f %f %f)",u,it->from,it->with,it->vector(0),it->vector(1),it->vector(2));
			++u;
		}
	}
	fprintf(stdout,"......... end of data  ...........");

}

/*
void outputSingleCircularRegion(CircularRegion *circle){
	
	NSWarnLog(@"...............................................");
	NSWarnLog(@"circular region [%d / %d]:",circle->id, circle->index);
	NSWarnLog(@"vector: (%f %f %f) normal: (%f %f %f)",circle->vector.vector[0],circle->vector.vector[1],circle->vector.vector[2],circle->normal.vector[0],circle->normal.vector[1],circle->normal.vector[2]);
	NSWarnLog(@"length: %f(%f), g: %f, a: %f, openingAngle: %f, sphereRadius: %f",circle->vector.length,circle->normal.length,circle->g, circle->a, circle->openingAngle, circle->sphereRadius); 
	
}
*/


void Tesselation::filterIntersectionPoints(std::vector<CircularRegion> &circles, int except){
	std::list<IntersectionPoint>::iterator it;
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



void Tesselation::clearFlags(std::vector<CircularRegion> &circles){
	std::list<IntersectionPoint>::iterator it;
	int i;	
	for(i=0;i<circles.size();i++){
		it = circles[i].forwardIntersections.begin();
		for(it = circles[i].forwardIntersections.begin(); it != circles[i].forwardIntersections.end(); ++it){
			it->visited=false;
		}
	}
}


std::vector<CircularRegion>* Tesselation::deepCopy(std::vector<CircularRegion> &circles){
	int i,j;
	std::list<IntersectionPoint>::iterator it;
	std::vector<CircularRegion>* newCircles;
	newCircles = new std::vector<CircularRegion>();
	
	CircularRegion c;
	for(i=0;i<circles.size();i++){
		c.id=circles[i].id;
		c.index=circles[i].index;
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


void Tesselation::outputGaussBonnetPath(std::list<IntersectionPoint*> &points){
	std::list<IntersectionPoint*>::iterator it;
	int i;
	
	for(it=points.begin(), i=0; it!=points.end(); ++it, ++i)
		fprintf(stdout,"GBPATH[%d] %d - %d - %d", i, (*it)->from, (*it)->id, (*it)->with);
}



void Tesselation::prepareCircularRegions(std::vector<CircularRegion> &circles, std::vector<CircularRegion> **newCircles){
	*newCircles = deepCopy(circles);
	
	filterEmptyCircularRegions(**newCircles);
	clearFlags(**newCircles);
	reindexCircularRegions(**newCircles);
	clearFlags(**newCircles);
}


void Tesselation::insertFakeIntersectionPoints(std::vector<CircularRegion> &circles){
	Vector v,v2,o,o2;
	int k;
	IntersectionPoint p;
	
	//we have to find circles that did not intersect with other circles and we have to insert fake intersectionpoints for them
	for(k=0;k<circles.size();k++){
		if(!circles[k].intersect){
			//NSWarnLog(@"INSERTING FAKE INTERSECTIONPOINTS");
			v = randu<vec>(3);
			v2 = circles[k].normal * circles[k].g;
			o = cross(v2,v);
			o = o / norm(o,2);
			p.vector = v2+o;
			p.with=k;
			p.from=k;
			circles[k].forwardIntersections.push_back(p);
			
			o2 = cross(v2,o);
			o2 = o2 / norm(o2,2);
			p.vector=v2+o2;
			circles[k].forwardIntersections.push_back(p);
			
			p.vector = v2 + (o * -1);
			circles[k].forwardIntersections.push_back(p);
			
			p.vector = v2 + (o2 * -1);
			circles[k].forwardIntersections.push_back(p);
			
		}
	}
	
}



void Tesselation::buildGaussBonnetPath(Vector &origin, double radius, std::vector<vec> &atoms, std::vector<int> &indizee, std::vector<double> &radii, std::vector<CircularRegion> &circles){
	IntersectionPoint p;
	int i,k,j;
	IntersectionPair points;
	double angle;
	Vector v,v2,o,o2;
	
	srand(2);


	makeCircularRegions(origin, radius, atoms, indizee, radii, circles);
	for(i=0;i<circles.size();i++){
		determineProjection(origin, radius, circles[i]);
	}
	filterCircularRegions(radius, circles);
	
	clearFlags(circles);
	reindexCircularRegions(circles);

	p.visited = false;
	p.flagged = false;
	for(k=0;k<circles.size();k++)
		for(j=k+1;j<circles.size();j++){
			if(k!=j){
				//determine if there will be intersections
				angle = getAngleBetweenNormals(circles[k].normal,circles[j].normal);
				if(angle < circles[k].openingAngle + circles[j].openingAngle)
					if(angle + circles[k].openingAngle > circles[j].openingAngle && angle + circles[j].openingAngle > circles[k].openingAngle){
						points = determineIntersectionPoints(radius, circles[k], circles[j]);
						//NSWarnLog(@"POINTS: [%f %f %f], [%f %f %f]",points.k_j.vector[0],points.k_j.vector[1],points.k_j.vector[2],points.j_k.vector[0],points.j_k.vector[1],points.j_k.vector[2]);
						p.vector = points.k_j;
						p.with = j;
						p.from = k;
						circles[k].forwardIntersections.push_back(p);	
						circles[k].intersect = true;
						p.vector = points.j_k;
						p.with = k;
						p.from = j;
						circles[j].forwardIntersections.push_back(p);						
						circles[j].intersect = true;
					}
			}
		}

	//NSWarnLog(@"filtering intersectionpoints...");
	filterIntersectionPoints(circles,-1);
	//NSWarnLog(@"filtering circular regions...");

	//we're not filtering empty regions at this point, we need all regions later on
	//filterEmptyCircularRegions(circles);
	//because we're not deleting, we don't need reindexing
	//NSWarnLog(@"reindexing circular regions...");
	//reindexCircularRegions(circles);
	//NSWarnLog(@"clearing flags...");
	clearFlags(circles);
	//NSWarnLog(@"gauss bonnet path is ready...");
	
	insertFakeIntersectionPoints(circles);
	
	
	
}




void Tesselation::harvestIntersectionPoints(std::vector<CircularRegion> &circles, std::vector<vec> &intersections){
	int i;
	std::list<IntersectionPoint>::iterator it;

	for(i=0;i<circles.size();i++)
		for(it = circles[i].forwardIntersections.begin(); it != circles[i].forwardIntersections.end(); ++it){
			//NSWarnLog(@"POINT: [%f %f %f]",it->vector.vector[0],it->vector.vector[1],it->vector.vector[2]);
			intersections.push_back(it->vector);
			}
}


bool Tesselation::hasUnflaggedIntersectionPoints(CircularRegion &circle, IntersectionPoint **ip){
	std::list<IntersectionPoint>::iterator it;
	for(it = circle.forwardIntersections.begin(); it != circle.forwardIntersections.end(); ++it){
		if(!it->visited){
			*ip=&*it;
			return true;
		}
	}
	return false;
}


std::list<IntersectionPoint*>* Tesselation::retrieveIntersections(CircularRegion &circle){
	std::list<IntersectionPoint*>* res;
	res = new std::list<IntersectionPoint*>();
	
	std::list<IntersectionPoint>::iterator it;
	for(it = circle.forwardIntersections.begin(); it != circle.forwardIntersections.end(); ++it){
		//if(!it->visited){
			//save a pointer to the intersectionpoint
			res->push_back(&(*it));
		//}
	}
	return res;
	
}


void Tesselation::showIntersections(std::list<IntersectionPoint*> &intersections){
	std::list<IntersectionPoint*>::iterator it;
	
	for(it=intersections.begin();it!=intersections.end();++it)
		fprintf(stdout,"INTERSECTION %d-%d",(*it)->from, (*it)->with);
}

/*
 * path and res needs to be deleted somewhere
 */
std::vector<std::list<IntersectionPoint*>*>*  Tesselation::harvestGaussBonnetPaths(std::vector<CircularRegion> &circles){
	std::vector<std::list<IntersectionPoint*>*>* res;
	std::list<IntersectionPoint*>* path;
	std::list<IntersectionPoint*>* intersectionsTest;
	std::list<IntersectionPoint*>* intersections;
	std::list<IntersectionPoint*>::iterator it_i, minIt;

	Vector v,ortho,test,reference;
	IntersectionPoint  *p, *p_prev;
	CircularRegion circle;
	bool done;
	double minAngle,angle;
	int i,j;
	
	res = new std::vector<std::list<IntersectionPoint*>*>();
	path = new std::list<IntersectionPoint*>();
	//the whole thing can actually be emtpy, in that case we better get outta here
	if(circles.size()==0) return res;
	//outputCircularRegions(circles);
	
	clearFlags(circles);
	
	
	for(j=0;j<circles.size();++j){
		if(hasUnflaggedIntersectionPoints(circles[j],&p_prev)){
			//NSWarnLog(@"USING SEED %d, %d-%d",j, p_prev->from, p_prev->with);
			//p_prev->visited=true;
			done=false;
			while(!done){
				//NSWarnLog(@"STARTING ROUND %d",p_prev->with);
				circle = circles[p_prev->with];
				//NSWarnLog(@"GETTING INTERSECTIONS FOR: %d",p_prev->with);
				intersections = retrieveIntersections(circles[p_prev->with]);
				
				if(intersections->size() > 1){
					showIntersections(*intersections);
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
						if(circle.form == CONCAVE){
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
			path = new std::list<IntersectionPoint*>();
		}
		delete intersectionsTest;
	}
	return res;
}




/*
void splitCircularRegions(Vector3D *origin, double radius, std::vector<CircularRegion> *circles, CircularRegion *splitter, std::vector<CircularRegion> **inCircles, std::vector<CircularRegion> **offCircles){
	//insert new points
	CircularRegion c;
	int i;
	double angle;
	IntersectionPair points;
	IntersectionPoint p;
	std::list<IntersectionPoint>::iterator it;
	//copy
	*inCircles = deepCopy(circles);
	*offCircles = deepCopy(circles);
	
	for(i=0;i<circles->size();i++){
		(**inCircles)[i].forwardIntersections.clear();
		(**offCircles)[i].forwardIntersections.clear();
	}
	splitter->id = circles->size();
	//NSWarnLog(@"SPLITTER ID: %d",splitter.id);
	splitter->form = CONVEX;
	(*inCircles)->push_back(*splitter);
	splitter->form = CONCAVE;
	(*offCircles)->push_back(*splitter);
	
	//outputCircularRegions(inCircles);
	
	
	p.flagged=false;

	for(i=0;i<circles->size();i++){
		for(it = (*circles)[i].forwardIntersections.begin(); it != (*circles)[i].forwardIntersections.end(); ++it){
			angle = getAngle(&it->vector,&splitter->normal);
			//NSWarnLog(@"COMPARE: %f %f",splitter.openingAngle, angle);
			if(splitter->openingAngle > angle){
				p = *it;
				(**inCircles)[i].forwardIntersections.push_back(p);
				//NSWarnLog(@"%d %d %f INCIRCLE %d",i,it->with,angle,p.with);
			}
			else{
				p = *it;
				(**offCircles)[i].forwardIntersections.push_back(p);
				//NSWarnLog(@"%d %d %f OFFCIRCLE %d",i,it->with,angle,p.with);
			}
			
		}
		


		//determine if there will be intersections with the splitter
		angle = getAngleBetweenNormals(&(*circles)[i].normal,&splitter->normal);
		
		//NSWarnLog(@"compare: (%f %f %f), (%f %f %f) a:%f oi: %f os: %f", splitter.normal.vector[0], splitter.normal.vector[1], splitter.normal.vector[2], circles[i].normal.vector[0], circles[i].normal.vector[1], circles[i].normal.vector[2],angle,circles[i].openingAngle,splitter.openingAngle);
		//NSWarnLog(@"radius: %f", radius);
		//outputSingleCircularRegion(splitter);
		//outputSingleCircularRegion(circles[i]);
		if(angle < (*circles)[i].openingAngle + splitter->openingAngle)
			if(angle + (*circles)[i].openingAngle > splitter->openingAngle && angle + splitter->openingAngle > (*circles)[i].openingAngle){
				points = determineIntersectionPoints(origin, radius, &((*circles)[i]), splitter);
				//NSWarnLog(@"POINTS: [%f %f %f], [%f %f %f]",points.k_j.vector[0],points.k_j.vector[1],points.k_j.vector[2],points.j_k.vector[0],points.j_k.vector[1],points.j_k.vector[2]);
				p.vector = points.k_j;
				p.with = splitter->id;
				p.from = (*circles)[i].id;
				(**offCircles)[i].forwardIntersections.push_back(p);
				
				p.vector = points.j_k;
				p.with = splitter->id;
				p.from = (*circles)[i].id;
				(**inCircles)[i].forwardIntersections.push_back(p);						
								
				p.vector = points.j_k;
				p.with = (*circles)[i].id;
				p.from = splitter->id;
				(**offCircles)[splitter->id].forwardIntersections.push_back(p);
				
				p.vector = points.k_j;
				p.with = (*circles)[i].id;
				p.from = splitter->id;
				(**inCircles)[splitter->id].forwardIntersections.push_back(p);						
				
				//NSWarnLog(@"SPLITTING INTERSECTION (%f %f %f) (%f %f %f)",points.k_j.vector[0],points.k_j.vector[1],points.k_j.vector[2],points.j_k.vector[0],points.j_k.vector[1],points.j_k.vector[2]);
			}
	}
	
	//outputCircularRegions(inCircles);
	
	filterIntersectionPoints(origin, radius, *inCircles,splitter->id);
	filterIntersectionPoints(origin, radius, *offCircles,-1);
	//outputCircularRegions(inCircles);


	filterEmptyCircularRegions(*inCircles);
	filterEmptyCircularRegions(*offCircles);
	//outputCircularRegions(inCircles);

	clearFlags(*inCircles);
	clearFlags(*offCircles);

	reindexCircularRegions(*inCircles);
	reindexCircularRegions(*offCircles);

	clearFlags(*inCircles);
	clearFlags(*offCircles);
	//outputCircularRegions(inCircles);
	
		
	
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
	if(g<0) res.form=CONVEX;
	else res.form=CONCAVE;
	//this particular sphereradius is just a construct for later calculations and does represent anything physical 
	res.sphereRadius = sqrt((K->length-g)*(K->length-g) + a * a);

	return res;
}
*/
