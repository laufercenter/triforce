#include "multiLayeredDepthBuffer.h"




MultiLayeredDepthBuffer::MultiLayeredDepthBuffer(Depth3D &data, Data1D &occludedDistribution, Data1D &exposedDistribution, int mode){
	this->data=data;
	this->len = data.parameter0Dim;
	this->mode=mode;
	//create buffer
	dbuffer.resize(len, DepthBufferLine());
	dmode.resize(len, SCANLINE_EMPTY);
	orientationAxis=Vector(3);
	orientationAxisAuxiliary=Vector(3);
	orientationPlaneNormal=Vector(3);
	
	this->occludedDistribution = occludedDistribution;
	this->exposedDistribution = exposedDistribution;
	
	
	if(mode==0){
		orientationAxis(0)=1;
		orientationAxis(1)=0;
		orientationAxis(2)=0;
		
		orientationAxisAuxiliary(0)=0;
		orientationAxisAuxiliary(1)=1;
		orientationAxisAuxiliary(2)=0;
		
		orientationPlaneNormal(0)=0;
		orientationPlaneNormal(1)=0;
		orientationPlaneNormal(2)=1;
	}
	else{
		orientationAxis(0)=0;
		orientationAxis(1)=1;
		orientationAxis(2)=0;
		
		orientationAxisAuxiliary(0)=0;
		orientationAxisAuxiliary(1)=0;
		orientationAxisAuxiliary(2)=1;
		
		orientationPlaneNormal(0)=1;
		orientationPlaneNormal(1)=0;
		orientationPlaneNormal(2)=0;
	}
	
test=false;
}




DepthBufferLine::iterator MultiLayeredDepthBuffer::increaseLineInterator(DepthBufferLine::iterator it, DepthBufferLine &line){
	++it;
	if(it == line.end()) it=line.begin();
		
	
	return it;
	
}

DepthBufferLine::iterator MultiLayeredDepthBuffer::decreaseLineInterator(DepthBufferLine::iterator it, DepthBufferLine &line){	
	if(it == line.begin()) it=line.end();
	--it;
		
	
	return it;
}




ScanlineMode MultiLayeredDepthBuffer::insertIntoLineBuffer(DepthBufferLine &line, double front, double back){
	DepthBufferLine::iterator it0, it1, it, it_next;
	bool deleteBack, deleteFront;
	pair<double,LineType> p0, p1;
	pair<DepthBufferLine::iterator, bool> r0,r1;
	
	
	
	p0.first=front;
	p0.second=FRONT;
	p1.first=back;
	p1.second=BACK;
	
	r0 = line.insert(p0);
	r1 = line.insert(p1);
	
	it0=r0.first;
	it1=r1.first;
	
	
	//prune linebuffer
	it=decreaseLineInterator(it0,line);
	if(it->second==FRONT) deleteFront=true;
	it=increaseLineInterator(it0,line);
	if(it->second==BACK && it!=it1) deleteFront=true;
	
	it=increaseLineInterator(it1,line);
	if(it->second==BACK) deleteBack=true;
	it=decreaseLineInterator(it1,line);
	if(it->second==FRONT && it!=it0) deleteBack=true;
	
	
	it=increaseLineInterator(it0,line);
	while(it!=it1){
		it_next=increaseLineInterator(it,line);
		//printf("deleting %f-%d\n",it->first,it->second);
		line.erase(it);
		it=it_next;
	}
	
	if(deleteFront){
		//printf("deleting %f-%d\n",it0->first,it0->second);
		line.erase(it0);
	}
	if(deleteBack){
		//printf("deleting %f-%d\n",it1->first,it1->second);
		line.erase(it1);
	}
	
	if(line.size()==0) return SCANLINE_FULL;
	else return SCANLINE_PARTIAL;
		
}




bool MultiLayeredDepthBuffer::wouldChangeLineBuffer(DepthBufferLine &line, double front, double back, bool &frontChanges, bool &backChanges){
	DepthBufferLine::iterator it0, it1, it, it_next;
	bool wouldChange;
	pair<double,LineType> p0, p1;
	pair<DepthBufferLine::iterator, bool> r0,r1;
	
	p0.first=front;
	p0.second=FRONT;
	p1.first=back;
	p1.second=BACK;
	
	r0 = line.insert(p0);
	r1 = line.insert(p1);
	
	it0=r0.first;
	it1=r1.first;
	
	wouldChange=false;
	frontChanges=false;
	backChanges=false;
// 	printf("CHECK: %f %f\n",front,back);
	
	if(line.size()==0){
		wouldChange=true;
		frontChanges=true;
		backChanges=true;
	}
	
	if(!r0.second){
		wouldChange=true;
		frontChanges=true;
// 		printf("CHECK 0\n");
	}
	if(!r1.second){
		wouldChange=true;
		backChanges=true;
// 		printf("CHECK 1\n");
	}
	
	it=decreaseLineInterator(it0,line);
	if(it->second==BACK && it!=it1){
		wouldChange=true;
		frontChanges=true;
// 		printf("CHECK 2 %d %f\n",it->second,it->first);
	}
	it=increaseLineInterator(it0,line);
	if(it->second==FRONT && it!=it1){
		wouldChange=true;
		frontChanges=true;
// 		printf("CHECK 3 %d %f\n",it->second,it->first);
	}
	
	it=increaseLineInterator(it1,line);
	if(it->second==FRONT && it!=it0){
		wouldChange=true;
		backChanges=true;
// 		printf("CHECK 4 %d %f\n",it->second,it->first);
	}
	it=decreaseLineInterator(it1,line);
	if(it->second==BACK && it!=it0){
		wouldChange=true;
		backChanges=true;
// 		printf("CHECK 5 %d %f\n",it->second,it->first);
	}
	
	it=increaseLineInterator(it0,line);
	while(it!=it1){
		it_next=increaseLineInterator(it,line);
		wouldChange=true;
// 		frontChanges=true;
// 		backChanges=true;
		it=it_next;
	}
	

	if(r0.second) line.erase(it0);
	if(r1.second) line.erase(it1);
	
	return wouldChange;
		
}



Vector MultiLayeredDepthBuffer::normalise(Vector x){
	double l;
	Vector v(3);
	l = norm(x,2);
	v = x/l;
	return v;
	
}





int MultiLayeredDepthBuffer::sgn(double d){
	if(d>=0) return 1;
	else return -1;
}

void MultiLayeredDepthBuffer::addSphere(Vector &v, double lambda, bool invert, double &kappa, double &psi, int index){
	Vector c;
	double x;
	double front,back;
	DepthInformation dat;
	Vector n,n2;
	double s;
	
	
	
	
	psi=acos(dot(v,orientationAxis));
	
	
	
	
	//construct normal
	n2=cross(v, orientationAxis);
	n=normalise(cross(orientationAxis,n2));
	s = sgn(dot(n,orientationPlaneNormal));
	
	kappa = acos(dot(n,orientationAxisAuxiliary));
	if(s<0) kappa=2*M_PI-kappa;
	
	
	//printf("ADDING SPHERE %d: lambda:%f psi:%f kappa:%f cart(%f %f %f) form: %d\n",index,lambda, psi, kappa, v(0), v(1), v(2), invert);
	//printf("n2: (%f %f %f)\n", n2(0), n2(1), n2(2));
	//printf("n: (%f %f %f)\n", n(0), n(1), n(2));

	dat = data.getFloorScanlines(kappa, psi, lambda, invert);
	
	for(int i=0; i<len; ++i){
		//printf("line[%d]...",i);
		
		if(dmode[i]!=SCANLINE_FULL){
			if(dat.mode[i]==SCANLINE_FULL){
				dmode[i]=SCANLINE_FULL;
				dbuffer[i].clear();
				//printf("F");
			}
			else if(dat.mode[i]==SCANLINE_PARTIAL){
				dmode[i]=insertIntoLineBuffer(dbuffer[i], dat.scanline0[i], dat.scanline1[i]);
				//printf("%f %f",dat.scanline0[i], dat.scanline1[i]);
				
			}
			//else printf("E");
		}
		else{
			//printf("#");
			
		}
		//printf("\n");
	}
	
	//print();

	
}


bool MultiLayeredDepthBuffer::passesBuffer(double kappa, double psi, double lambda, bool invert, int index, vector<Vector> &exposedVectors){
	DepthInformation dat;
	bool pass;
	dat = data.getCeilScanlines(kappa, psi, lambda, invert);
	pass=scanBuffer(kappa, index, exposedVectors, dat);
	if(invert) pass=true;
	return pass;
}


bool MultiLayeredDepthBuffer::getSplitterExposedVectors(vector<Vector> &exposedVectors){
	DepthInformation dat;
	bool res;
	dat.scanline0.resize(len,0);
	dat.scanline1.resize(len,M_PI);
	dat.mode.resize(len,SCANLINE_PARTIAL);
	res=scanBuffer(M_PI, -1, exposedVectors, dat);
	return res;
}


bool MultiLayeredDepthBuffer::scanBuffer(double kappa, int index, vector<Vector> &exposedVectors, DepthInformation &dat){
	Vector c;
	double x;
	double front,back;
	double s;
	double offset=0;
	bool wouldChange;
	bool frontChanges, backChanges;
	DepthBufferLine::iterator it;
	vector<DepthBufferProjection> projections[2];
	vector<DepthBufferProjection> segment;
	vector<vector<DepthBufferProjection> > segments;
	DepthBufferProjection pr;
	double g_prev, k_prev;
	
	if(index==2) test=true;
	
	int ic=524;
	
	//if(invert) return true;
	
	wouldChange=false;
	
	for(int i=0; i<len; ++i){
// 		printf("line %d %d (%f %f)\n",i,dat.mode[i], dat.scanline0[i]-offset, dat.scanline1[i]+offset);
		if(index==ic){
			printf("line[%d]: ",i);
		}
		if(dmode[i]==SCANLINE_PARTIAL){
			if(dat.mode[i]==SCANLINE_FULL){
				wouldChange=true;
				dmode[i]=SCANLINE_FULL;
				for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
					//exposedVectors.push_back(convertToCartesian(data.getHeaderParameter0Cell(i), it->first));
					pr.g = data.getHeaderParameter0Cell(i);
					pr.kappa = it->first;
					printf("PR: %d\n",it->second);
					projections[it->second].push_back(pr);
				}
				
				
				if(index==ic){
					printf("F:");
				}
				
			}
			else if(dat.mode[i]==SCANLINE_PARTIAL){
				if(index==ic){
					printf("P: %f %f",dat.scanline0[i], dat.scanline1[i]);
				}
				
				if(wouldChangeLineBuffer(dbuffer[i], dat.scanline0[i]-offset, dat.scanline1[i]+offset, frontChanges, backChanges)){
					wouldChange=true;
					if(frontChanges){
						//exposedVectors.push_back(convertToCartesian(data.getHeaderParameter0Cell(i),dat.scanline0[i]-offset));
						pr.g = data.getHeaderParameter0Cell(i);
						pr.kappa = dat.scanline0[i]-offset;
						projections[0].push_back(pr);
						
					}
					if(backChanges){
						//exposedVectors.push_back(convertToCartesian(data.getHeaderParameter0Cell(i),dat.scanline1[i]+offset));
						pr.g = data.getHeaderParameter0Cell(i);
						pr.kappa = dat.scanline1[i]+offset;
						projections[1].push_back(pr);
					}
					
				}
			}
			else{
				if(index==ic){
					printf("E: ");
				}
			}
		}
		else if(dmode[i]==SCANLINE_EMPTY){
			if(dat.mode[i]!=SCANLINE_EMPTY){
				if(dat.mode[i]!=SCANLINE_FULL){
					//exposedVectors.push_back(convertToCartesian(data.getHeaderParameter0Cell(i),dat.scanline0[i]-offset));
					//exposedVectors.push_back(convertToCartesian(data.getHeaderParameter0Cell(i),dat.scanline1[i]+offset));
					pr.g = data.getHeaderParameter0Cell(i);
					pr.kappa = dat.scanline0[i]-offset;
					projections[0].push_back(pr);
						
					//exposedVectors.push_back(convertToCartesian(data.getHeaderParameter0Cell(i),dat.scanline1[i]+offset));
					pr.g = data.getHeaderParameter0Cell(i);
					pr.kappa = dat.scanline1[i]+offset;
					projections[1].push_back(pr);
					
				}
				wouldChange=true;
			}
		}
		else{
			if(dat.mode[i]==SCANLINE_FULL){
			}
			else{
			}
			
			if(index==ic){
				if(dat.mode[i]==SCANLINE_FULL){
					printf("+: ");
				}
				else{
					printf("#:");
				}
			}
		}
		if(index==ic){
			printf("\n");
		}
		
	}
	test=false;
	
	
	//divide the list into separate segment candidates
	for(unsigned int z=0; z<=1; ++z){
		segment.clear();
		g_prev=-1;
		k_prev=-1;
		for(unsigned int i=0; i<projections[z].size(); ++i){
			if(i!=0){
				if(abs(projections[z][i].g-g_prev)>1){
					segments.push_back(segment);
					segment.clear();
				}
			}
				
			g_prev=projections[z][i].g;
			k_prev=projections[z][i].kappa;
			segment.push_back(projections[z][i]);
			
		}
		if(segment.size()>0) segments.push_back(segment);
	}
	
	//now reduce projections to just one per segment
	for(unsigned int i=0; i<segments.size(); ++i){
		exposedVectors.push_back(convertToCartesian(segments[i][segments[i].size()/2]));
	}
	
	
	
	return wouldChange;
}



void MultiLayeredDepthBuffer::print(){
	DepthBufferLine::iterator it;
	vector<int> line;
	int pos;
	printf("----------------------------------------------------------------\n");
	
	for(int i=0; i<len; ++i){
		printf("[%d]:\t",i);
		for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
			printf("(%f %d) ",it->first,it->second);
		}
		printf("\n");
	}
	
	
	for(int i=0; i<len; ++i){
		printf("[%d]:\t",i);
		line.clear();
		line.resize(len+1,-1);
		for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
			pos=it->first*(double)len/(2.0*M_PI);
			line[pos]+=it->second+1;
		}
		
		for(int j=0; j<len+1; ++j){
			if(dmode[i]==SCANLINE_FULL) printf("#");
			else if(line[j]>=0) printf("%d",line[j]);
			else printf("*");
		}
		printf("\n");
	}
	printf("----------------------------------------------------------------\n");
		
}

void MultiLayeredDepthBuffer::getDepthBounds(int p, double d, DepthBufferLine::iterator &it0, DepthBufferLine::iterator &it1){
	it1 = dbuffer[p].upper_bound(d);
	if(it1==dbuffer[p].end()) it1 = dbuffer[p].begin();
	it0=decreaseLineInterator(it1,dbuffer[p]);
	
}



double MultiLayeredDepthBuffer::probe(Vector &projection, bool &isExtended){
	DepthBufferLine::iterator it0, it1;
	int p;
	double d;
	double d0,d1;
	double r0,r1;
	double dist;
	double m;
	p = data.closestGridPoint(projection(0));
	d = projection(1);
	
	isExtended=dmode[p]==SCANLINE_EXTENDED;
	
	
	//printf("PROJECTION: %f (%d), %f\n",projection(0),p,d);
	
	if(dmode[p]==SCANLINE_EMPTY){
		//printf("scanline empty\n");
		return 0;
	}
	if(dbuffer[p].size()==0){
		//printf("scanline of size 0\n");
		return M_PI;
	}
	
	getDepthBounds(p, d, it0, it1);
	
	if(it0->second==BACK){
		//printf("in exposed area\n");
		//return 0;
	}
	
	
	
	d0=it0->first;
	d1=it1->first;
	
	r0=min(min(abs(d-d0), abs(d-(d0-2*M_PI))), abs(d-(d0+2*M_PI)));
	r1=min(min(abs(d-d1), abs(d-(d1-2*M_PI))), abs(d-(d1+2*M_PI)));
	
	
	//printf("SUB %f, %f %f\n",d, it0->first, it1->first);
	//printf("SUB2 %f %f %f\n", abs(d-d0), abs(d-(d0-2*M_PI)), abs(d-(d0+2*M_PI)));
	//printf("SUB3 %f %f %f\n", abs(d-d1), abs(d-(d1-2*M_PI)), abs(d-(d1+2*M_PI)));
	
	m = min(r0,r1);
	
	if(dmode[p]==SCANLINE_EXTENDED){
		m = sqrt(m*m + extensionAngle[p]);
	}
	
	return m;
}


bool MultiLayeredDepthBuffer::isWithinNumericalLimits(double x, double l){
	if(abs(x-l)<=0.0001) return true;
	else return false;
}


Vector MultiLayeredDepthBuffer::project(bool splitter, bool hemflip, Vector n, double kappa, double PHI, double psi, double lambda, bool invert, Vector &projection0, Vector &projection1){
	Vector v(3),v2(3);
	double g,g2;
	double d,d2,d3;
	Vector n2,n3;
	int s;
	Vector ex(3),ey(3),ez(3);
	
	
	ex(0)=1;
	ex(1)=0;
	ex(2)=0;
	
	ey(0)=0;
	ey(1)=1;
	ey(2)=0;
	
	ez(0)=0;
	ez(1)=0;
	ez(2)=1;
	
	PHI=-PHI;
	
	
	//printf("BURY: %d %d %f %f %f %f\n",hemflip, invert, kappa, PHI, psi, lambda);
	
	
	//if(invert) PHI=-PHI;
	if(hemflip){
		if(splitter){
			PHI=-PHI;
		}
		else{
			if(invert){
				PHI=-(M_PI-PHI);
			}
			else{
				PHI=-(M_PI-PHI);
				
			}
		}
	}
		
	
	

	//calculated projected coordinates
	v(0) = cos(lambda)*cos(psi) + cos(PHI)*sin(lambda)*sin(psi);
	v(1) = sin(kappa)*sin(lambda)*sin(PHI) + cos(kappa) * (-cos(PHI)*cos(psi)*sin(lambda) + cos(lambda)*sin(psi));
	v(2) = -cos(kappa)*sin(lambda)*sin(PHI) + sin(kappa) * (-cos(PHI)*cos(psi)*sin(lambda) + cos(lambda)*sin(psi));
	
	//printf("PIP: (%f %f %f)\n",v(0),v(1),v(2));
	
	g = v(0);
	//d = asin(v(2));
	//if(d<0) d+=2*M_PI;
	
	
	//construct normal
	n2=cross(v, ex);
	n3=normalise(cross(ex,n2));
	s = sgn(dot(n3,ez));
	
	d = acos(dot(n3,ey));
	if(s<0) d=2*M_PI-d;
	
	
	projection0=Vector(2);
	projection0(0) = g;
	projection0(1) = d;
	
	
	
	//v2(0) = cos(lambda)*cos(psi) + cos(PHI)*sin(lambda)*sin(psi);
	//v2(1) = -cos(PHI)*cos(psi)*sin(lambda) + cos(lambda)*sin(psi);
	//v2(2) = -sin(lambda)*sin(PHI);
	
	
	
	//construct normal
	n2=cross(v, ey);
	n3=normalise(cross(ey,n2));
	s = sgn(dot(n3,ex));
	
	d2 = acos(dot(n3,ez));
	if(s<0) d2=2*M_PI-d;
	
	
	
	g2 = v(2);
	//d2 = asin(v(0));
	//if(d2<0) d2+=2*M_PI;
	
	projection1=Vector(2);
	projection1(0) = g2;
	projection1(1) = d2;
	
	
	return v;
	
}


void MultiLayeredDepthBuffer::extend(){
	extensionDistance.resize(len,len);
	extensionAngle.resize(len,0);
	DepthBufferLine line;
	int c;
	int upperBound, lowerBound;
	bool hasSASA;
	double step;
	
	upperBound = len-1;
	lowerBound = 0;
	hasSASA=false;
	
	//first, look for upper and lower line bounds
	for(int i=0; i<len; ++i){
		if(dmode[i]==SCANLINE_PARTIAL){
			hasSASA=true;
			upperBound=min(i,upperBound);
			lowerBound=max(i,lowerBound);
		}
		if(dmode[i]==SCANLINE_EMPTY){
			hasSASA=true;
			upperBound=min(i,upperBound);
			lowerBound=max(i,lowerBound);
		}
	}
	
	//printf("BOUNDS: %d %d\n",upperBound, lowerBound);
	
	if(!hasSASA) return;
	
	
	//forward
	line=dbuffer[lowerBound];
	
	c=len-lowerBound;
	for(int i=0; i<len; ++i){
		if(dmode[i]==SCANLINE_PARTIAL){
			c=0;
			line = dbuffer[i];
		}
		if(dmode[i]==SCANLINE_FULL){
			c++;
			dbuffer[i]=line;
		}
		extensionDistance[i]=c;
	}
	
	

	//backward
	line=dbuffer[upperBound];
	c=upperBound;
	for(int i=len-1; i>=0; --i){
		if(dmode[i]==SCANLINE_PARTIAL){
			c=0;
			line = dbuffer[i];
		}
		if(dmode[i]==SCANLINE_FULL){
			dmode[i]=SCANLINE_EXTENDED;
			c++;
			if(extensionDistance[i]>c){
				dbuffer[i]=line;
				extensionDistance[i]=c;
			}
		}
	}
	

	step=M_PI/len;
	
	for(int i=0; i<len; ++i){
		extensionAngle[i] = (extensionDistance[i]*step)*(extensionDistance[i]*step);
	}
	

	
	
	
}


bool MultiLayeredDepthBuffer::retrieveLimitingInterfaces(vector<LimitingInterface> &limitList){
	bool limitFound;
	int p;
	double g;
	bool hasSASA;
	LimitingInterface limit;
	
	
	if(mode==0) limit.psi=0;
	else limit.psi=M_PI/2;
	
	limitList.clear();

	
	//forward search
	limitFound=false;
	p=0;
	for(int i=0; i<len && !limitFound; ++i){
		if(dmode[i]!=SCANLINE_FULL){
			limitFound=true;
			p=i-1;
		}
	}
	if(limitFound){
		hasSASA=true;
		if(p>=0){
			g=data.getHeaderParameter0Cell(p);
			if(g>0){
				limit.g = g;
				limit.v = orientationAxis;
				limit.flip=true;
			}
			else{
				limit.g = -g;
				limit.v = -orientationAxis;
				limit.flip=false;
			}
			limitList.push_back(limit);
		}
		
		
		
		limitFound=false;
		p=0;
		for(int i=len-1; i>=0 && !limitFound; --i){
			if(dmode[i]!=SCANLINE_FULL){
				limitFound=true;
				p=i+1;
			}
		}
		if(limitFound){
			if(p<len){
				g=data.getHeaderParameter0Cell(p);
				if(g>0){
					limit.g = g;
					limit.v = orientationAxis;
					limit.flip=false;
				}
				else{
					limit.g = -g;
					limit.v = -orientationAxis;
					limit.flip=true;
				}
				limitList.push_back(limit);
			}
			
			
		}
	}
	else hasSASA=false;
	
	return hasSASA;
	
}






void MultiLayeredDepthBuffer::startNewCycle(){
	exposed=1;
	occluded=1;
	prior_occluded=occludedDistribution.getAuxiliaryCell(0);
	prior_exposed=exposedDistribution.getAuxiliaryCell(0);
	printf("priors: %f %f\n",prior_exposed, prior_occluded);
}

bool MultiLayeredDepthBuffer::isExposed(double x){
	unsigned int p;
	double l;
	
	exposedDistribution.closestGridPoint(x,p,l);
	
	printf("adding %f [%d]\n",x,p);
	
	if(exposedDistribution.getDataCell(p)>occludedDistribution.getDataCell(p)) return true;
	else return false;
	
}


void MultiLayeredDepthBuffer::addProbe(double x){
	unsigned int p;
	double l;
	
	exposedDistribution.closestGridPoint(x,p,l);
	
	printf("adding %f [%d]\n",x,p);
	
	exposed*=exposedDistribution.getDataCell(p);
	occluded*=occludedDistribution.getDataCell(p);
	
	printf("prob: %f %f\n",exposed, occluded);
}

bool MultiLayeredDepthBuffer::isCycleExposed(){
	double e;
	e=exposedProbability();
	printf("endprob: %f %f\n",e, 1-e);
	if(e > (1-e)) return true;
	else return false;
}


double MultiLayeredDepthBuffer::exposedProbability(){
	double pe, po;
	pe = exposed * prior_exposed;
	po = occluded * prior_occluded;
	if(pe+po == 0) return 0;
	
	return pe/(pe+po);
}




Vector MultiLayeredDepthBuffer::convertToCartesian(DepthBufferProjection &pr){
	Vector v(3);
	double h;
	double g;
	double kappa;
	
	g=pr.g;
	kappa=pr.kappa;
	
	h=sqrt(1-g*g);
	if(mode==0){
		v(0)=g;
		v(1)=h*cos(kappa);
		v(2)=h*sin(kappa);
	}
	else{
		v(0)=h*sin(kappa);
		v(1)=g;
		v(2)=h*cos(kappa);
	}
	
//  	printf("[%d] g: %f kappa: %f\n",mode, g, kappa);
	//if(test) printf("VECTOR %f %f %f\n", v(0), v(1), v(2));
	return v;
		
}





