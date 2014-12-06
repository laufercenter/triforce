/**
	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
	Email: nils.drechsel@gmail.com
	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
*/

#include "multiLayeredDepthBuffer.h"


#define doublepi 6.283185307179586231996
#define halfpi 1.570796326794896557999


MultiLayeredDepthBuffer::MultiLayeredDepthBuffer(){
}


MultiLayeredDepthBuffer::MultiLayeredDepthBuffer(Depth3D &data, Data1D &occludedDistribution, Data1D &exposedDistribution, int mode){
	float g;
	
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
	
	//init g table
	
	gTable.resize(len,0);
	for(unsigned int i=0; i<len; ++i){
		g = data.getHeaderCell(0,i);
		gTable[i] = sqrt(1-g*g);
	}
	
	projections[0].reserve(len);
	projections[1].reserve(len);
	segments.reserve(len);
	segment.reserve(len);
	
	
	float c,h;
	c=0;
	for(unsigned int i=0; i<len; ++i){
		h=gTable[i];
		c+= doublepi*h;
	}
	C=1.0/c;
	
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




ScanlineMode MultiLayeredDepthBuffer::insertIntoLineBuffer(DepthBufferLine &line, float front, float back){
	DepthBufferLine::iterator it0, it1, it, it_next, it0bound, it1bound;
	pair<float,LineType> p0, p1;
	pair<DepthBufferLine::iterator, bool> r0,r1;
	bool frontChanges, backChanges, wouldClear;
	float x0;
	float frontBound,backBound;
	frontBound = front;
	backBound = back;
	

	frontChanges=false;
	backChanges=false;
	wouldClear=false;
	
	
	
	if(line.size()==0){
		p0.first=front;
		p0.second=FRONT;
		p1.first=back;
		p1.second=BACK;
		
		it0 = line.insert(p0);
		it1 = line.insert(p1);
 		return SCANLINE_PARTIAL;
	}
	
	it0bound = line.lower_bound(front);
	if(it0bound==line.end()) it0=line.begin();
	else it0=it0bound;
	
	it1bound = line.lower_bound(back);
	if(it1bound==line.end()) it1=line.begin();
	else it1=it1bound;
	

	//for sure it changes
	if(it0!=it1){
		if(it0->second==FRONT) frontChanges=true;
		if(it1->second==FRONT) backChanges=true;
	}
	//it still might change
	else{
		
		if(it0->second==FRONT){
			frontChanges=true;
			backChanges=true;
		}
		else{
			x0 = it0->first;
			if(front>=back) back+=doublepi;
			if(x0<front) x0+=doublepi;
			
			if(back>=x0){
				wouldClear=true;
			}
		}
	}
	
	
	if(wouldClear){
		line.clear();
 		return SCANLINE_FULL;
	}
	
	if(frontChanges){
		p0.first=frontBound;
		p0.second=FRONT;
		it0 = line.insert(it0bound,p0);
	}
	
	if(backChanges){
		p1.first=backBound;
		p1.second=BACK;
		//in certain occasions this is not optimal, because it0 might have changed it1's lower bound
		it1 = line.insert(it1bound,p1);
	}
	
	
	
	if(frontChanges) it=increaseLineInterator(it0,line);
	else it=it0;
	//printf("starting at %f-%d\n",it->first,it->second);
	while(it!=it1){
		it_next=increaseLineInterator(it,line);
		//printf("deleting %f-%d\n",it->first,it->second);
		line.erase(it);
		it=it_next;
	}
	
	
	if(line.size()==0) return SCANLINE_FULL;
	else return SCANLINE_PARTIAL;
	
		
}




bool MultiLayeredDepthBuffer::wouldChangeLineBuffer(DepthBufferLine &line, float front, float back, bool &frontChanges, bool &backChanges){
	DepthBufferLine::iterator it0, it1, it, it_next;
	bool wouldChange;
	float x0;
	
	
	
	frontChanges=false;
	backChanges=false;
	wouldChange=false;
	
	
	it0 = line.lower_bound(front);
	if(it0==line.end()) it0=line.begin();
	//if(it0->first!=front) it0=decreaseLineInterator(it0,line);
	
	it1 = line.lower_bound(back);
	if(it1==line.end()) it1=line.begin();
	

	//for sure it changes
	if(it0!=it1){
		wouldChange=true;
		if(it0->second==FRONT) frontChanges=true;
		if(it1->second==FRONT) backChanges=true;
	}
	//it still might change
	else{
		
		if(it0->second==FRONT){
			wouldChange=true;
			frontChanges=true;
			backChanges=true;
		}
		else{
			x0 = it0->first;
			if(front>=back) back+=doublepi;
			if(x0<front) x0+=doublepi;
			
			if(back>=x0){
				wouldChange=true;
			}
		}
	}
	

	return wouldChange;
		
}



Vector MultiLayeredDepthBuffer::normalise(Vector x){
	float l;
	Vector v(3);
	l = norm(x,2);
	v = x/l;
	return v;
	
}





int MultiLayeredDepthBuffer::sgn(float d){
	if(d>=0) return 1;
	else return -1;
}

void MultiLayeredDepthBuffer::addSphere(Vector &v, float lambda, bool invert, float &kappa, float &psi){
	Vector c;
	DepthInformation dat;
	Vector n,n2;
	float s;
	float dot_n_aux;
	
	
	
	psi=acos(dot(v,orientationAxis));
	
	
	
	
	//construct normal
	if(isWithinNumericalLimits(abs(dot(v,orientationAxis)),1.0)){
		kappa=0;
	}
	else{
		n2=cross(v, orientationAxis);
		n=normalise(cross(orientationAxis,n2));
		s = sgn(dot(n,orientationPlaneNormal));
		dot_n_aux=dot(n,orientationAxisAuxiliary);
		kappa = acos(dot_n_aux);
		if(s<0) kappa=doublepi-kappa;
	}
	
	
	
	
	
	
	//printf("ADDING SPHERE %d: lambda:%f psi:%f kappa:%f cart(%f %f %f) form: %d\n",index,lambda, psi, kappa, v(0), v(1), v(2), invert);
	//printf("n2: (%f %f %f)\n", n2(0), n2(1), n2(2));
	//printf("n: (%f %f %f)\n", n(0), n(1), n(2));

	dat = data.getFloorScanlines(kappa, psi, lambda, invert);
	
	for(unsigned int i=0; i<len; ++i){
		//printf("line[%d]...",i);
		
		if(dmode[i]!=SCANLINE_FULL){
			if(dat.mode[i]==SCANLINE_FULL){
				dmode[i]=SCANLINE_FULL;
				dbuffer[i].clear();
				//printf("F");
			}
			else if(dat.mode[i]==SCANLINE_PARTIAL){
				//printf("SCANLINE %f %f\n",dat.scanline0[i], dat.scanline1[i]);
				dmode[i]=insertIntoLineBuffer(dbuffer[i], dat.scanline0[i], dat.scanline1[i]);
				
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


bool MultiLayeredDepthBuffer::passesBuffer(float kappa, float psi, float lambda, bool invert, vector<Vector> &exposedVectors){
	DepthInformation dat;
	bool pass;
	dat = data.getCeilScanlines(kappa, psi, lambda, invert);
	pass=scanBuffer(kappa, exposedVectors, dat);
	if(invert) pass=true;
	return pass;
}


bool MultiLayeredDepthBuffer::getSplitterExposedVectors(vector<Vector> &exposedVectors){
	DepthInformation dat;
	bool res;
	dat.scanline0.resize(len,0);
	dat.scanline1.resize(len,M_PI);
	dat.mode.resize(len,SCANLINE_PARTIAL);
	res=scanBuffer(M_PI, exposedVectors, dat);
	return res;
}


bool MultiLayeredDepthBuffer::scanBuffer(float kappa, vector<Vector> &exposedVectors, DepthInformation &dat){
	Vector c;
	float offset=0;
	bool wouldChange;
	bool frontChanges, backChanges;
	DepthBufferLine::iterator it;
	DepthBufferCoordinate pr;
	float i_prev;
	
	
	projections[0].clear();
	projections[1].clear();
	segments.clear();
	segment.clear();
	
	//if(invert) return true;
	
	wouldChange=false;
	
	for(unsigned int i=0; i<len; ++i){
		if(dmode[i]==SCANLINE_PARTIAL){
			if(dat.mode[i]==SCANLINE_FULL){
				wouldChange=true;
				dmode[i]=SCANLINE_FULL;
				for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
					pr.i = i;
					pr.kappa = it->first;
					projections[it->second].push_back(pr);
				}
				
				
				
			}
			else if(dat.mode[i]==SCANLINE_PARTIAL){
				
				if(wouldChangeLineBuffer(dbuffer[i], dat.scanline0[i]-offset, dat.scanline1[i]+offset, frontChanges, backChanges)){
					wouldChange=true;
					if(frontChanges){
						pr.i = i;
						pr.kappa = dat.scanline0[i]-offset;
						projections[0].push_back(pr);
						
					}
					if(backChanges){
						pr.i = i;
						pr.kappa = dat.scanline1[i]+offset;
						projections[1].push_back(pr);
					}
					
				}
			}
			else{
			}
		}
		else if(dmode[i]==SCANLINE_EMPTY){
			if(dat.mode[i]!=SCANLINE_EMPTY){
				if(dat.mode[i]!=SCANLINE_FULL){
					pr.i = i;
					pr.kappa = dat.scanline0[i]-offset;
					projections[0].push_back(pr);
						
					pr.i = i;
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
			
		}
		
	}
	
	
	//divide the list into separate segment candidates
	for(unsigned int z=0; z<=1; ++z){
		segment.clear();
		i_prev=-1;
		for(unsigned int i=0; i<projections[z].size(); ++i){
			if(i!=0){
				if(abs(projections[z][i].i-i_prev)>1){
					segments.push_back(segment);
					segment.clear();
				}
			}
				
			i_prev=projections[z][i].i;
			segment.push_back(projections[z][i]);
			
		}
		if(segment.size()>0) segments.push_back(segment);
	}
	
	//now reduce projections to just one per segment
	for(unsigned int i=0; i<segments.size(); ++i){
		//exposedVectors.push_back(convertToCartesian(segments[i][segments[i].size()/2]));
		for(unsigned int j=0; j<segments[i].size(); ++j){
			exposedVectors.push_back(convertToCartesian(segments[i][j]));
		}
	}
	
	
	
	return wouldChange;
}



void MultiLayeredDepthBuffer::print(){
	DepthBufferLine::iterator it;
	vector<int> line;
	int pos;
	printf("----------------------------------------------------------------\n");
	
	for(unsigned int i=0; i<len; ++i){
		printf("[%d]:\t",i);
		for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
			printf("(%f %d) ",it->first,it->second);
		}
		printf("\n");
	}
	
	
	for(unsigned int i=0; i<len; ++i){
		printf("[%d]:\t",i);
		line.clear();
		line.resize(len+1,-1);
		for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
			pos=it->first*(float)len/(doublepi);
			line[pos]+=it->second+1;
		}
		
		for(unsigned int j=0; j<len+1; ++j){
			if(dmode[i]==SCANLINE_FULL) printf("#");
			else if(line[j]>=0) printf("%d",line[j]);
			else printf("*");
		}
		printf("\n");
	}
	printf("----------------------------------------------------------------\n");
		
}

void MultiLayeredDepthBuffer::getDepthBounds(int p, float d, DepthBufferLine::iterator &it0, DepthBufferLine::iterator &it1){
	it1 = dbuffer[p].upper_bound(d);
	if(it1==dbuffer[p].end()) it1 = dbuffer[p].begin();
	it0=decreaseLineInterator(it1,dbuffer[p]);
	
}





bool MultiLayeredDepthBuffer::isWithinNumericalLimits(float x, float l){
	if(abs(x-l)<=0.0001) return true;
	else return false;
}





void MultiLayeredDepthBuffer::startNewCycle(){
	exposed=1;
	occluded=1;
	prior_occluded=occludedDistribution.getAuxiliaryCell(0);
	prior_exposed=exposedDistribution.getAuxiliaryCell(0);
	//printf("priors: %f %f\n",prior_exposed, prior_occluded);
}

bool MultiLayeredDepthBuffer::isExposed(float x){
	unsigned int p;
	float l;
	
	exposedDistribution.closestGridPoint(x,p,l);
	
	//printf("adding %f [%d]\n",x,p);
	
	if(exposedDistribution.getDataCell(p)>occludedDistribution.getDataCell(p)) return true;
	else return false;
	
}


void MultiLayeredDepthBuffer::addProbe(float x){
	unsigned int p;
	float l;
	
	exposedDistribution.closestGridPoint(x,p,l);
	
	//printf("adding %f [%d]\n",x,p);
	
	exposed*=exposedDistribution.getDataCell(p);
	occluded*=occludedDistribution.getDataCell(p);
	
	//printf("prob: %f %f\n",exposed, occluded);
}

bool MultiLayeredDepthBuffer::isCycleExposed(){
	float e;
	e=exposedProbability();
	//printf("endprob: %f %f\n",e, 1-e);
	if(e > (1-e)) return true;
	else return false;
}


float MultiLayeredDepthBuffer::exposedProbability(){
	float pe, po;
	pe = exposed * prior_exposed;
	po = occluded * prior_occluded;
	if(pe+po == 0) return 0;
	
	return pe/(pe+po);
}




Vector MultiLayeredDepthBuffer::convertToCartesian(DepthBufferCoordinate &pr){
	Vector v(3);
	float h;
	int i;
	float g;
	float kappa;
	
	i=pr.i;
	g = data.getHeaderCell(0,i);
	kappa=pr.kappa;
	
	h=gTable[i];
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


float MultiLayeredDepthBuffer::exposedArea(){
	float area;
	DepthBufferLine::iterator it;
	bool started;
	float back;
	float lineRadiants;
	float h;
	float c;
	float p;
	
	area=0;
	back=0;
	
	for(unsigned int i=0; i<len; ++i){
		lineRadiants=0;
		if(dmode[i]==SCANLINE_EMPTY){
			lineRadiants=doublepi;
		}
		else if(dmode[i]==SCANLINE_PARTIAL){
			started=false;
			for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
				if(!started && it->second==FRONT) back = decreaseLineInterator(it,dbuffer[i])->first-doublepi;
				
				if(it->second==FRONT) lineRadiants+=it->first-back;
				if(it->second==BACK) back=it->first;
				started=true;
			}
		}
		p=lineRadiants/doublepi;
		h = gTable[i];
		c=doublepi*h*p;
		area += c*C;
	}
	
	return area;
}





