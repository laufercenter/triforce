#include "multiLayeredDepthBuffer.h"




MultiLayeredDepthBuffer::MultiLayeredDepthBuffer(Depth3D data){
	this->data=data;
	this->len = data.parameter0Dim;
	//create buffer
	dbuffer.resize(len, DepthBufferLine());
	dmode.resize(len, EMPTY);

}






void MultiLayeredDepthBuffer::insertIntoLineBuffer(DepthBufferLine &line, double front, double back){
	DepthBufferLine::iterator it0, it1, it, it_next;
	bool deleteBack, deleteFront;
	pair<double,LineType> p0, p1;
	
	p0.first=front;
	p0.second=FRONT;
	p1.first=back;
	p1.second=BACK;
	
	it0 = line.insert(p0);
	it1 = line.insert(p1);
	//prune linebuffer
	it = it0;
	if(it!=line.begin()){
		--it;
		if(it->second==FRONT) deleteFront=true;
	}
	
	it = it1;
	++it;
	if(it!=line.end()){
		if(it->second==BACK) deleteBack=true;
	}
	
	
	it = it0;
	++it;
	while(it!=it1){
		it_next=it;
		++it_next;
		line.erase(it);
		it=it_next;
	}
	
	if(deleteFront) line.erase(it0);
	if(deleteBack) line.erase(it1);
		
}




bool MultiLayeredDepthBuffer::wouldChangeLineBuffer(DepthBufferLine &line, double front, double back){
	DepthBufferLine::iterator it0, it1, it;
	bool wouldChange;
	pair<double,LineType> p0, p1;
	
	p0.first=front;
	p0.second=FRONT;
	p1.first=back;
	p1.second=BACK;
	
	it0 = line.insert(p0);
	it1 = line.insert(p1);
	
	wouldChange=false;
	
	it = it0;
	if(it!=line.begin() && line.size()!=0){
		--it;
		if(it->second==BACK) wouldChange=true;
	}
	else wouldChange=true;

	it = it1;
	++it;
	if(it!=line.end() && line.size()!=0){
		if(it->second==FRONT) wouldChange=true;
	}
	else wouldChange=true;

	line.erase(it0);
	line.erase(it1);
	
	return wouldChange;
		
}



Vector MultiLayeredDepthBuffer::normalise(Vector x){
	double l;
	Vector v(3);
	l = norm(x,2);
	v = x/l;
	return v;
	
}


void MultiLayeredDepthBuffer::addSphere(Vector tessellationPlaneNormal, Vector tessellationAxisAuxiliary, Vector tessellationAxis, Vector &v, double psi, double lambda){
	Vector c;
	double x;
	double front,back;
	DepthInformation dat;
	Vector n;
	double kappa;
	
	
	printf("ADDING SPHERE: %f sph(%f %f) cart(%f %f %f)\n",lambda, c(1), c(2), v(0), v(1), v(2));
	
	
	//construct normal
	n2=cross(v, tessellationAxis);
	n=normalise(cross(tessellationAxis,n2));
	s = sign(dot(n,tessellationPlaneNormal));
	
	kappa = acos(dot(n,tessellationAxisAuxiliary));
	if(s>0) kappa+=2*M_PI;
	

	dat = data.getScanlines(kappa, psi, lambda);
	
	for(int i=0; i<len; ++i){
		if(dmode[i]!=OCCUPIED){
			if(dat.mode[i]==OCCLUDED){
				dmode[i]=OCCUPIED;
			}
			else if(dat.mode[i]==PARTIAL){
				insertIntoLineBuffer(dbuffer[i], dat.scanline0, dat.scanline1);
			}
		}
	}
	
	print();
	
	
}



bool MultiLayeredDepthBuffer::passesBuffer(Vector &v, double lambda){
	Vector c;
	double x;
	double front,back;
	
	c = cartesian2spherical(v);
	
	
	for(int i=0; i<len; ++i){
		
	}
	
	return false;
}



void MultiLayeredDepthBuffer::print(){
	DepthBufferLine::iterator it;
	vector<int> line;
	int pos;
	printf("----------------------------------------------------------------\n");
	
	for(int i=0; i<len; ++i){
		for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
			printf("(%f %d) ",it->first,it->second);
		}
		printf("\n");
	}
	
	
	for(int i=0; i<len; ++i){
		line.clear();
		line.resize(len+1,-1);
		for(it=dbuffer[i].begin(); it!=dbuffer[i].end(); ++it){
			pos=it->first*(double)len/(2.0*M_PI);
			line[pos]=it->second;
		}
		
		for(int j=0; j<len+1; ++j){
			if(line[j]>=0) printf("%d",line[j]);
			else printf("*");
		}
		printf("\n");
	}
	printf("----------------------------------------------------------------\n");
		
}



























