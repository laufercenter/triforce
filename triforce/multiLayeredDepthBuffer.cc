#include "multiLayeredDepthBuffer.h"




MultiLayeredDepthBuffer::MultiLayeredDepthBuffer(int detail){
	this->len = detail;
	//create buffer
	dbuffer.resize(len, DepthBufferLine());

}



Vector MultiLayeredDepthBuffer::sphericalVector(double phi, double theta){
	Vector v(3);
	v(0) = 1;
	v(1) = theta;
	v(2) = phi;
	v = spherical2cartesian(v);
	
	return v;
}

Vector MultiLayeredDepthBuffer::spherical2cartesian(Vector s) {
	Vector v(3);
	v(0) = s(0) * cos(s(1)) * sin(s(2));
	v(1) = s(0) * sin(s(1)) * sin(s(2));
	v(2) = s(0) * cos(s(2));
	return v;
}


Vector MultiLayeredDepthBuffer::cartesian2spherical(Vector v) {
	Vector s(3);
	s(0) = norm(v,2);
	s(1) = atan(v(1)/v(0));
	s(2) = acos(v(3)/s(0));
	return s;
}

bool MultiLayeredDepthBuffer::circlef(double pos, Vector &c, double lambda, double &front, double &back){
	double y;
	double y_offset;
	double theta;
	double phi;
	double x;
	
	theta = c(1);
	phi = c(2);
	y_offset = phi;
	x = pos-theta;
	
	if(abs(x)>lambda) return false;
	
	y = sqrt(lambda*lambda - x*x);
	
	front = y_offset - y;
	back = y_offset + y;
	
	return true;
	
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
	
	if(it!=line.begin()){
		--it;
		if(it->second==BACK) wouldChange=true;
	}
	else wouldChange=true;

	it = it1;
	++it;
	if(it!=line.end()){
		if(it->second==FRONT) wouldChange=true;
	}
	else wouldChange=true;

	line.erase(it0);
	line.erase(it1);
	
	return wouldChange;
		
}


void MultiLayeredDepthBuffer::addSphere(Vector &v, double lambda){
	Vector c;
	double x;
	double front,back;
	
	lambda-=0.01;
	
	c = cartesian2spherical(v);
	
	for(int i=0; i<len; ++i){
		x = i*((2*M_PI) / len);
		if(circlef(x, c, lambda, front, back)){
			if(front < 0){
				insertIntoLineBuffer(dbuffer[i],0,back);
				insertIntoLineBuffer(dbuffer[i],2*M_PI+front,2*M_PI);
			}
			else if(back > 2*M_PI){
				insertIntoLineBuffer(dbuffer[i],front,2*M_PI);
				insertIntoLineBuffer(dbuffer[i],0,2*M_PI-back);
			}
			else{
				insertIntoLineBuffer(dbuffer[i],front,back);
			}
		}
	}
	
	
}



bool MultiLayeredDepthBuffer::passesBuffer(Vector &v, double lambda){
	Vector c;
	double x;
	double front,back;
	
	c = cartesian2spherical(v);
	
	
	for(int i=0; i<len; ++i){
		x = i*((2*M_PI) / len);
		if(circlef(x, c, lambda, front, back)){
			if(front < 0){
				if(wouldChangeLineBuffer(dbuffer[i],0,back)) return false;
				if(wouldChangeLineBuffer(dbuffer[i],2*M_PI+front,2*M_PI)) return false;
			}
			else if(back > 2*M_PI){
				if(wouldChangeLineBuffer(dbuffer[i],front,2*M_PI)) return false;
				if(wouldChangeLineBuffer(dbuffer[i],0,2*M_PI-back)) return false;
			}
			else{
				if(wouldChangeLineBuffer(dbuffer[i],front,back)) return false;
			}
		}
		
	}
	
	return true;
}































