#include "multiLayeredDepthBuffer.h"




MultiLayeredDepthBuffer::MultiLayeredDepthBuffer(Depth3D &data, int m){
	this->data=data;
	this->len = data.parameter0Dim;
	//create buffer
	dbuffer.resize(len, DepthBufferLine());
	dmode.resize(len, SCANLINE_EMPTY);
	tessellationAxis=Vector(3);
	tessellationAxisAuxiliary=Vector(3);
	tessellationPlaneNormal=Vector(3);
	
	if(m==0){
		tessellationAxis(0)=1;
		tessellationAxis(1)=0;
		tessellationAxis(2)=0;
		
		tessellationAxisAuxiliary(0)=0;
		tessellationAxisAuxiliary(1)=1;
		tessellationAxisAuxiliary(2)=0;
		
		tessellationPlaneNormal(0)=0;
		tessellationPlaneNormal(1)=0;
		tessellationPlaneNormal(2)=1;
	}
	else{
		tessellationAxis(0)=0;
		tessellationAxis(1)=1;
		tessellationAxis(2)=0;
		
		tessellationAxisAuxiliary(0)=0;
		tessellationAxisAuxiliary(1)=0;
		tessellationAxisAuxiliary(2)=1;
		
		tessellationPlaneNormal(0)=1;
		tessellationPlaneNormal(1)=0;
		tessellationPlaneNormal(2)=0;
	}
	

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
	
	
	p0.first=front;
	p0.second=FRONT;
	p1.first=back;
	p1.second=BACK;
	
	it0 = line.insert(p0);
	it1 = line.insert(p1);
	
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




bool MultiLayeredDepthBuffer::wouldChangeLineBuffer(DepthBufferLine &line, double front, double back){
	DepthBufferLine::iterator it0, it1, it, it_next;
	bool wouldChange;
	pair<double,LineType> p0, p1;
	
	p0.first=front;
	p0.second=FRONT;
	p1.first=back;
	p1.second=BACK;
	
	it0 = line.insert(p0);
	it1 = line.insert(p1);
	
	wouldChange=false;
	
	it=decreaseLineInterator(it0,line);
	if(it->second==BACK) wouldChange=true;
	it=increaseLineInterator(it0,line);
	if(it->second==FRONT) wouldChange=true;
	
	it=increaseLineInterator(it1,line);
	if(it->second==FRONT) wouldChange=true;
	it=decreaseLineInterator(it1,line);
	if(it->second==BACK) wouldChange=true;
	
	it=increaseLineInterator(it0,line);
	while(it!=it1){
		it_next=increaseLineInterator(it,line);
		wouldChange=true;
		it=it_next;
	}
	

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





int MultiLayeredDepthBuffer::sgn(double d){
	if(d>=0) return 1;
	else return -1;
}

void MultiLayeredDepthBuffer::addSphere(CircularInterface &circle){
	Vector c;
	double x;
	double front,back;
	DepthInformation dat;
	double kappa;
	double psi;
	double lambda;
	Vector v;
	CircularInterfaceForm form;
	Vector n,n2;
	double s;
	
	lambda=circle.lambda.rotation;
	v=circle.normal;
	form=circle.form;
	
	
	psi=acos(dot(v,tessellationAxis));
	
	
	
	
	//construct normal
	n2=cross(v, tessellationAxis);
	n=normalise(cross(tessellationAxis,n2));
	s = sgn(dot(n,tessellationPlaneNormal));
	
	kappa = acos(dot(n,tessellationAxisAuxiliary));
	if(s<0) kappa=2*M_PI-kappa;
	
	//printf("ADDING SPHERE %d: lambda:%f psi:%f kappa:%f cart(%f %f %f) form: %d\n",circle.index,lambda, psi, kappa, v(0), v(1), v(2), form);
	//printf("n2: (%f %f %f)\n", n2(0), n2(1), n2(2));
	//printf("n: (%f %f %f)\n", n(0), n(1), n(2));

	dat = data.getFloorScanlines(kappa, psi, lambda, form);
	
	for(int i=0; i<len; ++i){
		//printf("line[%d]...",i);
		
		if(dmode[i]!=SCANLINE_FULL){
			if(dat.mode[i]==SCANLINE_FULL){
				dmode[i]=SCANLINE_FULL;
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



bool MultiLayeredDepthBuffer::passesBuffer(CircularInterface &circle){
	Vector c;
	double x;
	double front,back;
	DepthInformation dat;
	double kappa;
	double psi;
	double lambda;
	Vector v;
	CircularInterfaceForm form;
	Vector n,n2;
	double s;
	double offset=0.01;
	
	int ic=38768767625;
	
	if(circle.form==CONCAVE) return true;
	
	lambda=circle.lambda.rotation;
	v=circle.normal;
	form=circle.form;
	
	
	psi=acos(dot(v,tessellationAxis));
	
	
	//construct normal
	n2=cross(v, tessellationAxis);
	n=normalise(cross(tessellationAxis,n2));
	s = sgn(dot(n,tessellationPlaneNormal));
	
	kappa = acos(dot(n,tessellationAxisAuxiliary));
	if(s<0) kappa=2*M_PI-kappa;
	

	dat = data.getCeilScanlines(kappa, psi, lambda, form);
	
	for(int i=0; i<len; ++i){
		
		if(circle.index==ic){
			printf("line[%d]: ",i);
		}
		if(dmode[i]!=SCANLINE_FULL){
			if(dat.mode[i]==SCANLINE_FULL){
				return true;
				dmode[i]=SCANLINE_FULL;
				if(circle.index==ic){
					printf("F:");
				}
				
			}
			else if(dat.mode[i]==SCANLINE_PARTIAL){
				if(circle.index==ic){
					printf("P: %f %f",dat.scanline0[i], dat.scanline1[i]);
				}
				
				if(wouldChangeLineBuffer(dbuffer[i], dat.scanline0[i]-offset, dat.scanline1[i]+offset)) return true;
			}
			else{
				if(circle.index==ic){
					printf("E: ");
				}
			}
		}
		else{
			if(dat.mode[i]==SCANLINE_FULL){
			}
			
			if(circle.index==ic){
				if(dat.mode[i]==SCANLINE_FULL){
					printf("+: ");
				}
				else{
					printf("#:");
				}
			}
		}
		if(circle.index==ic){
			printf("\n");
		}
		
	}
	
	
	return false;
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



























