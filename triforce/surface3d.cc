#include "surface3d.h"

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;

#define THRESHOLD_NUMERICAL 0.00001


Surface3D::Surface3D(Data3D* d)
{
	parameter0Dim=d->parameter0Dim;
	parameter1Dim=d->parameter1Dim;
	parameter2Dim=d->parameter2Dim;
	minParameter0=d->minParameter0;
	minParameter1=d->minParameter1;
	minParameter2=d->minParameter2;
	cellLengthParameter0=d->cellLengthParameter0;
	cellLengthParameter1=d->cellLengthParameter1;
	cellLengthParameter2=d->cellLengthParameter2;
	headerParameter0=d->headerParameter0;
	headerParameter1=d->headerParameter1;
	headerParameter2=d->headerParameter2;
	data=d->data;
	gradient=d->gradient;
	hessian=d->hessian;	
	auxiliary=d->auxiliary;

}




void Surface3D::surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths){
	VectorInt v2(3);
	VectorInt v(3);
	unsigned int discontinuity;
	double vpsi,vlambda;
	closestGridPoint(x, v, lengths);
	
	//printf("CLOSEST GRIDPOINT: %d %d %d\n",v(0),v(1),v(2));
	
	if(v(0)>=parameter0Dim || v(1)>=parameter1Dim || v(2)>=parameter2Dim) return;
	
	//decide on which side of the discontinuities x lies.
	
	if(x(1) < x(2)){
		discontinuity = 0;
	}
	else if(x(1)+x(2) < M_PI){
		discontinuity = 1;
	}
	else{
		discontinuity = 2;
	}
	
	
	
	
	//printf("DISCONTINUITY: %d\n",discontinuity);
	
	
	r.clear();
	for(unsigned int i=0;i<2;++i)
		for(unsigned int j=0;j<2;++j)
			for(unsigned int k=0;k<2;++k){
				v2(0)=v(0)+i;
				v2(1)=v(1)+j;
				v2(2)=v(2)+k;
				//printf("v2: %d %d %d, par: %d %d %d, bool: %d\n",v2(0),v2(1),v2(2),parameter0Dim,parameter1Dim,parameter2Dim, v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim);
				if(v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim){
					if(!isnan((*data)[v2(0)][v2(1)][v2(2)])){
						vpsi = (*headerParameter1)[v2(1)];
						vlambda = (*headerParameter2)[v2(2)];
						
						if(discontinuity==0){
							if(vpsi < vlambda) r.push_back(v2);
						}
						else if(discontinuity==1){
							if(vpsi >= vlambda && vpsi + vlambda < M_PI && !isWithinNumericalLimits(vpsi + vlambda,M_PI)) r.push_back(v2);
						}
						else{
							if(vpsi >= vlambda && vpsi + vlambda >= M_PI) r.push_back(v2);
						}
							//printf("ACCEPTED %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));
						
						//else	printf("REJECTED 0 %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));

					}
					//else	printf("REJECTED 1 %d %d %d\n",v2(0),v2(1),v2(2));

				}
				//else	printf("PRE-REJECTED 2 %d %d %d\n",v2(0),v2(1),v2(2));
			}
			
			
			
	if(r.size()==0){
		for(unsigned int i=0;i<2;++i)
			for(unsigned int j=0;j<2;++j)
				for(unsigned int k=0;k<2;++k){
					v2(0)=v(0)+i;
					v2(1)=v(1)+j;
					v2(2)=v(2)+k;
					if(v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim){
						if(!isnan((*data)[v2(0)][v2(1)][v2(2)])){
							//printf("RE-ACCEPTED %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));
								r.push_back(v2);
						}

					}
				}
		
	}

			
}
