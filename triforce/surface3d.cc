#include "surface3d.h"

#include <algorithm>
#include <string>
#include <limits>
#include <math.h>


using namespace std;
using namespace arma;
using namespace boost;

#define THRESHOLD_NUMERICAL 0.00001


Surface3D::Surface3D(Data3D* d): Data3D()
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

}




void Surface3D::surroundingPointsAndCellLengths(Vector &x, vector<VectorInt> &r, Vector &lengths){
	VectorInt v2(3);
	int i,j,k;
	VectorInt v(3);
	bool neg;
	
	closestGridPoint(x, v, lengths);
	
	if((*headerParameter1)[v(1)]+(*headerParameter2)[v(2)] < M_PI && !isWithinNumericalLimits((*headerParameter1)[v(1)]+(*headerParameter2)[v(2)],M_PI)) neg=true;
	else neg=false;
	
	
	//printf("CLOSEST GRIDPOINT: %d %d %d\n",v(0),v(1),v(2));
	
	
	r.clear();
	for(i=0;i<2;++i)
		for(j=0;j<2;++j)
			for(k=0;k<2;++k){
				v2(0)=v(0)+i;
				v2(1)=v(1)+j;
				v2(2)=v(2)+k;
				if(v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim){
					if(!isnan((*data)[v2(0)][v2(1)][v2(2)])){
						if(!neg || ((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)]<M_PI && !isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI))){
							r.push_back(v2);
							printf("ACCEPTED %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));
						}
						else	printf("REJECTED 0 %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));

					}
					else	printf("REJECTED 1 %d %d %d\n",v2(0),v2(1),v2(2));

				}
				else	printf("PRE-REJECTED 2 %d %d %d\n",v2(0),v2(1),v2(2));
			}
			
	if(r.size()==0){
		for(i=0;i<2;++i)
			for(j=0;j<2;++j)
				for(k=0;k<2;++k){
					v2(0)=v(0)+i;
					v2(1)=v(1)+j;
					v2(2)=v(2)+k;
					if(v2(0)<parameter0Dim && v2(1)<parameter1Dim && v2(2)<parameter2Dim){
						if(!isnan((*data)[v2(0)][v2(1)][v2(2)])){
							printf("RE-ACCEPTED %d %d %d (%f %f [%d/%d])\n",v2(0),v2(1),v2(2),(*headerParameter1)[v2(1)],(*headerParameter2)[v2(2)],neg,isWithinNumericalLimits((*headerParameter1)[v2(1)]+(*headerParameter2)[v2(2)],M_PI));
								r.push_back(v2);
						}

					}
				}
		
	}

			
}
