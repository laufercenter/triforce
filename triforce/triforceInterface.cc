#include "triforceInterface.h"

#include <algorithm>
#include <string>

using namespace std;


TriforceInterface::TriforceInterface(string path){
	
	df0 = new DataFile(path+"/dataConcave.dat");
	df1 = new DataFile(path+"/dataConcave0.dat");
	df2 = new DataFile(path+"/dataConcave1.dat");
	df3 = new DataFile(path+"/dataConcave2.dat");
	df4 = new DataFile(path+"/depthBuffer.dat");
	
	dat0 = df0->digest3DBinaryTable();
	dat1 = df1->digest3DBinaryTable();
	dat2 = df2->digest3DBinaryTable();
	dat3 = df3->digest3DBinaryTable();
	dat4 = df4->digest3DBinaryTable();
	
	surf0 = new Surface3D(dat0);
	surf1 = new Surface3D(dat1);
	surf2 = new Surface3D(dat2);
	surf3 = new Surface3D(dat3);
	depth0 = new Depth3D(dat4);
	
	interpolator0 = new Interpolation(surf0,TAYLOR_QUADRATIC);
	interpolator1 = new Interpolation(surf1,TAYLOR_QUADRATIC);
	interpolator2 = new Interpolation(surf2,TAYLOR_QUADRATIC);
	interpolator3 = new Interpolation(surf3,TAYLOR_QUADRATIC);
	

	integrator = new IntegratorTriforce(interpolator0, interpolator1, interpolator2, interpolator3);
}



double TriforceInterface::calculateSurfaceArea(Molecule &mol){
	Tessellation *t;
	double area;
	t = new Tessellation(mol,*depth0);
	t->build(true);
	area = integrator->integrate(&mol, t);
	delete(t);
	
	return area;
}
