#include "triforceInterface.h"

#include <algorithm>
#include <string>

using namespace std;


TriforceInterface::TriforceInterface(string path, unsigned int buffer, unsigned int slack){
	benchmark=Benchmark(string("Interface"));
	benchmark.start(string("loading tables"));
	df0 = new DataFile(path+"/dataConcave.dat");
	df1 = new DataFile(path+"/dataConcave0.dat");
	df2 = new DataFile(path+"/dataConcave1.dat");
	df3 = new DataFile(path+"/dataConcave2.dat");
	
	dat0 = df0->digest3DBinaryTable();
	dat1 = df1->digest3DBinaryTable();
	dat2 = df2->digest3DBinaryTable();
	dat3 = df3->digest3DBinaryTable();
	
	surf0 = new Surface3D(dat0);
	surf1 = new Surface3D(dat1);
	surf2 = new Surface3D(dat2);
	surf3 = new Surface3D(dat3);
	
	interpolator0 = new Interpolation(surf0,TAYLOR_QUADRATIC);
	interpolator1 = new Interpolation(surf1,TAYLOR_QUADRATIC);
	interpolator2 = new Interpolation(surf2,TAYLOR_QUADRATIC);
	interpolator3 = new Interpolation(surf3,TAYLOR_QUADRATIC);
	
	interpolators.clear();
	interpolators.push_back(interpolator0);
	interpolators.push_back(interpolator1);
	interpolators.push_back(interpolator2);
	interpolators.push_back(interpolator3);
	

	integrator = new IntegratorTriforce(interpolators);
	withDepthBuffer=false;
	
	if(buffer>0){
		withDepthBuffer=true;
		this->buffer=buffer;
		df4 = new DataFile(path+"/depthBuffer.dat");
		df5 = new DataFile(path+"/occludedDistribution.dat");
		df6 = new DataFile(path+"/exposedDistribution.dat");

		dat4 = df4->digest3DBinaryTable();
		dat5 = df5->digest1DBinaryTable();
		dat6 = df6->digest1DBinaryTable();

		depth0 = new Depth3D(dat4,slack);
	}
	benchmark.stop();
	
}




TriforceInterface::TriforceInterface(string path){
	benchmark=Benchmark(string("Interface"));
	benchmark.start(string("loading tables"));
	
	df0 = new DataFile(path+"/dispersionFieldEps.dat");
	df1 = new DataFile(path+"/dataConcave0.dat");
	df2 = new DataFile(path+"/dataConcave1.dat");
	df3 = new DataFile(path+"/dataConcave2.dat");
	
	dat0 = df0->digest3DBinaryTable();
	dat1 = df1->digest3DBinaryTable();
	dat2 = df2->digest3DBinaryTable();
	dat3 = df3->digest3DBinaryTable();
	
	surf0 = new Surface3D(dat0);
	surf1 = new Surface3D(dat1);
	surf2 = new Surface3D(dat2);
	surf3 = new Surface3D(dat3);
	
	interpolator0 = new Interpolation(surf0,TAYLOR_QUADRATIC);
	interpolator1 = new Interpolation(surf1,TAYLOR_QUADRATIC);
	interpolator2 = new Interpolation(surf2,TAYLOR_QUADRATIC);
	interpolator3 = new Interpolation(surf3,TAYLOR_QUADRATIC);
	
	interpolators.clear();
	interpolators.push_back(interpolator0);
	interpolators.push_back(interpolator1);
	interpolators.push_back(interpolator2);
	interpolators.push_back(interpolator3);
	

	integrator = new IntegratorTriforce(interpolators);
	withDepthBuffer=false;
	
	if(buffer>0){
		withDepthBuffer=true;
		this->buffer=buffer;
		df4 = new DataFile(path+"/depthBuffer.dat");
		df5 = new DataFile(path+"/occludedDistribution.dat");
		df6 = new DataFile(path+"/exposedDistribution.dat");

		dat4 = df4->digest3DBinaryTable();
		dat5 = df5->digest1DBinaryTable();
		dat6 = df6->digest1DBinaryTable();
	}
	benchmark.stop();
	
}


float TriforceInterface::calculateSurfaceArea(Molecule &mol){
	Tessellation *t;
	float area;
	if(withDepthBuffer) t = new Tessellation(mol,buffer,*depth0,*dat5,*dat6);
	else t = new Tessellation(mol);
	t->build(true);
	t->outputTessellation(string("patches.csv"));
	area = integrator->integrate(&mol, t);
	tessellationBenchmark=t->getBenchmark();
	delete(t);
	
	return area;
}



Benchmark TriforceInterface::getBenchmark(){
	return benchmark;
}


void TriforceInterface::printBenchmark(FILE* outputfile){
	fprintf(outputfile,"time statistics\n");
	benchmark.print(outputfile);
	tessellationBenchmark.print(outputfile);
	integrator->getBenchmark().print(outputfile);
	
}


