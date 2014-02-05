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
	
	interpolator0 = new Interpolation(dat0,TAYLOR_LINEAR);
	interpolator1 = new Interpolation(dat1,TAYLOR_LINEAR);
	interpolator2 = new Interpolation(dat2,TAYLOR_LINEAR);
	interpolator3 = new Interpolation(dat3,TAYLOR_LINEAR);
	
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
	
	interpolator0 = new Interpolation(dat0,TAYLOR_QUADRATIC);
	interpolator1 = new Interpolation(dat1,TAYLOR_QUADRATIC);
	interpolator2 = new Interpolation(dat2,TAYLOR_QUADRATIC);
	interpolator3 = new Interpolation(dat3,TAYLOR_QUADRATIC);
	
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
	float area;
	if(withDepthBuffer) t = new Tessellation(mol,buffer,*depth0,*dat5,*dat6);
	else t = new Tessellation(mol);
	t->build(true);
// 	t->outputTessellation(string("patches.csv"));
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




void TriforceInterface::minimise(Molecule &mol0, Molecule &mol1){
	float area0,area1;
	float step=0.001;
	Tessellation *tref;
	vector<vector<ForcesDT*> >forces;
	Vector f(3);
	vector<AreasDT*> areas0,areas1;
	float diff,total;
	if(withDepthBuffer){
		t = new Tessellation(mol0,buffer,*depth0,*dat5,*dat6);
	}
	else{
		t = new Tessellation(mol0);
	}

	
	for(unsigned int s=0; s<10; ++s){
		mol0.update();
		t->update();
		printf("BUILDING 0\n");
	
		t->build(true);
		area0 = integrator->integrate(&mol0, t);
		forces=mol0.fetchForcePointers();
		
		mol1.update();
		tref = new Tessellation(mol1);
		printf("BUILDING 1\n");
		tref->build(true);
		area1 = integrator->integrate(&mol1, tref);
		
		total=0;
		areas0=mol0.fetchAreaPointers();
		areas1=mol1.fetchAreaPointers();
		//compare
		for(unsigned int i=0; i<mol1.fetchCoordinates().size(); ++i){
			diff=*(areas0[i])- *(areas1[i]);
			total=total+diff*diff;
		}
		total=sqrt(total);
		delete(tref);
		
		for(unsigned int i=0; i<mol0.fetchCoordinates().size(); ++i){
			f(0)=*(forces[i][0]);
			f(1)=*(forces[i][1]);
			f(2)=*(forces[i][2]);
			f*=-step;
			mol0.perturbInternallyStoredAtomCoordinates(i,f);
			Vector vec(3);
			vec=mol0.getInternallyStoredAtomCoordinates(i);
			mol1.setInternallyStoredAtomCoordinates(i,vec);
		}
		
		
		
		printf("TOTAL DIFFERENCE: %f AREAS: %f %f\n",total,area0, area1);
	}
	tessellationBenchmark=t->getBenchmark();
	
}

