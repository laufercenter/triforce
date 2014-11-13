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
	df7 = new DataFile(path+"/barycentricWeights.dat");
	
	dat0 = df0->digest3DBinaryTable();
	dat1 = df1->digest3DBinaryTable();
	dat2 = df2->digest3DBinaryTable();
	dat3 = df3->digest3DBinaryTable();
	dat7 = df7->digest3DBinaryVectorialTable();
	
	
	interpolator0 = new Interpolation(dat0,TAYLOR_LINEAR);
	interpolator1 = new Interpolation(dat1,TAYLOR_LINEAR);
	interpolator2 = new Interpolation(dat2,TAYLOR_LINEAR);
	interpolator3 = new Interpolation(dat3,TAYLOR_LINEAR);
/*	
	interpolator0 = new InterpolationPolytopical(dat0,dat7);
	interpolator1 = new InterpolationPolytopical(dat1,dat7,interpolator0);
	interpolator2 = new InterpolationPolytopical(dat2,dat7,interpolator0);
	interpolator3 = new InterpolationPolytopical(dat3,dat7,interpolator0);
*/	
	
	
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
	df7 = new DataFile(path+"/barycentricWeights.dat");
	
	dat0 = df0->digest3DBinaryTable();
	dat1 = df1->digest3DBinaryTable();
	dat2 = df2->digest3DBinaryTable();
	dat3 = df3->digest3DBinaryTable();
	dat7 = df7->digest3DBinaryVectorialTable();
	
	
	interpolator0 = new Interpolation(dat0,TAYLOR_LINEAR);
	interpolator1 = new Interpolation(dat1,TAYLOR_LINEAR);
	interpolator2 = new Interpolation(dat2,TAYLOR_LINEAR);
	interpolator3 = new Interpolation(dat3,TAYLOR_LINEAR);
/*	
	interpolator0 = new InterpolationPolytopical(dat0,dat7);
	interpolator1 = new InterpolationPolytopical(dat1,dat7);
	interpolator2 = new InterpolationPolytopical(dat2,dat7);
	interpolator3 = new InterpolationPolytopical(dat3,dat7);
*/	
	
	interpolators.clear();
	interpolators.push_back(interpolator0);
	interpolators.push_back(interpolator1);
	interpolators.push_back(interpolator2);
	interpolators.push_back(interpolator3);
	

	integrator = new IntegratorTriforce(interpolators);
	withDepthBuffer=false;
	
	if(buffer>0){
		withDepthBuffer=true;
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
	//delete(t); //THIS NEEDS TO BE UN-DECOMMENTED
	
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


void TriforceInterface::printTessellation(FILE* outputfile){
	t->print(outputfile);
}



void TriforceInterface::printSurfaces(Molecule &mol, FILE* outputfile){
// 	t = new Tessellation(mol);
// 	t->build(false, false);
	t->print(outputfile);
	delete(t);
	
}






void TriforceInterface::minimise(Molecule &mol0){
	float area0,area1;
	float step=0.0001;
	Tessellation *tref;
	vector<vector<ForcesDT*> >forces;
	Vector f(3);
	vector<AreasDT*> areas0,areas1;
	float diff,total;
	vector<float> mag;
	vector<float> mag2;
	float med;
	float mad;
	float AR;
	if(withDepthBuffer){
		t = new Tessellation(mol0,buffer,*depth0,*dat5,*dat6);
	}
	else{
		t = new Tessellation(mol0);
	}

	printf("step\tarea\tarea2\n");
	
	for(unsigned int s=0; s<10000; ++s){


		
		mol0.update();
		t->update();
		t->build(true);
		area0 = integrator->integrate(&mol0, t);
		forces=mol0.fetchForcePointers();
		

		total=0;
		areas0=mol0.fetchAreaPointers();
		AR=0;
		mag.clear();
		for(unsigned int i=0; i<mol0.fetchCoordinates().size(); ++i){
			f(0)=*(forces[i][0]);
			f(1)=*(forces[i][1]);
			f(2)=*(forces[i][2]);
			mag.push_back(norm(f,2));
			AR+=*(areas0[i]);
		}
		std::sort(mag.begin(),mag.end());
		med=mag[(unsigned int)floor(mag.size()/2.0)];
		mag2.clear();
		for(unsigned int i=0; i<mol0.fetchCoordinates().size(); ++i){
			mag2.push_back(fabs(mag[i]-med));
		}
		std::sort(mag2.begin(),mag2.end());
		mad=mag2[(unsigned int)floor(mag2.size()/2.0)];
		
		
		
		
		for(unsigned int i=0; i<mol0.fetchCoordinates().size(); ++i){
			f(0)=*(forces[i][0]);
			f(1)=*(forces[i][1]);
			f(2)=*(forces[i][2]);
			
			//if((mag[i]-med)/mad >= 0.2) f=Vector(3).zeros();
			
			f*=-step;
			mol0.perturbInternallyStoredAtomCoordinates(i,f);
			Vector vec(3);
			vec=mol0.getInternallyStoredAtomCoordinates(i);
		
			
		}
		printf("%d\t%f\t%f\n",s,area0,AR);

		
			
	}
	tessellationBenchmark=t->getBenchmark();
	

	
}

