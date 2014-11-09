#include "tessellation.h"
#include <limits>
#include <boost/math/distributions/normal.hpp>


Tessellation::Tessellation(Molecule &m){
	molecule = m;
	hasDepthBuffer=false;
	benchmark=Benchmark(string("Tessellation"));

}


Tessellation::Tessellation(Molecule &m, unsigned int numbBuffer, Depth3D &depthData, Data1D &occludedDistribution, Data1D &exposedDistribution){
	molecule = m;
	this->depthData = depthData;
	this->occludedDistribution = occludedDistribution;
	this->exposedDistribution = exposedDistribution;
	hasDepthBuffer=true;
	this->numbBuffer=numbBuffer;
	
	benchmark=Benchmark(string("Tessellation"));


}



void Tessellation::update(){
	if(hasDepthBuffer)
		prevSasasForMolecule=SASAs(sasasForMolecule);
	molecule.update();
}	
/*	
	bool split;
	CircularInterfacesPerAtom circlesPerAtom;
	vector<int> neighbourlist;
	
	//always perform a split, no-split is not supported to this date
	split=true;
	
	
	molecule.update();
	atoms = molecule.fetchCoordinates();
	radii = molecule.fetchRadii();
	
	if(hasDepthBuffer)
		prevSasasForMolecule=SASAs(sasasForMolecule);
	
	
	
	sasasForMolecule.clear();
	//iterate over all atoms
	for(unsigned int i=0; i<atoms.size(); ++i){
		benchmark.start(string("initialisation"));
		
		neighbourlist = molecule.getNeighborListFor(i);

		benchmark.stop();
		benchmark.start(string("gauss-bonnet path"));
	
		
		buildGaussBonnetPath(i, atoms, radii, sasasForMolecule, split, neighbourlist, useDepthBuffer);
		
		benchmark.stop();
	}
	
	
	
}

*/




void Tessellation::build(bool useDepthBuffer, bool split){
	CircularInterfacesPerAtom circlesPerAtom;
	vector<int> neighbourlist;
	vector<int> closestNeighbours;
	
	//always perform a split, no-split is not supported to this date
	
	//molecule.update();
	atoms = molecule.fetchCoordinates();
	radii = molecule.fetchRadii();
	closestNeighbours = molecule.fetchClosestNeighbours();
	totalInterfaces=0;
	survivedInterfaces=0;
	
	//for(int j=0; j<=10; ++j){
	sasasForMolecule.clear();
	//iterate over all atoms and build the tessellation for each of them
	for(unsigned int i=0; i<atoms.size(); ++i){
		//printf("COORDS %f %f %f\n",atoms[i](0),atoms[i](1),atoms[i](2));
		benchmark.start(string("initialisation"));
		
		//int i=12;{Tessellation
		neighbourlist = molecule.getNeighborListFor(i);

		
		//printf("NEIGHBOURS RECEIVED: %d\n",neighbourlist.size());
		
		
		benchmark.stop();
		benchmark.start(string("gauss-bonnet path"));
	
		buildGaussBonnetPath(i, atoms, radii, sasasForMolecule, split, neighbourlist, closestNeighbours[i], useDepthBuffer);
		
		benchmark.stop();
		
		
		/*
		if(i==10){
		fprintf(stderr,"--GBATOM-- %d\n",i);
			for(int j=0; j<sasasForMolecule[i].sasas.size(); ++j){
				fprintf(stderr,"SASA %d (%d)\n",j,sasasForMolecule[i].sasas[j].hemisphere);
				outputGaussBonnetPath(sasasForMolecule[i].sasas[j]);
			}
		}
		*/
		//exit(0);
	}
	
	//printf("Interfaces: %d kept: %d\n",totalInterfaces,survivedInterfaces);
	
	//}
	
	
}






CircularInterfacesPerAtom Tessellation::coverHemisphere(TessellationAxis &tessellationAxis, CircularInterfacesPerAtom circles, CircularInterfaceForm form){
	CircularInterface C;
	Vector v(3);
	Matrix NullMatrix(3,3);
	NullMatrix.zeros();
	
	v = tessellationAxis.v;
	
	C.lambdaRotation.rotation = M_PI * 0.5;
	C.lambda.rotation = M_PI * 0.5;
	C.psi.rotation=0;
	C.g=0;
	C.normal = v;
	C.form = form;
	C.id = circles.size();
	C.valid = true;
	C.index = -1;
	C.hasDerivatives=true;
	C.extended=false;
	
	C.dmu_dx = NullMatrix;
	C.lambda.drotation_dxi = Vector(3).zeros();
	C.lambda.drotation_dxj = Vector(3).zeros();
	C.lambda.drotation_dxl = Vector(3).zeros();
	C.lambda.drotation_dxt = Vector(3).zeros();

	C.psi.drotation_dxi = Vector(3).zeros();
	C.psi.drotation_dxj = Vector(3).zeros();
	C.psi.drotation_dxl = Vector(3).zeros();
	C.psi.drotation_dxt = Vector(3).zeros();
	
	
	circles.push_back(C);
	
	return circles;
	
}




unsigned int Tessellation::coverHemisphere2(TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, CircularInterfaceForm form){
	CircularInterface C;
	Vector v(3);
	Matrix NullMatrix(3,3);
	NullMatrix.zeros();
	
	v = tessellationAxis.v;
	
	C.lambdaRotation.rotation = M_PI * 0.5;
	C.lambda.rotation = M_PI * 0.5;
	C.psi.rotation=0;
	C.g=0;
	C.normal = v;
	C.form = form;
	C.id = circles.size();
	C.valid = true;
	C.index = -1;
	C.hasDerivatives=true;
	C.extended=false;
	
	C.dmu_dx = NullMatrix;
	C.lambda.drotation_dxi = Vector(3).zeros();
	C.lambda.drotation_dxj = Vector(3).zeros();
	C.lambda.drotation_dxl = Vector(3).zeros();
	C.lambda.drotation_dxt = Vector(3).zeros();

	C.psi.drotation_dxi = Vector(3).zeros();
	C.psi.drotation_dxj = Vector(3).zeros();
	C.psi.drotation_dxl = Vector(3).zeros();
	C.psi.drotation_dxt = Vector(3).zeros();
	
	
	circles.push_back(C);
	
	return circles.size()-1;
}


void Tessellation::addLimitingInterface(LimitingInterface &limit, CircularInterfacesPerAtom &circles){
	CircularInterface C;
	Vector v(3);
	Matrix NullMatrix(3,3);
	NullMatrix.zeros();
	
	C.lambdaRotation.rotation = acos(limit.g);
	C.lambda.rotation = C.lambdaRotation.rotation;
	C.psi.rotation=limit.psi;
	C.g=limit.g;
	C.normal = limit.v;
	if(limit.flip) C.form=CONCAVE;
	else C.form = CONVEX;
	C.id = circles.size();
	C.valid = true;
	C.index = -1;
	C.hasDerivatives=true;
	C.extended=true;
	
	C.dmu_dx = NullMatrix;
	C.lambda.drotation_dxi = Vector(3).zeros();
	C.lambda.drotation_dxj = Vector(3).zeros();
	C.lambda.drotation_dxl = Vector(3).zeros();
	C.lambda.drotation_dxt = Vector(3).zeros();

	C.psi.drotation_dxi = Vector(3).zeros();
	C.psi.drotation_dxj = Vector(3).zeros();
	C.psi.drotation_dxl = Vector(3).zeros();
	C.psi.drotation_dxt = Vector(3).zeros();
	
	
	circles.push_back(C);
	
}


void Tessellation::setupTessellationAxis(TessellationAxis &tessellationAxis, Hemisphere hemisphere, int closestNeighbour, vector<Vector> &atoms, vector<float> &radii, Vector &origin, float radius, DerivativeMode dmode, CircularInterfacesPerAtom &circles){
	float d;
	Vector mu(3);
	Matrix Identity(3,3);
	Matrix dchi_dx(3,3);
	Identity.eye();
	dchi_dx.zeros();
	tessellationAxis.v=Vector(3);
	tessellationAxis.auxiliary=Vector(3);
	tessellationAxis.planeNormal=Vector(3);
	float g;
	float r_i,r_l;
	
	if(dmode==GENERAL_POSITION){
		tessellationAxis.v(0) = 1.0;
		tessellationAxis.v(1) = 0.0;
		tessellationAxis.v(2) = 0.0;
		tessellationAxis.dchi_dx=dchi_dx;
		tessellationAxis.index=-2;
		
		tessellationAxis.auxiliary(0)=0;
		tessellationAxis.auxiliary(1)=1;
		tessellationAxis.auxiliary(2)=0;
		
		tessellationAxis.planeNormal(0)=0;
		tessellationAxis.planeNormal(1)=0;
		tessellationAxis.planeNormal(2)=1;		
	}
	tessellationAxis.hemisphere=hemisphere;
	tessellationAxis.mode=dmode;
	
	if(dmode==ALIGNED){
		tessellationAxis.index=closestNeighbour;
		Vector cn(3);
		cn=atoms[closestNeighbour];
		mu=calculateInterfaceNormal(origin, atoms[closestNeighbour], d);
		tessellationAxis.v=mu;
		tessellationAxis.dchi_dx= (1.0/d)*(Identity - kron(mu,mu.t())).t();
		
		//printf("TAX: %f %f %f, d: %d, origin: %f %f %f, CN: %f %f %f\n",mu(0),mu(1),mu(2),d,origin(0),origin(1),origin(2),cn(0),cn(1),cn(2));
		r_i = radii[closestNeighbour];
		r_l = radius;
		g = (d * d + r_l * r_l - r_i * r_i ) / (2 * d * r_l);

		if(g<0){
			tessellationAxis.v=-tessellationAxis.v;
		}
	}
	
	if(hemisphere==BACKHEMISPHERE){
		tessellationAxis.v=-tessellationAxis.v;
		tessellationAxis.planeNormal=-tessellationAxis.planeNormal;
		tessellationAxis.dchi_dx=-tessellationAxis.dchi_dx;
	}
	
	
	
}


void Tessellation::cleanCircularIntersections(CircularInterfacesPerAtom &circles){
	vector<CircularInterface>::iterator it;
	float angle;
	bool erased;
	
	it = circles.begin();
	while(it != circles.end()){
		erased=false;
		if(it->form!=SPLITTER && it->intersectionBranches.size()==0){
			it = circles.erase(it);
			erased=true;
		}
		
		if(!erased) ++it;
	}
	
	
}


void Tessellation::buildGaussBonnetPath(int i, vector<Vector> &atoms, vector<float> &radii, SASAs &sasas, bool split, vector<int> &neighbourlist, int closestNeighbour, bool useDepthBuffer){
	CircularInterfacesPerAtom circles,precircles;
	CircularInterfacesPerAtom circlesFrontHemisphere;
	CircularInterfacesPerAtom circlesBackHemisphere;
	//CircularInterface *splitter0, *splitter1;
	bool pass0, pass1;
	
	bool tessellationComplete;
	printf("BUILDING PATH %d\n",i);
	
	TessellationAxis frontTessellationAxis;
	TessellationAxis backTessellationAxis;
	Vector origin;
	float radius;
	bool isInsideAnotherSphere;
	Vector v;
	MultiLayeredDepthBuffer depthBuffer0;
	MultiLayeredDepthBuffer depthBuffer1;
	unsigned int globalSegmentCounter;
	SASASegmentList sasa;
	
	tessellationComplete=false;
	
	if(hasDepthBuffer && prevSasasForMolecule.size()==0) useDepthBuffer=false;
	
	if(hasDepthBuffer && useDepthBuffer){
		if(numbBuffer>1)
			depthBuffer0 = MultiLayeredDepthBuffer(depthData, occludedDistribution, exposedDistribution, 0);
		depthBuffer1 = MultiLayeredDepthBuffer(depthData, occludedDistribution, exposedDistribution, 1);
	}
	
	
	origin = atoms[i];
	torigin = origin;
	radius = radii[i];
	tradius = radius;
	ti = i;
	
	
	
	globalSegmentCounter=0;
	
	

	
	
	srand(2);
	
	
	
	sasas.push_back(sasa);
	
	precircles.clear();
	circles.clear();
	precircles.reserve(neighbourlist.size()+1);
	circles.reserve(neighbourlist.size()+1);
	isInsideAnotherSphere = makeCircularInterfaces(i,origin, radius, atoms, radii, precircles, neighbourlist);
	
	if(isInsideAnotherSphere) return;
	
	
	for(unsigned int j=0;j<precircles.size();j++){
		determineProjection(origin, radius, precircles[j]);
	}
	
	totalInterfaces+=precircles.size();
	
	
	
	while(!tessellationComplete){
	
		
		if(hasDepthBuffer && useDepthBuffer){
			for(unsigned int j=0;j<precircles.size();j++){
				if(numbBuffer>1)
					depthBuffer0.addSphere(precircles[j].normal, precircles[j].lambda.rotation, precircles[j].form!=CONVEX, precircles[j].kappa[0], precircles[j].psi2[0]);
				depthBuffer1.addSphere(precircles[j].normal, precircles[j].lambda.rotation, precircles[j].form!=CONVEX, precircles[j].kappa[1], precircles[j].psi2[1]);
			}

			pass0=false;
			for(unsigned int j=0;j<precircles.size();j++){
				if(numbBuffer>1)
					pass0=depthBuffer0.passesBuffer(precircles[j].kappa[0], precircles[j].psi2[0], precircles[j].lambda.rotation, precircles[j].form!=CONVEX, precircles[j].exposedVectors);
				
				pass1=depthBuffer1.passesBuffer(precircles[j].kappa[1], precircles[j].psi2[1], precircles[j].lambda.rotation, precircles[j].form!=CONVEX, precircles[j].exposedVectors);
				
				
				if(pass0 || pass1){
					circles.push_back(precircles[j]);
				}
			
				
			}
		}
		else circles=precircles;
		
	// 	if(i==0){
	// 		FILE* file;
	// 		file = fopen ("segmentRays.csv","w");
	// 		
	// 		for(unsigned int j=0;j<precircles.size();j++){
	// 			for(unsigned int m=0; m<precircles[j].exposedVectors.size(); ++m){
	// 				fprintf(file,"%f %f %f 0\n", precircles[j].exposedVectors[m](0), precircles[j].exposedVectors[m](1), precircles[j].exposedVectors[m](2));
	// 			}
	// 		}
	// 	}
		
		
		benchmark.addQuantity("interfaces found",precircles.size());
		benchmark.addQuantity("interfaces kept",circles.size());
		
		
		if(circles.size()==0) return;
		
		depthBufferEstimatedArea=-1;
		if(hasDepthBuffer && useDepthBuffer){
			if(numbBuffer>1){
				depthBufferEstimatedArea=(depthBuffer0.exposedArea()+depthBuffer1.exposedArea())/2.0;
			}
			else{
				depthBufferEstimatedArea=depthBuffer1.exposedArea();
			}

			//this will delete very small areas. It is not necessary for correctness of the method, but speeds things up a bit more
			if(depthBufferEstimatedArea<0.01) return;
			
		}
		
		
		survivedInterfaces+=circles.size();
		
		
		//frontTessellationAxis = generalPosition(circles);
		//DerivativeMode dmode=ALIGNED;
		DerivativeMode dmode=GENERAL_POSITION;
		setupTessellationAxis(frontTessellationAxis, FRONTHEMISPHERE, closestNeighbour, atoms, radii, atoms[i], radii[i], dmode, circles);
		setupTessellationAxis(backTessellationAxis, BACKHEMISPHERE, closestNeighbour, atoms, radii, atoms[i], radii[i], dmode, circles);
			
		split=false;
		
		//there is room for optimisation here...
		if(split){
		
/*			circlesFrontHemisphere = coverHemisphere(frontTessellationAxis, circles, SPLITTER);
			circlesBackHemisphere = coverHemisphere(backTessellationAxis, circles, SPLITTER);
			
			
			
			splitter0 = &circlesFrontHemisphere[circlesFrontHemisphere.size()-1];
			splitter1 = &circlesBackHemisphere[circlesBackHemisphere.size()-1];
			
			if(hasDepthBuffer && useDepthBuffer){
				//we use only the buffer which is aligned with the auxiliary axis to maximize its scanning performance
				depthBuffer1.getSplitterExposedVectors(splitter0->exposedVectors);
				depthBuffer1.getSplitterExposedVectors(splitter1->exposedVectors);
			}
			
			//this will remove interaces that are completely inside other interfaces
			filterCircularInterfaces(circlesFrontHemisphere, false);
			filterCircularInterfaces(circlesBackHemisphere, false);
			
			
			reindexCircularInterfaces(circlesFrontHemisphere);
			reindexCircularInterfaces(circlesBackHemisphere);
			
			determineCircularIntersections(circlesFrontHemisphere);
			determineCircularIntersections(circlesBackHemisphere);
			
			bool tessellationComplete0, tessellationComplete1;
			//when we parallelize triforce, globalsegmentpointer needs to be made parallel as well somehow
			tessellationComplete0=buildIntersectionGraph(i, radius, frontTessellationAxis, circlesFrontHemisphere, sasas[sasas.size()-1], FRONTHEMISPHERE, string("gbonnet0.csv"), depthBuffer0, depthBuffer1, useDepthBuffer, split, globalSegmentCounter,true);
			tessellationComplete1=buildIntersectionGraph(i, radius, backTessellationAxis, circlesBackHemisphere, sasas[sasas.size()-1], BACKHEMISPHERE, string("gbonnet1.csv"), depthBuffer0, depthBuffer1, useDepthBuffer, split, globalSegmentCounter,true);
			
			tessellationComplete=tessellationComplete0 && tessellationComplete1;
			if(!tessellationComplete){
				sasas[sasas.size()-1].clear();
				useDepthBuffer=false;
			}
*/			
		}
		else{
			filterCircularInterfaces(circles, false);
			reindexCircularInterfaces(circles);
			determineCircularIntersections(circles, false);
			buildIntersectionGraphFirstPass(i, radius, frontTessellationAxis, circles);
			cleanCircularIntersections(circles);
			reindexCircularInterfaces2(circles);
			determinePsiRotations(frontTessellationAxis, circles);
			

printf("____________________________\n INITIAL ROUND\n");			
for(unsigned int j=0; j<circles.size(); ++j){
	printf("CIRCLE[%d] (%d)\n",circles[j].index,circles[j].form);
	for(IntersectionBranches::iterator it=circles[j].intersectionBranches.begin(); it!=circles[j].intersectionBranches.end(); ++it){
		if(it->second->direction==IN)
			printf(" -> %d [IN] %f\n", circles[it->second->id].index,it->first);
		else
			printf(" -> %d [OUT] %f\n", circles[it->second->id].index, it->first);
	}
}			
	
			
			//sortGaussBonnetPaths(i, radius, frontTessellationAxis, circles, sasas[sasas.size()-1], FRONTHEMISPHERE, string("gbonnet0.csv"), depthBuffer0, depthBuffer1, useDepthBuffer, split, globalSegmentCounter,true);
			
			
			copyIntersectionGraph(i, radius, backTessellationAxis, circles, circlesBackHemisphere);
			
			determinePsiRotations(backTessellationAxis, circlesBackHemisphere);
			
			
			
			
			
			coverHemisphere2(frontTessellationAxis, circles, SPLITTER);
			coverHemisphere2(backTessellationAxis, circlesBackHemisphere, SPLITTER);
			
			
			
			
			filterCircularInterfaces(circles, true);
			filterCircularInterfaces(circlesBackHemisphere, true);
			
			cleanCircularIntersections(circles);
			cleanCircularIntersections(circlesBackHemisphere);

			
			
			reindexCircularInterfaces2(circles);
			reindexCircularInterfaces2(circlesBackHemisphere);
			
			
			
			
			determineCircularIntersections(circles, true);
			determineCircularIntersections(circlesBackHemisphere, true);
			
			
			
			
			
			buildIntersectionGraphSplitterPass(i, radius, frontTessellationAxis, circles);
			splitterSanityCheck(circles);
			
			
printf("____________________________\nFRONT ROUND\n");			
for(unsigned int j=0; j<circles.size(); ++j){
	printf("CIRCLE[%d]\n",circles[j].index);
	for(IntersectionBranches::iterator it=circles[j].intersectionBranches.begin(); it!=circles[j].intersectionBranches.end(); ++it){
		if(it->second->direction==IN)
			printf(" -> %d [IN] %f\n", circles[it->second->id].index, it->first);
		else
			printf(" -> %d [OUT] %f\n", circles[it->second->id].index, it->first);
	}
}			
			
			
			buildIntersectionGraphSplitterPass(i, radius, backTessellationAxis, circlesBackHemisphere);
			splitterSanityCheck(circlesBackHemisphere);


printf("____________________________\nBACKFLIPPED ROUND\n");			
for(unsigned int j=0; j<circlesBackHemisphere.size(); ++j){
	printf("CIRCLE[%d]\n",circlesBackHemisphere[j].index);
	for(IntersectionBranches::iterator it=circlesBackHemisphere[j].intersectionBranches.begin(); it!=circlesBackHemisphere[j].intersectionBranches.end(); ++it){
		if(it->second->direction==IN)
			printf(" -> %d [IN] %f\n", circlesBackHemisphere[it->second->id].index,it->first);
		else
			printf(" -> %d [OUT] %f\n", circlesBackHemisphere[it->second->id].index,it->first);
	}
}			
			
			
			globalSegmentCounter=-2;

			
			buildIntersectionGraphArtificialPointsPass(i, radius, frontTessellationAxis, circles, sasas[sasas.size()-1], FRONTHEMISPHERE, globalSegmentCounter);
			buildIntersectionGraphArtificialPointsPass(i, radius, backTessellationAxis, circlesBackHemisphere, sasas[sasas.size()-1], BACKHEMISPHERE, globalSegmentCounter);
			

			//just to make sure there are no lonely branches floating around...
			reindexCircularInterfaces3(circles);
			reindexCircularInterfaces3(circlesBackHemisphere);

			
			
			buildIntersectionGraphCollectionPass(i, radius, frontTessellationAxis, circles, sasas[sasas.size()-1], FRONTHEMISPHERE, string("gbonnet0.csv"), depthBuffer0, depthBuffer1, useDepthBuffer, split, globalSegmentCounter,true);
			buildIntersectionGraphCollectionPass(i, radius, backTessellationAxis, circlesBackHemisphere, sasas[sasas.size()-1], BACKHEMISPHERE, string("gbonnet0.csv"), depthBuffer0, depthBuffer1, useDepthBuffer, split, globalSegmentCounter,true);

tessellationComplete=true;			
			
/*			
			filterCircularInterfaces(frontTessellationAxis, radius, circles);
			reindexCircularInterfaces(circles);
			
			determineCircularIntersections(circles);
			
			tessellationComplete=buildIntersectionGraph(i, radius, frontTessellationAxis, circles, sasas[sasas.size()-1], FRONTHEMISPHERE, string("gbonnet0.csv"), depthBuffer0, depthBuffer1, useDepthBuffer, split, globalSegmentCounter,false);
			splitTessellation(sasas[sasas.size()-1],frontTessellationAxis,backTessellationAxis, circles);
			addDerivatives(sasas[sasas.size()-1],circles, frontTessellationAxis, backTessellationAxis);
*/
		}
	}
	
/*	
	//quick print
	for(unsigned int j=0; j<sasas[sasas.size()-1].size(); ++j){
		SASASegment seg;
		seg=sasas[sasas.size()-1][j];
		
		
		
		printf("SASA: %d %d,%d,%d [%d,%d,%d] (%d) {%d} psi:%f lmbd:%f PHI0:%f PHI1:%f\n",i,seg.id0,seg.id1,seg.id2,seg.index0,seg.index1,seg.index2,seg.tessellationAxis.hemisphere, seg.form1, seg.psi.rotation,seg.lambda.rotation,seg.rotation0.rotation, seg.rotation1.rotation);
		
		
	}
*/	
	
}






Vector Tessellation::generalPosition(CircularInterfacesPerAtom &circles){
	Vector chi(3),y(3);
	chi(0)=1;
	chi(1)=0;
	chi(2)=0;
	y(0)=0;
	y(1)=1;
	y(2)=0;
	
	float d0;
	Vector ni;
	
	Vector n(3);
	n=Vector(3).zeros();
	for(unsigned int i=0; i<circles.size(); ++i){
		n+=circles[i].normal;
	}
	n=n/(float)circles.size();
	n=-n;
	chi=normalise(n);
	
	
	srand(50);
	bool degenerate;
	do{
		degenerate=false;
		for(unsigned int i=0; i<circles.size() && !degenerate; ++i){
			d0=abs(1-abs(norm_dot(chi,circles[i].normal)));
			if(d0<=THRESHOLD_GENERAL_POSITION){
				degenerate=true;
			}
		}
		if(degenerate){
			chi = randu<fvec>(3);
			chi(0) = chi(0)-0.5;
			chi(1) = chi(1)-0.5;
			chi(2) = chi(2)-0.5;
			chi=normalise(chi);
		}
	}while(degenerate);
	
	
	
	
	
	return chi;
}






SASAs &Tessellation::sasas(){
	return sasasForMolecule;
}



float Tessellation::vsign(float v)
{
	if(v>=0) return 1.0;
	else return -1.0;
}

float Tessellation::cot(float a){
	return 1.0/tan(a);
}

float Tessellation::csc(float a){
	return 1.0/sin(a);
}




float Tessellation::getAngleBetweenNormals(Vector &a, Vector &b){
	return acos(dot(a,b));
}


float Tessellation::getAngle(Vector &a, Vector &b){
	return acos(norm_dot(a,b));
}


bool Tessellation::isZero(float v){
	if(fabs(v) <= THRESHOLD_NUMERICAL) return true;
	else return false;
}

bool Tessellation::isInPositiveEpsilonRange(float v, float eps){
	if(eps-(v+THRESHOLD_NUMERICAL) <= 0) return true;
	else return false;
}

bool Tessellation::isWithinNumericalLimits(float x, float l){
	if(abs(x-l)<=THRESHOLD_NUMERICAL) return true;
	else return false;
}


bool Tessellation::isWithinStrongNumericalLimits(float x, float l){
	if(abs(x-l)<=THRESHOLD_STRONG_NUMERICAL) return true;
	else return false;
}


LambdaRotation Tessellation::calculateLambda(float d_i, float r_l, float r_i, Vector &mu_i){
	CircularInterfaceForm form;
	return calculateLambda(0,d_i, r_l, r_i, mu_i, form);
}


LambdaRotation Tessellation::calculateLambda(int index_i, float d_i, float r_l, float r_i, Vector &mu_i, CircularInterfaceForm &form){
	float g;
	LambdaRotation r;
	Vector mu_i_original;
	Vector dlambda_dxi(3), dlambda_dxl(3);
	Vector dd_dxi;
	
	mu_i_original = mu_i;
	g = (d_i * d_i + r_l * r_l - r_i * r_i ) / (2 * d_i * r_l);

	if(g<0){
		form=CONCAVE;
		g = abs(g);
		mu_i=-mu_i;
	}
	else form=CONVEX;
	
	r.rotation = acos(g);
	r.d_i = d_i;
	r.r_l = r_l;
	r.g = g;
	
	return r;
	
}



Rotation Tessellation::calculateLambdaDerivatives(LambdaRotation &r, CircularInterface &circle){
	float g;
	float dg_dd;
	float dlambda_dg;
	float q;
	Vector mu_i_original;
	Vector dlambda_dxi(3), dlambda_dxl(3);
	Vector dd_dxi;
	Rotation lambda;
	float d_i;
	float r_l;
	Vector mu_i;
	CircularInterfaceForm form;
	
	g = r.g;
	d_i = r.d_i;
	r_l = r.r_l;
	mu_i = circle.normal;
	form = circle.form;

	if(g > 1.0-MINISCULE) g = 1.0-MINISCULE;
	if(g < -1.0+MINISCULE) g = -1.0+MINISCULE;
	

	dlambda_dg = -1/sqrt(1-g*g);
	
	
	
	
	if(form==CONVEX){
		dg_dd = -g/d_i + 1/r_l;
		q = 1;
	}
	else{
		dg_dd = -(g/d_i + 1/r_l);
		q = -1;
	}
	

	
	dd_dxi = q*mu_i;
	
	
	dlambda_dxi = dlambda_dg * dg_dd * dd_dxi;
	//smoothing
	//if(norm(dlambda_dxi,2)>1.0) dlambda_dxi = normalise(dlambda_dxi);
	dlambda_dxl = -dlambda_dxi;
	
	lambda.rotation = r.rotation;
	lambda.drotation_dxi = dlambda_dxi;
	lambda.drotation_dxj = Vector(3).zeros();
	lambda.drotation_dxl = dlambda_dxl;
	lambda.drotation_dxt = Vector(3).zeros();
	
	
	
	
	
	return lambda;
	
	
}



void Tessellation::determineProjection(Vector &origin, float r_l, CircularInterface &circle){
	circle.lambdaRotation = calculateLambda(circle.index, circle.d, r_l, circle.sphereRadius, circle.normal, circle.form);
	circle.lambda.rotation = circle.lambdaRotation.rotation;
	circle.valid=true;
}






void Tessellation::calculateProjectionAndDerivatives(TessellationAxis &tessellationAxis, CircularInterface &circle){
	Matrix Identity(3,3);
	Matrix dmu_dx(3,3);
	Matrix fd_dmu_dx(3,3);
	Identity.eye();
	
	
	
	

	circle.lambda = calculateLambdaDerivatives(circle.lambdaRotation, circle);
	
	//derivatives
	if(circle.form==CONVEX)
		dmu_dx = (1.0/circle.d)*(Identity - kron(circle.normal,circle.normal.t())).t();
	else
		dmu_dx = -(1.0/circle.d)*(Identity - kron(circle.normal,circle.normal.t())).t();
	
	circle.dmu_dx = dmu_dx;
	
	
	
	
	circle.psi = calculatePsi(tessellationAxis, circle);
	
	circle.hasDerivatives=true;
	
	
	
}



void Tessellation::determinePsiRotations(TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles){
	for(unsigned i=0; i < circles.size(); ++i)
		circles[i].psi.rotation = getAngle(tessellationAxis.v,circles[i].normal);
}


IntersectionPair Tessellation::determineIntersectionPoints(float radius, CircularInterface &K, CircularInterface &J){
	IntersectionPair res;
	float g_k, g_j;
	float r_i;
	Vector mu_k, mu_j;
	float phi_kj, phi_jk;
	float tau_kj,tau_jk;
	Vector tmp1, tmp2;
	Vector eta_kj, eta_jk;
	Vector omega_kj, omega_jk;
	Vector p_kj, p_jk;
	float gamma_kj,gamma_jk;





	g_k = K.g;
	g_j = J.g;

	r_i = radius;


	mu_k = K.normal;
	mu_j = J.normal;

	phi_kj = acos(dot(mu_k,mu_j));
	phi_jk = acos(dot(mu_j,mu_k));

	tau_kj = (g_k - g_j * cos(phi_kj)) / (pow(sin(phi_kj),2));
	tau_jk = (g_j - g_k * cos(phi_jk)) / (pow(sin(phi_jk),2));


	eta_jk = eta_kj = (mu_k * tau_kj) + (mu_j * tau_jk);

	omega_kj = cross(mu_k, mu_j) / sin(phi_kj);
	omega_jk = cross(mu_j, mu_k) / sin(phi_jk);

	gamma_kj = sqrt(r_i * r_i - g_k*tau_kj - g_j*tau_jk);
	gamma_jk = sqrt(r_i * r_i - g_j*tau_jk - g_k*tau_kj);

	p_kj = (omega_kj * gamma_kj) + eta_kj;
	p_jk = (omega_jk * gamma_jk) + eta_jk;
	
	
	res.k_j = p_kj;
	res.j_k = p_jk;


	return res;
	

	
}


Vector Tessellation::calculateInterfaceNormal(const Vector &v_l, const Vector &v_i){
	return calculateInterfaceNormal(v_l,v_i);
}


Vector Tessellation::calculateInterfaceNormal(const Vector &v_l, const Vector &v_i, float &d){
	Vector v;
	v = v_i - v_l;
	d = norm(v,2);
	v = v/d;
	return v;
	
}


bool Tessellation::makeCircularInterfaces(int l,Vector &origin, float r_l, vector<fvec> &atoms, vector<float> &radii, vector<CircularInterface> &circles, vector<int> &neighbourlist){
	CircularInterface circle;
	float r_i;
	float d_i;
	Vector mu_i;
	vector<int>::iterator it;
	int j;
	

	for(it=neighbourlist.begin(); it!=neighbourlist.end(); ++it){
		j=*it;
		if(l != j){
			r_i = radii[j];
			
			mu_i = calculateInterfaceNormal(origin, atoms[j], d_i);
			
			//reject, if no intersection
			if(d_i < r_l + r_i && d_i+r_i > r_l && d_i+r_l > r_i){
				circle.id = circles.size();
				circle.normal = mu_i;
				circle.d = d_i;
				circle.sphereRadius = r_i;
				circle.intersect = false;
				circle.index=j;
				circle.hasDerivatives=false;
				circle.extended=false;
				circle.erased=false;
				circle.hemisphere=FRONTHEMISPHERE;
				circles.push_back(circle);
			}
			//in this case, the atom is completely inside another atom and has no SASA
			else if(d_i+r_l <= r_i){
				circles.clear();
				return true;
			}
			
		}
	}
	return false;
}



void Tessellation::depleteCircularInterfaces(TessellationAxis &tessellationAxis, float radius, vector<CircularInterface> &circles){
	circles.clear();
	TessellationAxis t;
	t=tessellationAxis;
	t.v = -tessellationAxis.v;
	circles = coverHemisphere(t, circles, CONCAVE);
}


/*
bool Tessellation::isOccludedBy(CircularInterface &a, CircularInterface &splitter){
	float angle;
	Vector n0,n1;
	
	n0 = a.normal;
	n1 = splitter.normal;
	
	angle = getAngle(n0,n1);
	
	if(a.form==CONVEX){
		//convex circle IT is inside of convex circle i
		if(splitter.form == CONVEX){
			
			if(a.lambda.rotation + angle-THRESHOLD_INTERFACE < splitter.lambda.rotation){
				return true;
			}
		}
		//convex circle is outside of concave circle i
		else{
			if(angle+THRESHOLD_INTERFACE-a.lambda.rotation >  splitter.lambda.rotation){
				return true;
			}
		}
		
	}
	else{
		//concave circle IT has a free area. This area is covered by convex circle i
		if(splitter.form == CONVEX){
			if(a.lambda.rotation + angle-THRESHOLD_INTERFACE < splitter.lambda.rotation){
				return true;
			}
			
		}
		else{
			//concave circle IT is completely inside the occlusion of concave circle i
			if(a.lambda.rotation > splitter.lambda.rotation + angle-THRESHOLD_INTERFACE){
				return true;
			}
			//concave circle IT has a free area. This area is covered by concave circle i
			else if(angle+THRESHOLD_INTERFACE > a.lambda.rotation + splitter.lambda.rotation){
				return true;
				
			}
			
		}

	}
	return false;
}
*/

int Tessellation::filterCircularInterfaces(vector<CircularInterface> &circles, bool splitterOnly){
	vector<CircularInterface>::iterator it;
	float angle;
	bool erased;
	Vector n0,n1;
	//sunsigned int c;
	vector<CircularInterface> circles2;
	
	//c=0;
	it = circles.begin();
	while(it != circles.end()){
		erased=false;
	
		
		n0 = it->normal;
		
		for(unsigned i=0;i<circles.size();i++){
			if(it->id != circles[i].id){// && !circles[i].erased){
				if(!splitterOnly || (splitterOnly && circles[i].form==SPLITTER)){
				
					n1 = circles[i].normal;
					
					angle = getAngle(n0,n1);
					
					if(it->form==CONVEX){
						//convex circle IT is inside of convex circle i
						if(circles[i].form == CONVEX){
							
							if(it->lambda.rotation + angle-THRESHOLD_INTERFACE < circles[i].lambda.rotation){
								it = circles.erase(it);
								it->erased=true;
								erased=true;
								break;
							}
						}
						//convex circle is outside of concave circle i
						else{
							if(angle+THRESHOLD_INTERFACE-it->lambda.rotation >  circles[i].lambda.rotation){
								it = circles.erase(it);
								it->erased=true;
								erased=true;
								break;
							}
						}
						
					}
					else{
						//concave circle IT has a free area. This area is covered by convex circle i
						if(circles[i].form == CONVEX){
							if(it->lambda.rotation + angle-THRESHOLD_INTERFACE < circles[i].lambda.rotation){
								circles.clear();
								return -1;
							}
							
						}
						else{
							//concave circle IT is completely inside the occlusion of concave circle i
							if(it->lambda.rotation > circles[i].lambda.rotation + angle-THRESHOLD_INTERFACE){
								it = circles.erase(it);
								it->erased=true;
								erased=true;
								break;
							}
							//concave circle IT has a free area. This area is covered by concave circle i
							else if(angle+THRESHOLD_INTERFACE > it->lambda.rotation + circles[i].lambda.rotation){
								circles.clear();
								return -1;
								
							}
							
						}

					}	
				}
				
						
					
			}
		}
		//if(!erased) ++c;
		//++it;
		if(!erased) ++it;
	}
	
	
	//if just the splitter is left, the hemisphere is completely free. In that case, the splitter is unnecessary
	//if(circles.size()==1 && circles[0].form==SPLITTER)
	//	circles.clear();
	
	
	/*
	circles2.reserve(c);
	for(unsigned int i=0; i<circles.size(); ++i){
		if(!circles[i].erased) circles2.push_back(circles[i]);
	}
	
	circles=circles2;
	*/
	
	
	return 0;

}




void Tessellation::outputGaussBonnetPath(SASAs &points){
	/*
	SASASegmentList::iterator it;
	int i;
	for(it=points.sasa.begin(), i=0; it!=points.sasa.end(); ++it, ++i)
		fprintf(stderr,"GBPATH[%d] %d - %d indexes: %d - %d form: %d %d\n", i, it->id0, it->id1, it->index0, it->index1, it->form0, it->form1);
	*/
}



void Tessellation::reindexCircularInterfaces(CircularInterfacesPerAtom &circles){
	for(unsigned i=0;i<circles.size();i++){
		circles[i].id = i;
	}
}

void Tessellation::reindexCircularInterfaces2(CircularInterfacesPerAtom &circles){
	map<int,int> cross;
	IntersectionBranches::iterator it,it2;
	CircularIntersections::iterator it3;
	bool erased;
	
	
	
	
	for(unsigned i=0;i<circles.size();i++){
		cross[circles[i].id]=i;
	}
	
	for(unsigned i=0;i<circles.size();i++){
		it3=circles[i].circularIntersections.begin();
		while(it3!=circles[i].circularIntersections.end()){
			erased=false;
			if(cross.find(it3->second.id)!=cross.end()){
				it3->second.id=cross[it3->second.id];
			}
			else{
				it3=circles[i].circularIntersections.erase(it3);
				erased=true;
			}
			if(!erased) it3++;
		}
	}
	
	
	for(unsigned i=0;i<circles.size();i++){
		circles[i].id = i;
		it = circles[i].intersectionBranches.begin();
		
		//if the branch is pointing to a not-deleted circle, we update its id, or remove it in the other case
		while(it != circles[i].intersectionBranches.end()){
			erased=false;
			if(cross.find(it->second->id)!=cross.end()){
				it->second->id=cross[it->second->id];
				it->second->PHI.omega.id_i=cross[it->second->PHI.omega.id_i];
				it->second->PHI.omega.id_j=cross[it->second->PHI.omega.id_j];
				it->second->PHI.eta.id_i=cross[it->second->PHI.eta.id_i];
				it->second->PHI.eta.id_j=cross[it->second->PHI.eta.id_j];
			}
			else{
				//it's not good to delete the intersection points, but we have to give them a different id
				it->second->id=i;
				it->second->flagged=true;
			}
			if(!erased) ++it;
		}
	}
}

void Tessellation::reindexCircularInterfaces3(CircularInterfacesPerAtom &circles){
	map<int,int> cross;
	IntersectionBranches::iterator it,it2;
	bool erased;
	
	
	for(unsigned i=0;i<circles.size();i++){
		cross[circles[i].id]=i;
	}
	
	for(unsigned i=0;i<circles.size();i++){
		it = circles[i].intersectionBranches.begin();
		
		//if the branch is pointing to a not-deleted circle, we update its id, or remove it in the other case
		while(it != circles[i].intersectionBranches.end()){
			erased=false;
			if(cross.find(it->second->id)!=cross.end()){
			}
			else{
				it2=it;
				++it2;
				circles[i].intersectionBranches.erase(it);
				erased=true;
				it=it2;
			}
			if(!erased) ++it;
		}
	}
}



void Tessellation::insertArtificialIntersectionPoints(CircularInterface &I, TessellationAxis &tessellationAxis, Hemisphere hemisphere, SASASegmentList &sasa, unsigned int &globalSegmentCounter){
	Rotation f,f_reverse;
	SASASegment sasaSegment;
	
	
	
	
	//this will create two points on the border of the circular region on opposite sides.

	
	f.rotation = -M_PI/2;
	f.drotation_dxi=Vector(3).zeros();
	f.drotation_dxj=Vector(3).zeros();
	f.drotation_dxl=Vector(3).zeros();
	f.drotation_dxt=Vector(3).zeros();

	f_reverse.rotation = M_PI/2;
	f_reverse.drotation_dxi=Vector(3).zeros();
	f_reverse.drotation_dxj=Vector(3).zeros();
	f_reverse.drotation_dxl=Vector(3).zeros();
	f_reverse.drotation_dxt=Vector(3).zeros();
	
	
	
	
	calculateProjectionAndDerivatives(tessellationAxis, I);
	
	sasaSegment.i = globalSegmentCounter;
	sasaSegment.i2 = 0;
	sasaSegment.artificial=true;
	
	//this will create two improper intersection-points, but in the end, it's the same result as with proper ones..
	sasaSegment.v0=I.normal;
	sasaSegment.v1=I.normal;
	sasaSegment.weight=1;
	sasaSegment.weight1=1;
	
	sasaSegment.id0 = I.id;
	sasaSegment.id1 = I.id;
	sasaSegment.id2 = I.id;
	
	sasaSegment.index0 = I.index;
	sasaSegment.index1 = I.index;
	sasaSegment.index2 = I.index;
	sasaSegment.lambda = I.lambda;
	
	
	sasaSegment.psi = I.psi;
	sasaSegment.normalForCircularInterface = I.normal;
	sasaSegment.form0 = I.form;
	sasaSegment.form1 = I.form;
	sasaSegment.form2 = I.form;
	
	sasaSegment.tessellationAxis = tessellationAxis;
	sasaSegment.hemisphere = hemisphere;
	
	
	sasaSegment.rotation0 = f;
	sasaSegment.rotation1 = f_reverse;
	sasaSegment.depthBufferEstimatedArea=depthBufferEstimatedArea;
	
	sasa.push_back(sasaSegment);
	
	sasaSegment.i2 = 1;
	sasaSegment.rotation1 = f;
	sasaSegment.rotation0 = f_reverse;
	
	sasa.push_back(sasaSegment);
	
	globalSegmentCounter--;
	
			
}




/*


void Tessellation::insertArtificialIntersectionPoints(CircularInterface &I, IntersectionGraph &intersectionGraph, TessellationAxis &tessellationAxis){
	Vector v,v2,o,o2;
	int k;
	IntersectionNode a,b;
	IntersectionAddress addressA,addressB;
	Vector x0,x1;
	Rotation f,f_reverse;
	int q=1;
	
	
	
	//this will create two points on the border of the circular region on opposite sides.
	//if(tessellationAxis(0)<0) q=-q;
	//if(I.form==CONVEX) q=-q;

	
	if(q<0){
		f.rotation = M_PI/2;
		f_reverse.rotation = -M_PI/2;
	}
	else{
		f.rotation = -M_PI/2;
		f_reverse.rotation = M_PI/2;
	}
	
	a.id0 = addressA.id0 = I.id;
	a.id1 = addressA.id1 = -I.id;
	a.pointsTo0 = a.prev0 = -I.id;
	a.pointsTo1 = a.prev1 = I.id;
	a.vector = x0;
	a.rotation0=f;
	a.rotation1=f;
	a.visited=false;
	a.rotation0.drotation_dxi=Vector(3).zeros();
	a.rotation0.drotation_dxj=Vector(3).zeros();
	a.rotation0.drotation_dxl=Vector(3).zeros();
	a.rotation1.drotation_dxi=Vector(3).zeros();
	a.rotation1.drotation_dxj=Vector(3).zeros();
	a.rotation1.drotation_dxl=Vector(3).zeros();

	b.id0 = addressB.id0 = -I.id;
	b.id1 = addressB.id1 = I.id;
	b.pointsTo0 = b.prev0 = I.id;
	b.pointsTo1 = b.prev1 = -I.id;
	b.vector = x1;
	b.rotation0=f_reverse;
	b.rotation1=f_reverse;
	b.visited=false;
	b.rotation0.drotation_dxi=Vector(3).zeros();
	b.rotation0.drotation_dxj=Vector(3).zeros();
	b.rotation0.drotation_dxl=Vector(3).zeros();
	b.rotation1.drotation_dxi=Vector(3).zeros();
	b.rotation1.drotation_dxj=Vector(3).zeros();
	b.rotation1.drotation_dxl=Vector(3).zeros();
	
	
	intersectionGraph[addressA]=a;
	intersectionGraph[addressB]=b;
	
			
}
*/





int Tessellation::sgn(float d){
	if(d>=0) return 1;
	else return -1;
}
/*
void Tessellation::doIntersect(CircularInterface &k, CircularInterface &j){
	float angle;
	RhoContainer c;
	
	angle = getAngle(k.normal,j.normal);
	if(angle+THRESHOLD_INTERFACE < k.lambda.rotation + j.lambda.rotation){
		if(angle-THRESHOLD_INTERFACE + k.lambda.rotation > j.lambda.rotation && angle-THRESHOLD_INTERFACE + j.lambda.rotation > k.lambda.rotation){
			return true
		}
	}
	return false;
	
}
*/

void Tessellation::determineCircularIntersections(CircularInterfacesPerAtom &circles, bool splitterOnly){
	float angle;
	RhoContainer c;
	unsigned int start;
	//this doesn't really need to be optimised. The reason is, that the number of excluded interface-interface intersections is much smaller and at most equal to the number of included intersections
	for(unsigned k=0;k<circles.size();k++){
		start=k+1;
		if(splitterOnly && circles[k].form!=SPLITTER) continue;
		if(splitterOnly) start=0;
		for(unsigned j=start;j<circles.size();j++){
			if(k!=j){
				//determine if there will be intersections
				angle = getAngle(circles[k].normal,circles[j].normal);
				if(angle+THRESHOLD_INTERFACE < circles[k].lambda.rotation + circles[j].lambda.rotation){
					if(angle-THRESHOLD_INTERFACE + circles[k].lambda.rotation > circles[j].lambda.rotation && angle-THRESHOLD_INTERFACE + circles[j].lambda.rotation > circles[k].lambda.rotation){
						c.rho=angle;
						c.id = k;
						circles[j].circularIntersections[circles[k].index]=c;
						//circles[j].circularIntersections[k]=c;
						c.id = j;
						circles[k].circularIntersections[circles[j].index]=c;
						//circles[k].circularIntersections[j]=c;
					}
				}
			}
		}
	}
		
	
}




/**
 * n is a normal vector to the plane connecting integration origin, origin and center of circular region
 * o is a vector pointing in direction of integration origin but intersects with the interface of circular region
 * a is one of the PHI vectors from the center of circular region to its interface
 * 
 * function will return the angle between a and the plane between 0 and PI (instead of 0 and PI/2)
 * 
 */


float Tessellation::complLongAngle(Vector &vi, Vector &vj, Vector &vk){
	float v;
	int s;
	Vector nij(3), nik(3), vij(3), vik(3);
	Vector ni(3);
	
	
	vij = vj - vi;
	nij = cross(vi, vij);
	nij = nij/norm(nij,2);
	vik = vk - vi;
	nik = cross(vi, vik);
	nik = nik/norm(nik,2);
	
	ni = cross(vi,nij);
	
	
	
	v = acos(dot(nij,nik));
	s = sgn(dot(ni,nik));
	
	if(s>0) v = -v;
	
	//addition
	if(v<0) v = 2*M_PI+v;
	   
	return v;
		
}




/*
 * 
[7] (7,1) d: 0.725943 (OUT)
[7] (5,7) d: 0.769139 (IN)
[7] (7,3) d: 1.506200 (OUT)
[7] (7,4) d: 2.223748 (OUT)
[7] (4,7) d: 2.265115 (IN)
[7] (3,7) d: 3.938589 (IN)
[7] (1,7) d: 4.631484 (IN)
[7] (7,5) d: 4.641597 (OUT)

 * 
 * */

/*
float Tessellation::angularInterface(Vector &x0, Vector &v, Vector &p0, Vector &p1){
	Vector vx(3);
	Vector vp0(3);
	Vector vp1(3);
	float eta;
	int s;
	
	
	
	
	vx = x0-v;
	vp0 = p0-v;
	vp1 = p1-v;

	eta = getAngle(vp1,vx);
	
	
	s = sgn(dot(vp0,vx));
	
	if(s > 0) eta = -eta;
	
	return eta;
	
}

void Tessellation::measurementPoints(Vector &p0, Vector &p1, TessellationAxis &tessellationAxis, CircularInterface &I){
	Vector v(3);
	Vector o(3);
	Vector vx(3);
	Vector s(3), s2(3);
	
	o = tessellationAxis;
	v = I.normal * I.g;
	
	
	
	if(isWithinNumericalLimits(dot(o,I.normal),1.0) ){
		p1 = Vector(3);
		p1(0) = 0;
		p1(1) = -1;
		p1(2) = 0;
		
		p0 = cross(p1,-o);
		
	}
	else 	if(isWithinNumericalLimits(dot(o,I.normal),-1.0) ){
		p1 = Vector(3);
		p1(0) = 0;
		p1(1) = 1;
		p1(2) = 0;
		
		p0 = cross(p1,-o);

		
	}
	else{
		

		
		s = cross(v,o);
		s = I.a*s / norm(s,2);
		p0 = v+s;
		
		s2 = -cross(v,s);
		s2 = I.a*s2 / norm(s2,2);
		p1 = v+s2;
		
		if(dot(tessellationAxis,I.normal)<0){
			//p0 = v-s;
			//p1 = v-s2;
		}
		
		
	}
	
}

Interfaces Tessellation::angularInterfaces(Vector &x0, Vector &x1, TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J){
	Vector p0(3), p1(3);
	Interfaces res;
	Vector v(3);
	Vector t;
	int q;
	
	measurementPoints(p0,p1,tessellationAxis,I);
	v = I.normal * I.g;
	
	
	
	
	q = 1;
	if(I.form != CONVEX) q*=-1;
	if(J.form != CONVEX) q*=-1;
	
	
	if(q>0){
		res.out = angularInterface(x0,v,p0,p1);
		res.vectorOut = x0;
		res.in = angularInterface(x1,v,p0,p1);
		res.vectorIn = x1;
	}
	else{
		res.in = angularInterface(x0,v,p0,p1);
		res.vectorIn = x0;
		res.out = angularInterface(x1,v,p0,p1);
		res.vectorOut = x1;
	}
	
	
	
	return res;
	
}






PHIContainer Tessellation::retrieveInterfaces(TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J, float dij, float radius){
	float gij;
	float eta;
	float mu;
	float entryPoint, exitPoint;
	IntersectionPair ip;
	Vector x0(3), x1(3);
	Interfaces intf;
	Interfaces intf1;
	float baseAngleIJ;
	float rotationalPartIJ;
	float baseAngleJI;
	float rotationalPartJI;
	PHIContainer p;
	Vector Null(3);
	Null.zeros();
	
	
		ip = determineIntersectionPoints(radius, I, J);
		x0 = ip.k_j;
		x1 = ip.j_k;
		
		
		intf1 = angularInterfaces(x0,x1, tessellationAxis, I, J);
		
		

	p.in.rotation = intf1.in;
	p.in.drotation_dxi = Null;
	p.in.drotation_dxj = Null;
	p.in.drotation_dxl = Null;
	p.in.vector=intf1.vectorIn;
	p.out.rotation = intf1.out;
	p.out.drotation_dxi = Null;
	p.out.drotation_dxj = Null;
	p.out.drotation_dxl = Null;
	p.out.vector=intf1.vectorOut;
	

		
	return p;
	
	
}




float Tessellation::rotationalAngle(TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J){
	Vector ex(3);
	Vector ez(3);
	Vector n0(3);
	Vector n1(3);
	Vector n2(3);
	float sn;
	float a;
	float rho;
	
	ez(0)=0;
	ez(1)=0;
	ez(2)=1;
	float d;
	float s0,s1;

	
	d = dot(tessellationAxis,I.normal);
	
	if(isWithinNumericalLimits(d,1.0)){
		n0 = Vector(3);
		n0(0) = 0;
		n0(1) = 0;
		n0(2) = 1;
	}
	else if(isWithinNumericalLimits(d,-1.0)){
		n0 = Vector(3);
		n0(0) = 0;
		n0(1) = 0;
		n0(2) = -1;
	}
	else n0 = cross(tessellationAxis,I.normal);
	
	n2 = cross(I.normal,n0);
	
	n1 = cross(I.normal, J.normal);
	
	s0 = sgn(norm_dot(n1,n2));
	a = asin(norm_dot(n0,n1));
	s1 = sgn(a);
	
	//if(sn<0 && a>=0) a = M_PI-a;
	//else if(sn<0 && a<0) a = -M_PI-a;
	
	a = -s0*(s1*(1-s0)*M_PI/2 - a);


	
	
	return a;
	
}

*/

Matrix Tessellation::matrixCross(Matrix &m, Vector &v){
	Vector c(3);
	Vector r(3);
	Matrix res(3,3);
	int i,j;
	//produce a matrix in which each column is the crossproduct with the corresponding column in m and v
	
	for(i=0;i<3;++i){
		for(j=0;j<3;++j){
			c(j) = m(j,i);			
		}
		r = cross(c,v);
		for(j=0;j<3;++j){
			res(j,i) = r(j);			
		}
	}
	
	return res;
}

Vector Tessellation::normalise(Vector x){
	float l;
	Vector v(3);
	l = norm(x,2);
	v = x/l;
	return v;
	
}

Vector Tessellation::normalise(Vector x, float &l){
	Vector v(3);
	l = norm(x,2);
	v = x/l;
	return v;
	
}





float Tessellation::sacos(Vector &a, Vector &b){
	float la = norm(a,2);
	float lb = norm(b,2);
	return acos(dot(a,b)/(la*lb));
	
}


float Tessellation::sdot(Vector &a, Vector &b){
	float la = norm(a,2);
	float lb = norm(b,2);
	return dot(a,b)/(la*lb);
	
}

float Tessellation::l(Vector &a){
	float s = 0;
	for(int i=0; i<3; ++i)
		s += a(i)*a(i);
	return sqrt(s);
}


OmegaRotation Tessellation::calculateOmega(TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J, RhoContainer &rhoContainer){
	return calculateOmega(tessellationAxis, I.normal, J.normal, I.id, J.id, I.form, J.form, rhoContainer,I.index);

}


OmegaRotation Tessellation::calculateOmega(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, RhoContainer &rhoContainer, int index_i){
	return calculateOmega(tessellationAxis, mu_i, mu_j, 0, 0, CONVEX, CONVEX, rhoContainer, index_i);
}



OmegaRotation Tessellation::calculateOmega(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer, int index_i){
	Vector ex(3);
	Vector ez(3);
	Vector ey(3);
	Vector ni(3),ni2(3);
	Vector nij(3);
	Vector nii(3);
	float varpi;
	Matrix Identity(3,3),reverseIdentity(3,3);
	Identity.eye();
	reverseIdentity.eye();
	reverseIdentity = reverseIdentity * -1;
	OmegaRotation omega;
	float dot_ni_nij;
	bool problem=false;
	
	ez(0)=0;
	ez(1)=0;
	ez(2)=1;
	ey(0)=0;
	ey(1)=1;
	ey(2)=0;
	float dotChi_mui, dotmui_muj;
	float s0,s2;
	Vector pp(3),p(3),p2(3);
	float di;
	float dij;
	Vector vi,vij;
	
	dotChi_mui = dot(tessellationAxis.v,mu_i);
	//printf("TAX: %f %f %f\n",tessellationAxis.v(0),tessellationAxis.v(1),tessellationAxis.v(2));
	
	s2 = 1;
	if(isWithinNumericalLimits(dotChi_mui,1.0) || form_i==SPLITTER || (tessellationAxis.index==index_i && tessellationAxis.hemisphere==FRONTHEMISPHERE)){
		/*pp(0) = 0;
		pp(1) = 0;
		pp(2) = 1;
		
		p = normalise(cross(pp,mu_i));
		*/
		
		p(0) = 0;
		p(1) = 0;
		p(2) = 1;
		
		ni=normalise(cross(tessellationAxis.v,p),di);
		di=1;
		vi=cross(tessellationAxis.v,p);
		
		problem=true;
		
		
		
		if(dot(ni,mu_j)<0)
			s2 = -1;
		
	}
	else if(isWithinNumericalLimits(dotChi_mui,-1.0) || (tessellationAxis.index==index_i && tessellationAxis.hemisphere==BACKHEMISPHERE)){
		p(0) = 0;
		p(1) = -1;
		p(2) = 0;
		
		ni=normalise(cross(tessellationAxis.v,p),di);
		vi=cross(tessellationAxis.v,p);
		di=1;
		if(dot(ni,mu_j)<0)
			s2 = -1;
		
		
	}
	else{
		p = mu_i;
		ni = normalise(cross(tessellationAxis.v,p),di);
		vi=cross(tessellationAxis.v,p);
		
		
	}
		
	
	
	dotmui_muj = rhoContainer.dot_mui_muj;
	if(isWithinNumericalLimits(dotmui_muj,1.0)){
		p2(0) = 0;
		p2(1) = -1;
		p2(2) = 0;
		nij = normalise(cross(tessellationAxis.v,p2),dij);
		vij=cross(tessellationAxis.v,p2);
		dij=1;
		//exit(-2);
	}
	else if(isWithinNumericalLimits(dotmui_muj,-1.0)){
		p2(0) = 0;
		p2(1) = 1;
		p2(2) = 0;
		nij = normalise(cross(tessellationAxis.v,p2),dij);
		vij=cross(tessellationAxis.v,p2);
		dij=1;
		//exit(-3);
	}
	else{
		nij = normalise(cross(mu_i, mu_j),dij);
		vij=cross(mu_i, mu_j);
	}
	
	
	
	
	dot_ni_nij = dot(ni,nij);
	

	
	
	if(dot_ni_nij > 1.0) dot_ni_nij=1.0;
	if(dot_ni_nij < -1.0) dot_ni_nij=-1.0;
	
	varpi = acos(dot_ni_nij);
	
	s0 = -sgn(dot(nij,tessellationAxis.v));
	
	if(problem){
		s0=-1;
	}
	
	
	
	omega.rotation = s2 * s0 * varpi;
	
	
	if(isnan(omega.rotation)){
		fprintf(stdout,"omega evaluated to nan (%f %f %f %f %f %f) %d %d\n",ni(0),ni(1),ni(2),nij(0),nij(1),nij(2),tessellationAxis.index,index_i);
		printf("extra: %f\n",dotChi_mui);
		printf("tax: %f %f %f\n",tessellationAxis.v(0),tessellationAxis.v(1),tessellationAxis.v(2));
		printf("mui: %f %f %f\n",mu_i(0),mu_i(1),mu_i(2));
		exit(-2);
	}
	
	
	omega.di = di;
	omega.dij = dij;
	omega.ni = ni;
	omega.nij = nij;
	omega.vi = vi;
	omega.vij = vij;
	omega.dot_ni_nij = dot_ni_nij;
	omega.id_i = id_i;
	omega.id_j = id_j;
	omega.s0 = s0;
	omega.s2 = s2;
	
	return omega;
	
}

	
Rotation Tessellation::calculateOmegaDerivatives(OmegaRotation &r, CircularInterfacesPerAtom &circles, TessellationAxis &tessellationAxis){
	Matrix dmui_dxl, dmuj_dxl;
	float domega_dvarpi;
	Vector dvarpi_dnij, dvarpi_dni;
	Matrix dvi_dmui, dvij_dmui, dvij_dmuj;
	Matrix dni_dvi, dnij_dvij;
	Vector domega_dxi, domega_dxj, domega_dxl, domega_dxt;
	Vector domeganew_dxi, domeganew_dxj, domeganew_dxl, domeganew_dxt;
	Vector dextra_dxi, dextra_dxj, dextra_dxl, dextra_dxt;
	Matrix dvi_dchi, dchi_dxl, dchi_dxt;
	
	float di,dij;
	Vector ni(3);
	Vector nij(3);
	Matrix Identity(3,3),reverseIdentity(3,3);
	Identity.eye();
	reverseIdentity.eye();
	reverseIdentity = reverseIdentity * -1;
	Rotation omega;
	float dot_ni_nij;
	Vector mu_i;
	Vector mu_j;
	Matrix dmui_dxi;
	Matrix dmuj_dxj;
	int id_i;
	int id_j;
	float s0,s2;
	float denominator;
	Vector chi,chi1,chi2;
	Vector vi,vij;
	
	Vector dextra_dniold, dextra_dni, extra_term_dxi, extra_term_dxl, extra_term_dxt, n2;
	Matrix dniold_dvi, dviold_dmui;
	bool extra;
	Vector a0,a1,b0,b1,b2;
	float dot_ni_niold,dot_n2_nij,dot_n2_ni;
	int q0,q1;
	float a_ni_nij, a_ni_niold;
	float q2,q3,q4,q5;
	
	
	
	
	
	
	
	extra=false;
	q2=q3=q4=q5=0;
	
	omega.rotation=r.rotation;
	omega.drotation_dxi = Vector(3).zeros();
	omega.drotation_dxj = Vector(3).zeros();
	omega.drotation_dxl = Vector(3).zeros();
	omega.drotation_dxt = Vector(3).zeros();
	
	
	
	di = r.di;
	dij = r.dij;
	ni = r.ni;
	nij = r.nij;
	vi = r.vi;
	vij = r.vij;
	dot_ni_nij = r.dot_ni_nij;
	
	id_i = r.id_i;
	id_j = r.id_j;
	
	s0=r.s0;
	s2=r.s2;
	
	
	mu_i = circles[id_i].normal;
	mu_j = circles[id_j].normal;
	dmui_dxi = circles[id_i].dmu_dx;
	dmuj_dxj = circles[id_j].dmu_dx;
	
	
	
	if(abs(dot_ni_nij) >= 0.80 && circles[id_i].form!=SPLITTER && circles[id_j].form!=SPLITTER){
		//printf("TESSELLATION AXIS TRANSFORMATION %d %f %f\n", tessellationAxis.hemisphere, s0, s2);
		extra=true;
		dniold_dvi = (1.0/(di))*(Identity - kron(ni,ni.t())).t();
		
		
		chi1=Vector(3);
		chi1(0)=0;
		chi1(1)=-1;
		chi1(2)=0;
		chi2=Vector(3);
		chi2(0)=0;
		chi2(1)=0;
		chi2(2)=-1;
		
		chi = cross(chi1,tessellationAxis.v);
		ni = normalise(cross(chi,mu_i),di);
		dot_ni_nij=norm_dot(ni,nij);
		
		if(abs(dot(chi,mu_i))>0.90 || abs(dot_ni_nij)>0.80){
			chi=cross(chi2,tessellationAxis.v);
			//printf("DOUBLE CROSS\n");
		}
		
		

		ni = normalise(cross(chi,mu_i),di);
		
		
		
		vi = cross(chi,mu_i);
		dot_ni_nij=norm_dot(ni,nij);
		
		dot_ni_niold = dot(r.ni,ni);

		n2 = normalise(cross(r.ni,mu_i));
		dot_n2_nij = dot(n2,nij);
		dot_n2_ni = dot(n2,ni);
		q0 = sgn(dot_n2_nij);
		q1 = sgn(dot_n2_ni);
		
		dextra_dniold = -ni/(sqrt(1-dot(ni,r.ni)*dot(ni,r.ni)));
		dextra_dni = -r.ni/(sqrt(1-dot(ni,r.ni)*dot(ni,r.ni)));
		dviold_dmui=-matrixCross(Identity,tessellationAxis.v);
		
		a_ni_nij=acos(dot_ni_nij);
		a_ni_niold=acos(dot_ni_niold);
		
		q2=1;
		q3=1;
		q4=1;
		q5=1;
		
		//flip
		if(q0!=q1){
			a_ni_nij=M_PI-a_ni_nij;
			a_ni_niold=M_PI-a_ni_niold;
			q2=-1.0;
			q3=-1.0;
		}
		
		float A;
		if(acos(r.dot_ni_nij) > a_ni_niold)
			A=a_ni_nij+a_ni_niold;
		else{
			if(a_ni_niold>a_ni_nij){
				A=a_ni_niold-a_ni_nij;
				q4=-q4;
			}
			else{
				A=a_ni_nij-a_ni_niold;
				q5=-q5;
			}
		}
		
		a_ni_nij=acos(dot_ni_nij);
		a_ni_niold=acos(dot_ni_niold);
		
		
		A=0.5*M_PI*(q4-q2*q4+q5-q3*q5) + q2*q4*a_ni_nij + q3*q5*a_ni_niold;
		if(abs(acos(r.dot_ni_nij)- A)>0.1){
			printf("Alternative PHI derivative convergence error: %f %f\n",acos(r.dot_ni_nij),A);
			exit(-2);
		}
		
		
	}
	else chi=tessellationAxis.v;
	
	
	
	
	
	
	dmui_dxl = -dmui_dxi;
	dmuj_dxl = -dmuj_dxj;
	
	dvi_dchi = Matrix(3,3).zeros();
	dchi_dxl = Matrix(3,3).zeros();
	dchi_dxt = Matrix(3,3).zeros();
	
	
	if(s2>0){
		dvi_dmui=-matrixCross(Identity,chi);
		
		
		
	}
	else{
		dvi_dmui = Matrix(3,3).zeros();
		
	}
	//dvi_dmui=-matrixCross(Identity,chi);
	
	
	if(tessellationAxis.mode==ALIGNED){
		dvi_dchi = -matrixCross(Identity,mu_i);
		dchi_dxl = tessellationAxis.dchi_dx;
		dchi_dxt = -tessellationAxis.dchi_dx*s0; //I actually didn't look into why I need to multiply by s0, but it doesn't work without..
	}
	
	
	
	dvij_dmui = matrixCross(Identity,mu_j);
	dvij_dmuj = matrixCross(reverseIdentity,mu_i);
	
	dni_dvi = (1.0/(di))*(Identity - kron(ni,ni.t())).t();
	dnij_dvij = (1.0/(dij))*(Identity - kron(nij,nij.t())).t();
	
	
	domega_dvarpi = s0;
	
		
		
	if(extra){
		a0 = (dextra_dni.t() * dni_dvi * dvi_dmui * dmui_dxi).t();
		b0 = (dextra_dniold.t() * dniold_dvi * dviold_dmui * dmui_dxi).t();
		
		a1 = (dextra_dni.t() * dni_dvi * dvi_dmui * dmui_dxl).t();
		b1 = (dextra_dniold.t() * dniold_dvi * (dviold_dmui * dmui_dxl + dvi_dchi * dchi_dxl)).t();

		//a2 = (dextra_dni.t() * dni_dvi * dvi_dmui * dmui_dxl).t();
		b2 = (dextra_dniold.t() * dniold_dvi * (dvi_dchi * dchi_dxt)).t();

		
		extra_term_dxi =  (a0+b0);
		extra_term_dxl =  (a1+b1);
		extra_term_dxt = b2;
	}
		
	if(dot_ni_nij >= 1.0) dot_ni_nij = 1.0-MINISCULE;
	if(dot_ni_nij <= -1.0) dot_ni_nij = -1.0+MINISCULE;
	
	

	denominator = (sqrt(1-dot_ni_nij*dot_ni_nij));
	
	dvarpi_dni = -(nij) / denominator;
	dvarpi_dnij = -(ni) / denominator;
	

	domega_dxi = s2* domega_dvarpi * ((dvarpi_dni.t() * dni_dvi * dvi_dmui * dmui_dxi).t() + (dvarpi_dnij.t() * dnij_dvij * dvij_dmui * dmui_dxi).t() );
	domega_dxj = s2* domega_dvarpi * (dvarpi_dnij.t() * dnij_dvij * dvij_dmuj * dmuj_dxj).t();
	domega_dxl = s2* domega_dvarpi * ( (dvarpi_dni.t() * dni_dvi * (dvi_dmui * dmui_dxl + dvi_dchi * dchi_dxl)).t() +  (dvarpi_dnij.t() * (dnij_dvij * dvij_dmui * dmui_dxl + dnij_dvij * dvij_dmuj * dmuj_dxl) ).t() );
	domega_dxt = s2* domega_dvarpi * ( (dvarpi_dni.t() * dni_dvi * (dvi_dchi * dchi_dxt)).t());

	
	
	
	if(extra){
		domega_dxi =  ((dvarpi_dni.t() * dni_dvi * dvi_dmui * dmui_dxi).t() + (dvarpi_dnij.t() * dnij_dvij * dvij_dmui * dmui_dxi).t() );
		domega_dxj =  (dvarpi_dnij.t() * dnij_dvij * dvij_dmuj * dmuj_dxj).t();
		//domega_dxl =  ( (dvarpi_dni.t() * dni_dvi * (dvi_dmui * dmui_dxl + dvi_dchi * dchi_dxl)).t() +  (dvarpi_dnij.t() * (dnij_dvij * dvij_dmui * dmui_dxl + dnij_dvij * dvij_dmuj * dmuj_dxl) ).t() );
		domega_dxl =  ( (dvarpi_dni.t() * dni_dvi * (dvi_dmui * dmui_dxl)).t() +  (dvarpi_dnij.t() * (dnij_dvij * dvij_dmui * dmui_dxl + dnij_dvij * dvij_dmuj * dmuj_dxl) ).t() );
		//domega_dxt =  ( (dvarpi_dni.t() * dni_dvi * (dvi_dchi * dchi_dxt)).t());
		domega_dxt = Vector(3).zeros();
		
		domeganew_dxi = q2*q4*domega_dxi;
		domeganew_dxj = q2*q4*domega_dxj;
		domeganew_dxl = q2*q4*domega_dxl;
		domeganew_dxt = q2*q4*domega_dxt;
		
		
		
		domega_dxi= s2*domega_dvarpi*(q2*q4*domega_dxi + q3*q5*extra_term_dxi);
		domega_dxl= s2*domega_dvarpi*(q2*q4*domega_dxl + q3*q5*extra_term_dxl);
		domega_dxj= s2*domega_dvarpi*(q2*q4*domega_dxj);
		domega_dxt= s2*domega_dvarpi*(q2*q4*domega_dxt) + q3*q5*extra_term_dxt;
		
		dextra_dxi = q3*q5*extra_term_dxi;
		dextra_dxj = Vector(3).zeros();
		dextra_dxl = q3*q5*extra_term_dxl;
	}
	
	omega.drotation_dxi = domega_dxi;
	omega.drotation_dxj = domega_dxj;
	omega.drotation_dxl = domega_dxl;
	omega.drotation_dxt = domega_dxt;
	
	
	
	
	
	
	
	
	


	
	return omega;
	
}




float Tessellation::calculateRho(Vector &mu_i, Vector &mu_j){
	return calculateRho(mu_i, mu_j, false);
}


float Tessellation::calculateRho(Vector &mu_i, Vector &mu_j, bool derivatives){
	float dot_IJ;
	float rho;
	dot_IJ = norm_dot(mu_i,mu_j);
	rho = acos(dot_IJ);
	return rho;
}


EtaRotation Tessellation::calculateEta(TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J, RhoContainer &rhoContainer){
	return calculateEta(tessellationAxis, I.normal, J.normal, I.lambda, J.lambda, I.id, J.id, I.form, J.form, rhoContainer);
}


EtaRotation Tessellation::calculateEta(Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, RhoContainer &rhoContainer){
	TessellationAxis tessellationAxis;
	tessellationAxis.v=Vector(3);
	tessellationAxis.v.zeros();
	return calculateEta(tessellationAxis, mu_i, mu_j, lambda_i, lambda_j, 0, 0, CONVEX, CONVEX, rhoContainer);
}


EtaRotation Tessellation::calculateEta(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer){
	float rho;
	float dot_IJ;
	float sig0,sig1,sig2;
	EtaRotation eta;
	float csc_lambda_i;
	float cos_lambda_i;
	float cot_lambda_i;
	float csc_rho;
	float cos_rho;
	float cot_rho;
	
	dot_IJ = norm_dot(mu_i,mu_j);
	//dot_IJ = rhoContainer.dot_mui_muj;
	
	rho = rhoContainer.rho;
	
	
	
	if(isWithinNumericalLimits(rho,0)){
		eta.rotation = 0;
	}
	else{
		csc_lambda_i = 1.0/sin(lambda_i.rotation);
		cos_lambda_i = cos(lambda_i.rotation);
		cot_lambda_i = cos_lambda_i * csc_lambda_i;

		csc_rho = 1.0/sin(rho);
		cos_rho = cos(rho);
		cot_rho = cos_rho * csc_rho;
		
		
		//sig0 = cot(lambda_i.rotation) * cot(rho);
		//sig1 = cos(lambda_j.rotation)*csc(lambda_i.rotation)*csc(rho);
		sig0 = cot_lambda_i * cot_rho;
		sig1 = cos(lambda_j.rotation)*csc_lambda_i*csc_rho;
		sig2 = sig0-sig1;
		
// 		printf("sigs: %f %f %f %f %f %f\n",sig0,sig1,sig2, lambda_i.rotation, lambda_j.rotation, rho);
		
		
		eta.rotation = acos(sig2);
		
		eta.sig0 = sig0;
		eta.sig1 = sig1;
		eta.sig2 = sig2;
		eta.dot_IJ = dot_IJ;
		eta.rho = rho;

		
		
		
	}
	
	if(isnan(eta.rotation)){
		fprintf(stdout,"eta evaluated to nan\n");
		fprintf(stdout,"sig: %f %f %f\n",sig0,sig1,sig2);
		fprintf(stdout,"angles: %f %f %f\n",lambda_i.rotation, lambda_j.rotation, rho);
		//exit(-2);
	}
	
	
	eta.id_i = id_i;
	eta.id_j = id_j;
	
	
	return eta;
}
	
	
	
Rotation Tessellation::calculateEtaDerivatives(EtaRotation &r, CircularInterfacesPerAtom &circles){
					       
	Vector drho_dmui(3), drho_dmuj(3);
	Vector drho_dxi(3), drho_dxj(3);
	Vector drho_dxl(3);
	float sig0,sig1,sig2;
	float deta_dlambdai, deta_dlambdaj, deta_drhoij;
	Vector deta_dxi(3), deta_dxj(3), deta_dxl(3), deta_dxt(3);
	Rotation eta;
	Matrix dmui_dxl, dmuj_dxl;
	Vector mu_i;
	Vector mu_j;
	Rotation lambda_i;
	Rotation lambda_j;
	Matrix dmui_dxi;
	Matrix dmuj_dxj;
	int id_i;
	int id_j;
	float dot_IJ;
	float rho;
	
	sig0 = r.sig0;
	sig1 = r.sig1;
	sig2 = r.sig2;
	dot_IJ = r.dot_IJ;
	rho = r.rho;
	id_i = r.id_i;
	id_j = r.id_j;
	
	mu_i = circles[id_i].normal;
	mu_j = circles[id_j].normal;
	lambda_i = circles[id_i].lambda;
	lambda_j = circles[id_j].lambda;
	dmui_dxi = circles[id_i].dmu_dx;
	dmuj_dxj = circles[id_j].dmu_dx;
	

	eta.rotation = r.rotation;
	
	
	eta.drotation_dxi = Vector(3).zeros();
	eta.drotation_dxj = Vector(3).zeros();
	eta.drotation_dxl = Vector(3).zeros();
	eta.drotation_dxt = Vector(3).zeros();
	

	if(isWithinNumericalLimits(rho,0)){
		deta_dxi = Vector(3).zeros();
		deta_dxj = Vector(3).zeros();
		deta_dxl = Vector(3).zeros();
		deta_dxt = Vector(3).zeros();
	}
	else{
		
		if(dot_IJ > 1.0-MINISCULE) dot_IJ = 1.0-MINISCULE;
		if(dot_IJ < -1.0+MINISCULE) dot_IJ = -1.0+MINISCULE;
		
	
		drho_dmui = -mu_j/(sqrt(1-dot_IJ*dot_IJ));
		drho_dmuj = -mu_i/(sqrt(1-dot_IJ*dot_IJ));
		
		dmui_dxl = -dmui_dxi;
		dmuj_dxl = -dmuj_dxj;
		
		
		drho_dxi = (drho_dmui.t() * dmui_dxi).t();
		drho_dxj = (drho_dmuj.t() * dmuj_dxj).t();
		
		drho_dxl = ((drho_dmui.t() * dmui_dxl).t() + (drho_dmuj.t() * dmuj_dxl).t());
		
		if(sig2 > 1.0-MINISCULE) sig2 = 1.0-MINISCULE;
		if(sig2 < -1.0+MINISCULE) sig2 = -1.0+MINISCULE;
		
// 		printf("rotation: %f %f %f %f %f\n",lambda_i.rotation, lambda_j.rotation, sig0, sig1, sig2);
		
		deta_dlambdai = csc(lambda_i.rotation) * (sig0/cos(lambda_i.rotation) - sig1*cos(lambda_i.rotation)) / sqrt(1-sig2*sig2);
		deta_dlambdaj =  -sin(lambda_j.rotation)*csc(lambda_i.rotation)*csc(rho) / sqrt(1-sig2*sig2);
		deta_drhoij = csc(rho) * (sig0/cos(rho) - sig1*cos(rho)) / sqrt(1-sig2*sig2);
		
		
		
		deta_dxi = deta_dlambdai * lambda_i.drotation_dxi + deta_drhoij * drho_dxi;
		deta_dxj = deta_dlambdaj * lambda_j.drotation_dxi + deta_drhoij * drho_dxj;
		deta_dxl = deta_dlambdai * lambda_i.drotation_dxl + deta_dlambdaj * lambda_j.drotation_dxl + deta_drhoij * drho_dxl;
		
		
		eta.drotation_dxi = deta_dxi;
		eta.drotation_dxj = deta_dxj;
		eta.drotation_dxl = deta_dxl;
		eta.drotation_dxt = Vector(3).zeros();
	}
	
	
	
	return eta;
}

PHIContainer Tessellation::calculatePHI(TessellationAxis &tessellationAxis, CircularInterface &I, CircularInterface &J, float radius, RhoContainer &rhoContainer){
	return calculatePHI(tessellationAxis, I.normal, J.normal, I.lambda, J.lambda, I.id, J.id, I.form, J.form, rhoContainer, I.index);
	
}


PHIContainer Tessellation::calculatePHI(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, RhoContainer &rhoContainer, int index_i){
	return calculatePHI(tessellationAxis, mu_i, mu_j, lambda_i, lambda_j, 0, 0, CONVEX, CONVEX, rhoContainer, index_i);
}


PHIContainer Tessellation::calculatePHI(TessellationAxis &tessellationAxis, Vector &mu_i, Vector &mu_j, Rotation &lambda_i, Rotation &lambda_j, int id_i, int id_j, CircularInterfaceForm form_i, CircularInterfaceForm form_j, RhoContainer &rhoContainer, int index_i){
	PHIContainer p,p2,p3;
	
	EtaRotation eta;
	OmegaRotation omega;
	int q;
	

	eta = calculateEta(tessellationAxis, mu_i, mu_j, lambda_i, lambda_j, id_i, id_j, form_i, form_j, rhoContainer);
	omega = calculateOmega(tessellationAxis, mu_i, mu_j, id_i, id_j, form_i, form_j, rhoContainer, index_i);
	
	
	
	p.out.rotation = omega.rotation + eta.rotation;
	if(p.out.rotation > M_PI)
		p.out.rotation = -2*M_PI + p.out.rotation;
	
	p.in.rotation = omega.rotation - eta.rotation;
	if(p.in.rotation < -M_PI)
		p.in.rotation = 2*M_PI + p.in.rotation;

	
	
	p.in.omega = omega;
	p.in.eta = eta;
	p.out.omega = omega;
	p.out.eta = eta;

	p.out.s = 1;
	p.in.s = -1;
	
	
	
	q = 1;
	if(form_i != CONVEX) q*=-1;
	if(form_j != CONVEX) q*=-1;
	
	
	
	
	if(q<0){
		p3=p;
		p.in=p3.out;
		p.out=p3.in;
	}

	
	


	
	return p;
	
	
	 
}


Rotation Tessellation::calculatePHIDerivatives(PHIRotation &r, CircularInterfacesPerAtom &circles, TessellationAxis &tessellationAxis){
	Rotation omega,eta;
	Rotation PHI;
	int s;
	
	
	
	omega = calculateOmegaDerivatives(r.omega, circles, tessellationAxis);
	eta = calculateEtaDerivatives(r.eta, circles);
	
	PHI.rotation = r.rotation;
	s=r.s;
		
	PHI.drotation_dxi = s*eta.drotation_dxi + omega.drotation_dxi;
	PHI.drotation_dxj = s*eta.drotation_dxj + omega.drotation_dxj;
	PHI.drotation_dxl = s*eta.drotation_dxl + omega.drotation_dxl;
	PHI.drotation_dxt = s*eta.drotation_dxt + omega.drotation_dxt;
	
	
	
	
	/*
	
	printf("=====================================================================\n");
	printf("TAX index: %d\n",tessellationAxis.index);
	//printf("TAX: %f %f %f\n",tessellationAxis(0),tessellationAxis(1),tessellationAxis(2));
	
	
	
	
	
		float d, g_chi, r_l, r_chi;
		Vector ai(3), aj(3), al(3), x(3), xp(3), xn(3), achi(3), chi(3);
		Vector nip,nin,nj;
		Vector njp,njn,ni;
		int indexi, indexl, indexj,idi,idj;
		float ri, rl, rj;
		Rotation lambdaip;
		Rotation lambdain;
		Rotation lambdaj;
		float dip,din,dj;
		Vector fd_dPHIij_dxi_in(3);
		Vector fd_dPHIij_dxi_out(3);
		Vector fd_detaij_dxi_in(3);
		Vector fd_detaij_dxi_out(3);
		Vector fd_domegaij_dxi_in(3);
		Vector fd_domegaij_dxi_out(3);
		Vector fd_dPHIij_dxj_in(3);
		Vector fd_dPHIij_dxj_out(3);
		Vector fd_detaij_dxj_in(3);
		Vector fd_detaij_dxj_out(3);
		Vector fd_domegaij_dxj_in(3);
		Vector fd_domegaij_dxj_out(3);
		Vector fd_dPHIij_dxl_in(3);
		Vector fd_dPHIij_dxl_out(3);
		Vector fd_detaij_dxl_in(3);
		Vector fd_detaij_dxl_out(3);
		Vector fd_domegaij_dxl_in(3);
		Vector fd_domegaij_dxl_out(3);
		Vector fd_dPHIij_dxt_in(3);
		Vector fd_dPHIij_dxt_out(3);
		Vector fd_detaij_dxt_in(3);
		Vector fd_detaij_dxt_out(3);
		Vector fd_domegaij_dxt_in(3);
		Vector fd_domegaij_dxt_out(3);
		PHIContainer PHIp, PHIn;
		Rotation lambdajp;
		Rotation lambdajn;
		Rotation lambdai;
		float djp,djn,di;
		RhoContainer rhop,rhon;
		float etap,etan,omegap,omegan;
		TessellationAxis taxp;	
		TessellationAxis taxn,tax;
		Vector chip,chin;

		r_l = tradius;
		r_chi = radii[tessellationAxis.index];

		
		//dphi_dxi
		indexi = circles[r.eta.id_i].index;
		indexj = circles[r.eta.id_j].index;
		
		printf("INDEX: %d %d %d\n",ti,indexi,indexj);
			printf("PC_dxi: (%f %f %f)\n",PHI.drotation_dxi(0),PHI.drotation_dxi(1),PHI.drotation_dxi(2));
			printf("PC_dxj: (%f %f %f)\n",PHI.drotation_dxj(0),PHI.drotation_dxj(1),PHI.drotation_dxj(2));
			printf("PC_dxl: (%f %f %f)\n",PHI.drotation_dxl(0),PHI.drotation_dxl(1),PHI.drotation_dxl(2));
		
		
		idi=r.eta.id_i;
		idj=r.eta.id_j;
		
			al = torigin;
			rl = tradius;
			
			if(tessellationAxis.mode==ALIGNED){
				achi = atoms[tessellationAxis.index];
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE){
					chi=calculateInterfaceNormal(al, achi, d);
				}
				else{
					chi=-calculateInterfaceNormal(al, achi, d);
				}
				
				g_chi = (d * d + r_l * r_l - r_chi * r_chi ) / (2 * d * r_l);

				if(g_chi<0){
					chi=-chi;
				}
			}
			else{
				chi=tessellationAxis.v;
			}
			
			
			if(indexi<0){
				ni = chi;
				lambdai.rotation = M_PI/2.0;

			}
			else{
				ai = atoms[indexi];
				ri = radii[indexi];
				
				ni = calculateInterfaceNormal(al, ai, di);
				lambdai.rotation = calculateLambda(di, rl, ri, ni).rotation;

			}			
			
			if(indexj<0){
				nj = chi;
				lambdaj.rotation = M_PI/2.0;
				
			}
			else{
				aj = atoms[indexj];
				rj = radii[indexj];
				nj = calculateInterfaceNormal(al, aj, dj);
				
				lambdaj.rotation = calculateLambda(dj, rl, rj, nj).rotation;
			}
			
			
						
			rhop.rho = calculateRho(ni,nj);
			tax=tessellationAxis;
			tax.mode=ALIGNED;
			taxp.v=chi;
			taxp=tessellationAxis;
			taxp.mode=ALIGNED;
			taxn.mode=ALIGNED;
			PHIp = calculatePHI(taxp, ni, nj, lambdai, lambdaj, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi);
			etap = calculateEta(taxp, ni, nj, lambdai, lambdaj, idi, idj, circles[idi].form, circles[idj].form, rhop).rotation;
			omegap = calculateOmega(taxp, ni, nj, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi).rotation;

			
			printf("PRECOMPARE: %f / %f %f\n",PHI.rotation, PHIp.in.rotation, PHIp.out.rotation);
			printf("PRECOMPARE e: %f / %f\n",eta.rotation, etap);
			printf("PRECOMPARE o: %f / %f\n",omega.rotation, omegap);
			
			
		if(indexi>=0){
			if(indexi<0){
				ai=atoms[tessellationAxis.index];
			}
			else{
				ai=atoms[indexi];
			}
			
			x=ai;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				
				if(indexi >=0 ){
					ri = radii[indexi];
					
					if(indexj<0){
						nip = calculateInterfaceNormal(al, xp, dip);
						nin = calculateInterfaceNormal(al, xn, din);
						nj = chi;
						
						lambdaip.rotation = calculateLambda(dip, rl, ri, nip).rotation;
						lambdain.rotation = calculateLambda(din, rl, ri, nin).rotation;
						lambdaj.rotation = M_PI/2.0;
						
					}
					else{
						aj = atoms[indexj];
						rj = radii[indexj];
						
						nip = calculateInterfaceNormal(al, xp, dip);
						nin = calculateInterfaceNormal(al, xn, din);
						nj = calculateInterfaceNormal(al, aj, dj);
						
						lambdaip.rotation = calculateLambda(dip, rl, ri, nip).rotation;
						lambdain.rotation = calculateLambda(din, rl, ri, nin).rotation;
						lambdaj.rotation = calculateLambda(dj, rl, rj, nj).rotation;
					}
					
					chip=chi;
					chin=chi;
				}
				else{
				
					if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
						chip=calculateInterfaceNormal(al, xp, d);
					else
						chip=-calculateInterfaceNormal(al, xp, d);
					
					if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
						chin=calculateInterfaceNormal(al, xn, d);
					else
						chin=-calculateInterfaceNormal(al, xn, d);
					
					
					if(g_chi<0){
						chip=-chip;
						chin=-chin;
					}						
					
					aj = atoms[indexj];
					rj = radii[indexj];
					
					
					nip = chip;
					nin = chin;
					nj = calculateInterfaceNormal(al, aj, dj);
					
					lambdaip.rotation = M_PI/2.0;
					lambdain.rotation = M_PI/2.0;
					lambdaj.rotation = calculateLambda(dj, rl, rj, nj).rotation;

						
				}
				
				
				
				rhop.rho = calculateRho(nip,nj);
				rhon.rho = calculateRho(nin,nj);
				
				taxp.v=chip;
				taxn.v=chin;
				
				
				
				PHIp = calculatePHI(taxp, nip, nj, lambdaip, lambdaj, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi);
				PHIn = calculatePHI(taxn, nin, nj, lambdain, lambdaj, idi, idj, circles[idi].form, circles[idj].form, rhon, indexi);
				etap = calculateEta(taxp, nip, nj, lambdaip, lambdaj, idi, idj, circles[idi].form, circles[idj].form, rhop).rotation;
				etan = calculateEta(taxn, nin, nj, lambdain, lambdaj, idi, idj, circles[idi].form, circles[idj].form, rhon).rotation;
				omegap = calculateOmega(taxp, nip, nj, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi).rotation;
				omegan = calculateOmega(taxn, nin, nj, idi, idj, circles[idi].form, circles[idj].form, rhon, indexi).rotation;
				
				fd_dPHIij_dxi_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
				fd_dPHIij_dxi_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
				fd_detaij_dxi_out(j) = (etap - etan) / (2*FD);
				fd_detaij_dxi_in(j) = (etap - etan) / (2*FD);
				fd_domegaij_dxi_out(j) = (omegap - omegan) / (2*FD);
				fd_domegaij_dxi_in(j) = (omegap - omegan) / (2*FD);
			}
			printf("RAW: PHI: %f %f / %f %f OMEGA: %f %f\n",PHIp.out.rotation, PHIn.out.rotation,PHIp.in.rotation, PHIn.in.rotation, omegap, omegan);
			printf("dPHIij_dxi: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",PHI.drotation_dxi(0),PHI.drotation_dxi(1),PHI.drotation_dxi(2), fd_dPHIij_dxi_in(0), fd_dPHIij_dxi_in(1), fd_dPHIij_dxi_in(2), fd_dPHIij_dxi_out(0), fd_dPHIij_dxi_out(1), fd_dPHIij_dxi_out(2));
			printf("detaij_dxi: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",eta.drotation_dxi(0),eta.drotation_dxi(1),eta.drotation_dxi(2), fd_detaij_dxi_in(0), fd_detaij_dxi_in(1), fd_detaij_dxi_in(2), fd_detaij_dxi_out(0), fd_detaij_dxi_out(1), fd_detaij_dxi_out(2));
			printf("domegaij_dxi: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",omega.drotation_dxi(0),omega.drotation_dxi(1),omega.drotation_dxi(2), fd_domegaij_dxi_in(0), fd_domegaij_dxi_in(1), fd_domegaij_dxi_in(2), fd_domegaij_dxi_out(0), fd_domegaij_dxi_out(1), fd_domegaij_dxi_out(2));
			for(int i=0; i<3; ++i){
 				if(abs(PHI.drotation_dxi(i) - fd_dPHIij_dxi_in(i)) > FDT && abs(PHI.drotation_dxi(i) - fd_dPHIij_dxi_out(i)) > FDT) printf("ERROR\n");//printf("DIVERGENCE\n");
				if(abs(omega.drotation_dxi(i) - fd_domegaij_dxi_in(i)) > FDT && abs(omega.drotation_dxi(i) - fd_domegaij_dxi_out(i)) > FDT) printf("ERROR\n");//printf("DIVERGENCE\n");
			}
			
		
		}
		
		
		if(indexj>=0){
		
		//dphi_dxj
		if(indexj<0){
			aj=atoms[tessellationAxis.index];
		}
		else{
			aj=atoms[indexj];
		}
		al = torigin;
		rl = tradius;
		
		
		
		
		x=aj;
		for(int j=0; j<3; ++j){
			xp=x;
			xp(j) = x(j) + FD;
			xn=x;
			xn(j) = x(j) - FD;
			
				
		if(indexj>=0){
			rj = radii[indexj];
			
			if(indexi<0){
				njp = calculateInterfaceNormal(al, xp, djp);
				njn = calculateInterfaceNormal(al, xn, djn);
				ni = chi;
				
				

				
				lambdajp.rotation = calculateLambda(djp, rl, rj, njp).rotation;
				lambdajn.rotation = calculateLambda(djn, rl, rj, njn).rotation;
				lambdai.rotation = M_PI/2.0;

			}
			else{
				ai = atoms[indexi];
				ri = radii[indexi];
				
				njp = calculateInterfaceNormal(al, xp, djp);
				njn = calculateInterfaceNormal(al, xn, djn);
				ni = calculateInterfaceNormal(al, ai, di);
				
				

				
				lambdajp.rotation = calculateLambda(djp, rl, rj, njp).rotation;
				lambdajn.rotation = calculateLambda(djn, rl, rj, njn).rotation;
				lambdai.rotation = calculateLambda(di, rl, ri, ni).rotation;

			}
			chip=chi;
			chin=chi;
		}
		else{
			
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
					chip=calculateInterfaceNormal(al, xp, d);
				else
					chip=-calculateInterfaceNormal(al, xp, d);
				
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
					chin=calculateInterfaceNormal(al, xn, d);
				else
					chin=-calculateInterfaceNormal(al, xn, d);
				
				
				if(g_chi<0){
					chip=-chip;
					chin=-chin;
				}						
								
				ai = atoms[indexi];
				ri = radii[indexi];
				
				njp = chip;
				njn = chin;
				ni = calculateInterfaceNormal(al, ai, di);
				
				

				
				lambdajp.rotation = M_PI/2.0;
				lambdajn.rotation = M_PI/2.0;
				lambdai.rotation = calculateLambda(di, rl, ri, ni).rotation;
			
		}
						
			taxp.v=chip;
			taxn.v=chin;
				
				
				rhop.rho = calculateRho(ni,njp);
				rhon.rho = calculateRho(ni,njn);
				
				
				PHIp = calculatePHI(taxp, ni, njp, lambdai, lambdajp, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi);
				PHIn = calculatePHI(taxn, ni, njn, lambdai, lambdajn, idi, idj, circles[idi].form, circles[idj].form, rhon, indexi);
				etap = calculateEta(taxp, ni, njp, lambdai, lambdajp, idi, idj, circles[idi].form, circles[idj].form, rhop).rotation;
				etan = calculateEta(taxn, ni, njn, lambdai, lambdajn, idi, idj, circles[idi].form, circles[idj].form, rhon).rotation;
				omegap = calculateOmega(taxp, ni, njp, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi).rotation;
				omegan = calculateOmega(taxn, ni, njn, idi, idj, circles[idi].form, circles[idj].form, rhon, indexi).rotation;
				
				fd_dPHIij_dxj_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
				fd_dPHIij_dxj_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
				fd_detaij_dxj_out(j) = (etap - etan) / (2*FD);
				fd_detaij_dxj_in(j) = (etap - etan) / (2*FD);
				fd_domegaij_dxj_out(j) = (omegap - omegan) / (2*FD);
				fd_domegaij_dxj_in(j) = (omegap - omegan) / (2*FD);
			}
			printf("dPHIij_dxj: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",PHI.drotation_dxj(0),PHI.drotation_dxj(1),PHI.drotation_dxj(2), fd_dPHIij_dxj_in(0), fd_dPHIij_dxj_in(1), fd_dPHIij_dxj_in(2), fd_dPHIij_dxj_out(0), fd_dPHIij_dxj_out(1), fd_dPHIij_dxj_out(2));
			printf("detaij_dxj: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",eta.drotation_dxj(0),eta.drotation_dxj(1),eta.drotation_dxj(2), fd_detaij_dxj_in(0), fd_detaij_dxj_in(1), fd_detaij_dxj_in(2), fd_detaij_dxj_out(0), fd_detaij_dxj_out(1), fd_detaij_dxj_out(2));
			printf("domegaij_dxj: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",omega.drotation_dxj(0),omega.drotation_dxj(1),omega.drotation_dxj(2), fd_domegaij_dxj_in(0), fd_domegaij_dxj_in(1), fd_domegaij_dxj_in(2), fd_domegaij_dxj_out(0), fd_domegaij_dxj_out(1), fd_domegaij_dxj_out(2));
			for(int i=0; i<3; ++i){
 				if(abs(PHI.drotation_dxj(i) - fd_dPHIij_dxj_in(i)) > FDT && abs(PHI.drotation_dxj(i) - fd_dPHIij_dxj_out(i)) > FDT) printf("ERROR\n");//printf("DIVERGENCE\n");
				if(abs(omega.drotation_dxj(i) - fd_domegaij_dxj_in(i)) > FDT && abs(omega.drotation_dxj(i) - fd_domegaij_dxj_out(i)) > FDT) printf("ERROR\n");//printf("DIVERGENCE\n");
			}
		
			
		}
		
		//dphi_dxl
		if(true){
			al = torigin;
			rl = tradius;
			
			
			x=al;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				
				if(tessellationAxis.mode==ALIGNED){

					if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
						chip=calculateInterfaceNormal(xp, achi, d);
					else
						chip=-calculateInterfaceNormal(xp, achi, d);
					
					if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
						chin=calculateInterfaceNormal(xn, achi, d);
					else
						chin=-calculateInterfaceNormal(xn, achi, d);
					
					if(g_chi<0){
						chip=-chip;
						chin=-chin;
					}				
				}
				else{
					chip=tessellationAxis.v;
					chin=tessellationAxis.v;
				}
				
				
				
			if(indexi<0){
				nip=chip;
				nin=chin;
				lambdaip.rotation = M_PI/2.0;
				lambdain.rotation = M_PI/2.0;
				
			}
			else{
				ai = atoms[indexi];
				ri = radii[indexi];
				nip = calculateInterfaceNormal(xp, ai, dip);
				nin = calculateInterfaceNormal(xn, ai, din);
				lambdaip.rotation = calculateLambda(dip, rl, ri, nip).rotation;
				lambdain.rotation = calculateLambda(din, rl, ri, nin).rotation;
			}
				
				
			if(indexj<0){
				njp=chip;
				njn=chin;
				lambdajp.rotation = M_PI/2.0;
				lambdajn.rotation = M_PI/2.0;
			}
			else{
				aj = atoms[indexj];
				rj = radii[indexj];
				njp = calculateInterfaceNormal(xp, aj, djp);
				njn = calculateInterfaceNormal(xn, aj, djn);
				lambdajp.rotation = calculateLambda(djp, rl, rj, njp).rotation;
				lambdajn.rotation = calculateLambda(djn, rl, rj, njn).rotation;
			}
				
				
			taxp.v=chip;
			taxn.v=chin;
				
// 				printf("XP: %f %f %f , %f %f %f\n",xp(0),xp(1),xp(2),xn(0),xn(1),xn(2));
				
				
				
				
// 				printf("RHO: %f %f\n",rhop.rho, rhon.rho);
				
				
				
				rhop.rho = calculateRho(nip,njp);
				rhon.rho = calculateRho(nin,njn);
				
				
// 				printf("LAMBDA: %f %f %f %f\n",lambdaip.rotation,lambdain.rotation,lambdajp.rotation,lambdajp.rotation);
				
				
				PHIp = calculatePHI(taxp, nip, njp, lambdaip, lambdajp, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi);
				PHIn = calculatePHI(taxn, nin, njn, lambdain, lambdajn, idi, idj, circles[idi].form, circles[idj].form, rhon, indexi);
				etap = calculateEta(taxp, nip, njp, lambdaip, lambdajp, idi, idj, circles[idi].form, circles[idj].form, rhop).rotation;
				etan = calculateEta(taxn, nin, njn, lambdain, lambdajn, idi, idj, circles[idi].form, circles[idj].form, rhon).rotation;
				omegap = calculateOmega(taxp, nip, njp, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi).rotation;
				omegan = calculateOmega(taxn, nin, njn, idi, idj, circles[idi].form, circles[idj].form, rhon, indexi).rotation;
				
// 				printf("eta: %f %f\n",etap,etan);
// 				printf("omega: %f %f\n",omegap,omegan);
				
				
				fd_dPHIij_dxl_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
				fd_dPHIij_dxl_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
				fd_detaij_dxl_out(j) = (etap - etan) / (2*FD);
				fd_detaij_dxl_in(j) = (etap - etan) / (2*FD);
				fd_domegaij_dxl_out(j) = (omegap - omegan) / (2*FD);
				fd_domegaij_dxl_in(j) = (omegap - omegan) / (2*FD);
			}
			printf("dPHIij_dxl: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",PHI.drotation_dxl(0),PHI.drotation_dxl(1),PHI.drotation_dxl(2), fd_dPHIij_dxl_in(0), fd_dPHIij_dxl_in(1), fd_dPHIij_dxl_in(2), fd_dPHIij_dxl_out(0), fd_dPHIij_dxl_out(1), fd_dPHIij_dxl_out(2));
			printf("detaij_dxl: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",eta.drotation_dxl(0),eta.drotation_dxl(1),eta.drotation_dxl(2), fd_detaij_dxl_in(0), fd_detaij_dxl_in(1), fd_detaij_dxl_in(2), fd_detaij_dxl_out(0), fd_detaij_dxl_out(1), fd_detaij_dxl_out(2));
			printf("domegaij_dxl: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",omega.drotation_dxl(0),omega.drotation_dxl(1),omega.drotation_dxl(2), fd_domegaij_dxl_in(0), fd_domegaij_dxl_in(1), fd_domegaij_dxl_in(2), fd_domegaij_dxl_out(0), fd_domegaij_dxl_out(1), fd_domegaij_dxl_out(2));
			for(int i=0; i<3; ++i){
 				if(abs(PHI.drotation_dxl(i) - fd_dPHIij_dxl_in(i)) > FDT && abs(PHI.drotation_dxl(i) - fd_dPHIij_dxl_out(i)) > FDT) printf("ERROR\n");//printf("DIVERGENCE\n");
				if(abs(omega.drotation_dxl(i) - fd_domegaij_dxl_in(i)) > FDT && abs(omega.drotation_dxl(i) - fd_domegaij_dxl_out(i)) > FDT) printf("ERROR\n");//printf("DIVERGENCE\n");
			}
		}
		
		
		
		
		
		if(tessellationAxis.mode==ALIGNED){
		//dphi_dxt
		if(true){
			al = torigin;
			rl = tradius;
			
			
			x=achi;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
					chip=calculateInterfaceNormal(al, xp, d);
				else
					chip=-calculateInterfaceNormal(al, xp, d);
				
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
					chin=calculateInterfaceNormal(al, xn, d);
				else
					chin=-calculateInterfaceNormal(al, xn, d);
				
				
				if(g_chi<0){
					chip=-chip;
					chin=-chin;
				}	
				
				
				
			if(indexi<0){
				nip=chip;
				nin=chin;
				lambdaip.rotation = M_PI/2.0;
				lambdain.rotation = M_PI/2.0;
				
			}
			else{
				ai = atoms[indexi];
				ri = radii[indexi];
				nip = calculateInterfaceNormal(al, ai, dip);
				nin = calculateInterfaceNormal(al, ai, din);
				lambdaip.rotation = calculateLambda(dip, rl, ri, nip).rotation;
				lambdain.rotation = calculateLambda(din, rl, ri, nin).rotation;
			}
				
				
			if(indexj<0){
				njp=chip;
				njn=chin;
				lambdajp.rotation = M_PI/2.0;
				lambdajn.rotation = M_PI/2.0;
			}
			else{
				aj = atoms[indexj];
				rj = radii[indexj];
				njp = calculateInterfaceNormal(al, aj, djp);
				njn = calculateInterfaceNormal(al, aj, djn);
				lambdajp.rotation = calculateLambda(djp, rl, rj, njp).rotation;
				lambdajn.rotation = calculateLambda(djn, rl, rj, njn).rotation;
			}
				
			
			taxp.v=chip;
			taxn.v=chin;
			
				
// 				printf("XP: %f %f %f , %f %f %f\n",xp(0),xp(1),xp(2),xn(0),xn(1),xn(2));
				
				
				
				
// 				printf("RHO: %f %f\n",rhop.rho, rhon.rho);
				
				
				
				rhop.rho = calculateRho(nip,njp);
				rhon.rho = calculateRho(nin,njn);
				
				
// 				printf("LAMBDA: %f %f %f %f\n",lambdaip.rotation,lambdain.rotation,lambdajp.rotation,lambdajp.rotation);
				
				
				PHIp = calculatePHI(taxp, nip, njp, lambdaip, lambdajp, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi);
				PHIn = calculatePHI(taxn, nin, njn, lambdain, lambdajn, idi, idj, circles[idi].form, circles[idj].form, rhon, indexi);
				etap = calculateEta(taxp, nip, njp, lambdaip, lambdajp, idi, idj, circles[idi].form, circles[idj].form, rhop).rotation;
				etan = calculateEta(taxn, nin, njn, lambdain, lambdajn, idi, idj, circles[idi].form, circles[idj].form, rhon).rotation;
				omegap = calculateOmega(taxp, nip, njp, idi, idj, circles[idi].form, circles[idj].form, rhop, indexi).rotation;
				omegan = calculateOmega(taxn, nin, njn, idi, idj, circles[idi].form, circles[idj].form, rhon, indexi).rotation;
				
// 				printf("eta: %f %f\n",etap,etan);
// 				printf("omega: %f %f\n",omegap,omegan);
				
				
				fd_dPHIij_dxt_out(j) = (PHIp.out.rotation - PHIn.out.rotation) / (2*FD);
				fd_dPHIij_dxt_in(j) = (PHIp.in.rotation - PHIn.in.rotation) / (2*FD);
				fd_detaij_dxt_out(j) = (etap - etan) / (2*FD);
				fd_detaij_dxt_in(j) = (etap - etan) / (2*FD);
				fd_domegaij_dxt_out(j) = (omegap - omegan) / (2*FD);
				fd_domegaij_dxt_in(j) = (omegap - omegan) / (2*FD);
			}
			printf("dPHIij_dxt: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",PHI.drotation_dxt(0),PHI.drotation_dxt(1),PHI.drotation_dxt(2), fd_dPHIij_dxt_in(0), fd_dPHIij_dxt_in(1), fd_dPHIij_dxt_in(2), fd_dPHIij_dxt_out(0), fd_dPHIij_dxt_out(1), fd_dPHIij_dxt_out(2));
			printf("detaij_dxt: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",eta.drotation_dxt(0),eta.drotation_dxt(1),eta.drotation_dxt(2), fd_detaij_dxt_in(0), fd_detaij_dxt_in(1), fd_detaij_dxt_in(2), fd_detaij_dxt_out(0), fd_detaij_dxt_out(1), fd_detaij_dxt_out(2));
			printf("domegaij_dxt: (%f %f %f) FD: (%f, %f, %f) FD: (%f %f %f)\n",omega.drotation_dxt(0),omega.drotation_dxt(1),omega.drotation_dxt(2), fd_domegaij_dxt_in(0), fd_domegaij_dxt_in(1), fd_domegaij_dxt_in(2), fd_domegaij_dxt_out(0), fd_domegaij_dxt_out(1), fd_domegaij_dxt_out(2));
			for(int i=0; i<3; ++i){
 				if(abs(PHI.drotation_dxt(i) - fd_dPHIij_dxt_in(i)) > FDT && abs(PHI.drotation_dxt(i) - fd_dPHIij_dxt_out(i)) > FDT) printf("ERROR\n");//printf("DIVERGENCE\n");
				if(abs(omega.drotation_dxt(i) - fd_domegaij_dxt_in(i)) > FDT && abs(omega.drotation_dxt(i) - fd_domegaij_dxt_out(i)) > FDT) printf("ERROR\n");//printf("DIVERGENCE\n");
			}
		}		
		}
	
	
	*/
	
	
	return PHI;
	
			
}









Rotation Tessellation::calculatePsi(TessellationAxis &tessellationAxis, CircularInterface &circle){
	Rotation psi;
	psi = calculatePsi(tessellationAxis, circle.normal, circle.dmu_dx, circle.form, circle.index);
	
	
	//printf("TAX: (%f %f %f) %d [%d : %d] %d\n",tessellationAxis.v(0),tessellationAxis.v(1),tessellationAxis.v(2),tessellationAxis.hemisphere,tessellationAxis.index, circle.index, circle.form);
	
	/*
	Vector ai,al,achi,x,mu,chi,v(3),chip,chin,xp,xn,mup,mun;
	float d,pos,neg;
	float r_i,r_l,r_chi,g_chi, g;
	if(circle.form != SPLITTER){
		if(tessellationAxis.mode==ALIGNED){
			r_i = radii[circle.index];
			r_l = tradius;
			r_chi = radii[tessellationAxis.index];
			
			
			
			ai = atoms[circle.index];
			al = torigin;
			achi = atoms[tessellationAxis.index];
			if(tessellationAxis.hemisphere==FRONTHEMISPHERE){
				chi=calculateInterfaceNormal(al, achi, d);
			}
			else{
				chi=-calculateInterfaceNormal(al, achi, d);
			}
			
			g_chi = (d * d + r_l * r_l - r_chi * r_chi ) / (2 * d * r_l);

			if(g_chi<0){
				chi=-chi;
			}
			
			
			


			//dxi
			x=ai;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				mup=calculateInterfaceNormal(al, xp, d);
				mun=calculateInterfaceNormal(al, xn, d);
				
				if(circle.form==CONCAVE){
					mup=-mup;
					mun=-mun;
				}
				
				pos=acos(norm_dot(mup,chi));
				neg=acos(norm_dot(mun,chi));
				
				v(j) = (pos - neg) / (2*FD);
			}
			
			
			printf("DXI: (%f %f %f) NUM: (%f %f %f)\n", psi.drotation_dxi(0), psi.drotation_dxi(1), psi.drotation_dxi(2), v(0), v(1), v(2));
			for(int i=0; i<3; ++i){
 				if(abs(psi.drotation_dxi(i) - v(i)) > FDT){
					printf("ERROR\n");
					exit(-1);
				}
			}
		
			
			
			//dxl
			x=al;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				mup=calculateInterfaceNormal(xp, ai, d);
				mun=calculateInterfaceNormal(xn, ai, d);
				
				if(circle.form==CONCAVE){
					mup=-mup;
					mun=-mun;
				}
				
				
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
					chip=calculateInterfaceNormal(xp, achi, d);
				else
					chip=-calculateInterfaceNormal(xp, achi, d);
				
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
					chin=calculateInterfaceNormal(xn, achi, d);
				else
					chin=-calculateInterfaceNormal(xn, achi, d);
				
				if(g_chi<0){
					chip=-chip;
					chin=-chin;
				}
				
				
				
				pos=acos(norm_dot(mup,chip));
				neg=acos(norm_dot(mun,chin));
				
				v(j) = (pos - neg) / (2*FD);
			}
			
			
			printf("DXL: (%f %f %f) NUM: (%f %f %f)\n", psi.drotation_dxl(0), psi.drotation_dxl(1), psi.drotation_dxl(2), v(0), v(1), v(2));
			for(int i=0; i<3; ++i){
 				if(abs(psi.drotation_dxl(i) - v(i)) > FDT){
					printf("ERROR\n");
					exit(-1);
				}
			}
			
			
			
			//dxchi
			x=achi;
			for(int j=0; j<3; ++j){
				xp=x;
				xp(j) = x(j) + FD;
				xn=x;
				xn(j) = x(j) - FD;
				
				
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
					chip=calculateInterfaceNormal(al, xp, d);
				else
					chip=-calculateInterfaceNormal(al, xp, d);
				
				if(tessellationAxis.hemisphere==FRONTHEMISPHERE)
					chin=calculateInterfaceNormal(al, xn, d);
				else
					chin=-calculateInterfaceNormal(al, xn, d);
				
				
				if(g_chi<0){
					chip=-chip;
					chin=-chin;
				}
				
				mup=calculateInterfaceNormal(al, ai, d);
				mun=calculateInterfaceNormal(al, ai, d);
				
				if(circle.form==CONCAVE){
					mup=-mup;
					mun=-mun;
				}
				
				
				pos=acos(norm_dot(mup,chip));
				neg=acos(norm_dot(mun,chin));
				
				v(j) = (pos - neg) / (2*FD);
			}
			
			
			printf("DXCHI: (%f %f %f) NUM: (%f %f %f)\n", psi.drotation_dxt(0), psi.drotation_dxt(1), psi.drotation_dxt(2), v(0), v(1), v(2));
			for(int i=0; i<3; ++i){
 				if(abs(psi.drotation_dxt(i) - v(i)) > FDT){
					printf("ERROR\n");
					exit(-1);
				}
			}
			
		}
	}
	*/
	
	return psi;
}

Rotation Tessellation::calculatePsi(TessellationAxis &tessellationAxis, Vector &mu_i){
	Matrix dmu_dx(3,3);
	Matrix dchi_dx(3,3);
	CircularInterfaceForm form;
	int id;
	
	dmu_dx.zeros();
	dchi_dx.zeros();
	form = CONVEX;
	id=0;
	
	return calculatePsi(tessellationAxis, mu_i, dmu_dx,form, id);
}

Rotation Tessellation::calculatePsi(TessellationAxis &tessellationAxis, Vector &mu_i, Matrix &dmu_dx, CircularInterfaceForm form, int index){
	Vector dpsi_dmui;
	Rotation r;
	Vector dpsi_dxi, dpsi_dxl, dpsi_dchi, dpsi_dxchi;
	
	Vector n_i, nn_i;
	Vector ey(3);
	Matrix I(3,3);
	I.eye();
	
	Matrix dni_dmui;
	Matrix dnni_dni;
	Vector dpsi_dnni;

	Matrix dni_dxi;
	Matrix dnni_dxi;
	float dot_mui_chi;
	
	
	r.drotation_dxi = Vector(3).zeros();
	r.drotation_dxj = Vector(3).zeros();
	r.drotation_dxl = Vector(3).zeros();
	r.drotation_dxt = Vector(3).zeros();
	
	//printf("ANGLE: %d %d\n",tessellationAxis.v.size(), mu_i.size());
	
	r.rotation = getAngle(tessellationAxis.v,mu_i);
	

	if(form != SPLITTER){
		
		
		dot_mui_chi = dot(mu_i,tessellationAxis.v);

		
		if(tessellationAxis.mode==ALIGNED){
			if(index==tessellationAxis.index){
				dpsi_dmui = Vector(3).zeros();
				dpsi_dxi = Vector(3).zeros();
				dpsi_dxchi = Vector(3).zeros();
				dpsi_dxl = Vector(3).zeros();
			}
			else{
				
				dpsi_dmui = -tessellationAxis.v/(sqrt(1-dot_mui_chi*dot_mui_chi));
				dpsi_dchi = -mu_i/(sqrt(1-dot_mui_chi*dot_mui_chi));
				
				
				
				dpsi_dxi = (dpsi_dmui.t() * dmu_dx).t();
				dpsi_dxchi = (dpsi_dchi.t() * tessellationAxis.dchi_dx).t();
				//smoothing
				//if(norm(dpsi_dxi,2)>1.0) dpsi_dxi = normalise(dpsi_dxi);
				
				dpsi_dxl = -dpsi_dxi-dpsi_dxchi;
			}
			
			r.drotation_dxi = dpsi_dxi;
			r.drotation_dxl = dpsi_dxl;
			r.drotation_dxt = dpsi_dxchi;
			
			if(isnan(r.drotation_dxi(0)) || isnan(r.drotation_dxi(1)) || isnan(r.drotation_dxi(2))) printf("DXI IS NAN\n");
			if(isnan(r.drotation_dxt(0)) || isnan(r.drotation_dxt(1)) || isnan(r.drotation_dxt(2))) printf("DXJ IS NAN\n");
			if(isnan(r.drotation_dxl(0)) || isnan(r.drotation_dxl(1)) || isnan(r.drotation_dxl(2))) printf("DXL IS NAN\n");
			
			if(isnan(r.drotation_dxi(0)) || isnan(r.drotation_dxi(1)) || isnan(r.drotation_dxi(2)) ||
			isnan(r.drotation_dxt(0)) || isnan(r.drotation_dxt(1)) || isnan(r.drotation_dxt(2)) || 
			isnan(r.drotation_dxl(0)) || isnan(r.drotation_dxl(1)) || isnan(r.drotation_dxl(2))){
				printf("TAX: %f %f %f\n",tessellationAxis.v(0),tessellationAxis.v(1),tessellationAxis.v(2));
				printf("mui: %f %f %f\n",mu_i(0),mu_i(1),mu_i(2));
				printf("dpsi_dmui: %f %f %f\n",dpsi_dmui(0),dpsi_dmui(1),dpsi_dmui(2));
				printf("dpsi_dchi: %f %f %f\n",dpsi_dchi(0),dpsi_dchi(1),dpsi_dchi(2));
				printf("index %d tax %d\n",index,tessellationAxis.index);
			}
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		
		}
		else{
			
			if(isWithinNumericalLimits(dot_mui_chi,1.0) || isWithinNumericalLimits(dot_mui_chi,-1.0)){
				dpsi_dmui = Vector(3).zeros();
				dpsi_dxi = Vector(3).zeros();
				dpsi_dxl = Vector(3).zeros();
			}
			else{
				
				//if(mu_i(0) > 1.0-MINISCULE) mu_i(0) = 1.0-MINISCULE;
				//if(mu_i(0) < -1.0+MINISCULE) mu_i(0) = -1.0+MINISCULE;
				
				if(dot_mui_chi > 1.0-MINISCULE) dot_mui_chi = 1.0-MINISCULE;
				if(dot_mui_chi < -1.0+MINISCULE) dot_mui_chi = -1.0+MINISCULE;
				
				dpsi_dmui = -tessellationAxis.v/(sqrt(1-dot_mui_chi*dot_mui_chi));
				
				
				dpsi_dxi = (dpsi_dmui.t() * dmu_dx).t();
				//smoothing
				//if(norm(dpsi_dxi,2)>1.0) dpsi_dxi = normalise(dpsi_dxi);
				
				dpsi_dxl = -dpsi_dxi;
			}
			
			
			r.drotation_dxi = dpsi_dxi;
			r.drotation_dxj = Vector(3).zeros();
			r.drotation_dxl = dpsi_dxl;
		}
	}
	else{
		//the splitter is immovable
		//r.drotation_dxi = Vector(3).zeros();
		//r.drotation_dxj = Vector(3).zeros();
		//r.drotation_dxl = Vector(3).zeros();
	}
	

	
	
	
	return r;
	
}




IntersectionBranches::iterator Tessellation::increaseBranchInterator(IntersectionBranches::iterator it, CircularInterface &circle){
	//IntersectionBranches* p;
	//p = it->second.body;
	if(circle.form == CONVEX){
		++it;
		if(it == circle.intersectionBranches.end()) it=circle.intersectionBranches.begin();
		
	}
	else{
		if(it == circle.intersectionBranches.begin()) it=circle.intersectionBranches.end();
		--it;
	}
	
	return it;
	
}

IntersectionBranches::iterator Tessellation::decreaseBranchInterator(IntersectionBranches::iterator it, CircularInterface &circle){
	//IntersectionBranches* p;
	//p = it->second.body;
	
	if(circle.form == CONVEX){
		if(it == circle.intersectionBranches.begin()) it=circle.intersectionBranches.end();
		--it;
		
	}
	else{
		++it;
		if(it == circle.intersectionBranches.end()) it=circle.intersectionBranches.begin();
	}
	
	return it;
}






void Tessellation::createIntersectionBranch(PHIContainer &PHII, CircularInterface &I, CircularInterface &J, RhoContainer &rho){
	
	pair<float, IntersectionBranch*> x;
	IntersectionBranches::iterator it0;
	IntersectionBranches::iterator it1;
	
	
	
	x.first = PHII.in.rotation;
	x.second = new IntersectionBranch();
	x.second->visited = 0;
	x.second->visited1 = 0;
	x.second->direction = IN;
	x.second->body = &I.intersectionBranches;
	x.second->id = J.id;
	x.second->flagged = false;
	x.second->PHI = PHII.in;
	x.second->rho=rho;
	x.second->weight=0;
	x.second->weight1=0;
	x.second->i=-1;
	it0 = I.intersectionBranches.insert(x);
	
	x.first = PHII.out.rotation;
	x.second = new IntersectionBranch();
	x.second->visited = 0;
	x.second->visited1 = 0;
	x.second->direction = OUT;
	x.second->body = &I.intersectionBranches;
	x.second->id = J.id;
	x.second->flagged = false;
	x.second->PHI = PHII.out;
	x.second->rho=rho;
	x.second->weight=0;
	x.second->weight1=0;
	x.second->i=-1;
	it1 = I.intersectionBranches.insert(x);

	
		
	
	
}

void Tessellation::createIntersectionBranch(PHIContainer &PHII, CircularInterface &I, CircularInterface &J, RhoContainer &rho, IntersectionBranch &b){
	
	pair<float, IntersectionBranch*> x;
	IntersectionBranches::iterator it0;
	IntersectionBranches::iterator it1;
	

	
	if(b.direction==IN){
		x.first = PHII.in.rotation;
		x.second = new IntersectionBranch();
		x.second->visited = 0;
		x.second->visited1 = 0;
		x.second->direction = IN;
		x.second->body = &I.intersectionBranches;
		x.second->id = J.id;
		x.second->flagged = false;
		x.second->PHI = PHII.in;
		x.second->rho=rho;
		x.second->v0=b.v0;
		x.second->v1=b.v1;
		x.second->weight=b.weight;
		x.second->weight1=b.weight1;
		x.second->i=b.i;
		
		it0 = I.intersectionBranches.insert(x);
	}
	
	if(b.direction==OUT){
		x.first = PHII.out.rotation;
		x.second = new IntersectionBranch();
		x.second->visited = 0;
		x.second->visited1 = 0;
		x.second->direction = OUT;
		x.second->body = &I.intersectionBranches;
		x.second->id = J.id;
		x.second->flagged = false;
		x.second->PHI = PHII.out;
		x.second->rho=rho;
		x.second->v0=b.v0;
		x.second->v1=b.v1;
		x.second->weight=b.weight;
		x.second->weight1=b.weight1;
		x.second->i=b.i;
		it1 = I.intersectionBranches.insert(x);
	}

	
		
	
	
}







void Tessellation::printBranch(const char* s, multimap<float, IntersectionBranch>::iterator &it){
	if(it->second.direction == OUT)
		printf("iterator %s: %d-%d (OUT)\n",s,it->second.node->id0, it->second.node->id1);
	else
		printf("iterator %s: %d-%d (IN)\n",s,it->second.node->id0, it->second.node->id1);
}

void Tessellation::printIntersectionGraph(IntersectionGraph &g,CircularInterfacesPerAtom &circles){
	IntersectionGraph::iterator it;
	
	for(it = g.begin(); it != g.end(); ++it){
		//fprintf(stderr,"NODE[%d,%d]: (%d,%d) -> (%d,%d)\n", circles[it->first.id0].index, circles[it->first.id1].index, circles[it->second.id0].index, circles[it->second.id1].index, circles[it->second.pointsTo0-1].index, circles[it->second.pointsTo1-1].index);
		fprintf(stderr,"NODE[%d,%d]: (%d,%d) -> (%d,%d)\n", it->first.id0, it->first.id1, it->second.id0, it->second.id1, it->second.pointsTo0, it->second.pointsTo1);
	}
}


float Tessellation::V2PHI(TessellationAxis &tessellationAxis, Hemisphere hemisphere, Vector v, Vector &normal, float g){
	Vector ctsl;
	Vector ncaux(3);
	Vector nc(3);
	Vector cv;
	float dot_nc_cv;
	float a;
	Vector p(3);
	
	


	
	if(isWithinNumericalLimits(dot(tessellationAxis.v,normal),1.0)){
		//printf("special treatment\n");
		nc(0) = 0;
		nc(1) = -1;
		nc(2) = 0;
		
		ncaux=cross(nc,tessellationAxis.v);
		
	}
	else{
		ctsl=tessellationAxis.v-normal;
		ncaux = cross(ctsl,normal);
		nc = cross(normal,ncaux);
	}
	
	
	cv=v-normal*g;
	
	dot_nc_cv = norm_dot(nc,cv);
	
	a = acos(dot_nc_cv);
	if(dot(cv,ncaux)<0) a = -a;
	
	return a;
	
	

	
}

void Tessellation::convertExposedVectors2PHIValues(TessellationAxis &tessellationAxis, Hemisphere hemisphere, CircularInterface &circle){
	for(unsigned int i=0; i<circle.exposedVectors.size(); ++i){
		//if(circle.form==SPLITTER || dot(circle.exposedVectors[i], tessellationAxis)>=THRESHOLD_NUMERICAL)
			circle.exposedPHI[hemisphere].push_back(V2PHI(tessellationAxis, hemisphere, circle.exposedVectors[i], circle.normal, circle.lambdaRotation.g));
	}
}




float Tessellation::exposition(Hemisphere hemisphere, IntersectionBranches::iterator it0, IntersectionBranches::iterator it1, CircularInterface &circle){
	float e;
	float m,mmax;
	float r0,r1;
	IntersectionBranches::iterator it;
	
	m=0;
	mmax=numeric_limits<float>::max();
	if(circle.exposedPHI[hemisphere].size()==0){
		//printf("BIIIIG PROBLEM %d\n",circle.index);
		return M_PI;
		
	}
	
	
	
	
	
	
	
	
	for(unsigned int i=0; i<circle.exposedPHI[hemisphere].size(); ++i){
			e = circle.exposedPHI[hemisphere][i];
			
			
			it = circle.intersectionBranches.upper_bound(e);
			if(it==circle.intersectionBranches.end())
				it=circle.intersectionBranches.begin();
			
			
			r0=abs(it0->first-e);
			if(r0>M_PI) r0=2*M_PI-r0;
			r1=abs(it1->first-e);
			if(r1>M_PI) r1=2*M_PI-r1;
			
			
			if(circle.form==CONVEX){
				if(it==it1) m=0;
				else{
					m=min(r0,r1);
				}
			}
			else{
				if(it==it0) m=0;
				else{
					m=min(r0,r1);
				}
			}
			
				
			
			mmax = min(mmax,m);
			
			
			
	}
	
	return mmax;
	
}





/**
 * This method visits all circular interfaces. Within each interface, it cycles through all intersections with other interfaces and
 * adds them to a structure called intersection graph.  This structure contains angular and connectivity information for the whole 
 * sphere (disregarding front hemisphere and back hemisphere boundaries)
 */
void Tessellation::buildIntersectionGraphFirstPass(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles){
	unsigned int j;
	CircularInterface *I, *J;
	map<int,RhoContainer>::iterator it_j;
	
	
	IntersectionBranches::iterator it;
	PHIContainer  PHII;
	
	bool erased;
	bool done;
	RhoContainer rhoContainer;
	int round;
	map<int,bool> activeBranches;
	
	
	//iterate through all circles and add them to the intersectiongraph, one by one
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		
		if(I->circularIntersections.size()==0) continue;
		
		//iterate through all intersections with other interfaces
		for(it_j=I->circularIntersections.begin(); it_j != I->circularIntersections.end(); ++it_j){
			j=it_j->second.id;
			rhoContainer = it_j->second;
			J = &circles[j];
			
			//retrieve external and internal interfaces and create a branch stored in PHII (which stands for PHI, circle I)
			PHII = calculatePHI(tessellationAxis, *I, *J, radius, rhoContainer);
			createIntersectionBranch(PHII, *I, *J,rhoContainer);
		}
		
		//now we have to do 2 rounds of eliminations
		round=0;
		activeBranches.clear();
		
		it = I->intersectionBranches.begin();
		done=false;
		while(!done){
			//for out directions, we have to check if the ip should be deleted before we add the interface to the active branches
			if(it->second->direction==OUT){
				++it->second->visited;
				if(activeBranches.size()>0) it->second->flagged=true;
				activeBranches[it->second->id]=true;				
			}
			//for in directions the other way around
			else{
				activeBranches.erase(it->second->id);
				if(activeBranches.size()>0) it->second->flagged=true;
			}
			
			it=increaseBranchInterator(it, *I);

			if(round==0 && it->second->visited==1) ++round;
			else if(round==1 && it->second->visited==2) done=true;

			
			
		}
			

		
		
		//The branches have been flagged to be deleted, here, we perform the actual deletion.
		it = I->intersectionBranches.begin();
		while(it != I->intersectionBranches.end()){
			erased=false;
			if(it->second->flagged){
				//dist0 = buriedness(hemisphere, *I, circles[it->second.id], it->second.PHI.rotation, buffer0, buffer1);
				//fprintf(stderr,"%f\n",dist0);
				it=I->intersectionBranches.erase(it);
				erased=true;
			}
			
			if(!erased) ++it;
				
		}
			
	}
}





/**
 * The splitter is a newly placed artificial interface which might, or might not, intersect with other interfaces.
 * We need to calculate its intersection points which is what this method is for.
 * */
void Tessellation::buildIntersectionGraphSplitterPass(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles){
	unsigned int j;
	CircularInterface *I, *J;
	map<int,RhoContainer>::iterator it_j;
	
	
	IntersectionBranches::iterator it;
	PHIContainer  PHII;
	
	bool erased;
	bool done;
	RhoContainer rhoContainer;
	int round;
	map<int,bool> activeBranches;
	unsigned int activeBranches2;
	bool iterateAll;
	
	
	//iterate through all circles and add them to the intersectiongraph, one by one
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		
	
		
		
		
		if(I->circularIntersections.size()==0) continue;
		
		
		
		
		it_j=I->circularIntersections.find(-1);
	
		iterateAll=false;
		if(I->form==SPLITTER){
			iterateAll=true;
			it_j=I->circularIntersections.begin();
		}
		
		if(it_j==I->circularIntersections.end()) continue;
			
		do{


			j=it_j->second.id;
			rhoContainer = it_j->second;
			J = &circles[j];
			
			//retrieve external and internal interfaces (respectively)
			//PHIJ = calculatePHI(tessellationAxis, *J, *I, radius, rhoContainer);
			PHII = calculatePHI(tessellationAxis, *I, *J, radius, rhoContainer);
			
			createIntersectionBranch(PHII, *I, *J,rhoContainer);
			
			++it_j;

		}
		while(iterateAll && it_j != I->circularIntersections.end());
		
	
		
		
		
		
		//now we have to do 2 rounds of eliminations
		round=0;
		activeBranches.clear();
		activeBranches2=0;
		it = I->intersectionBranches.begin();
		done=false;
		while(!done){
			//for out directions, we have to check if the ip should be deleted before we add the interface to the active branches
			if(it->second->direction==OUT){
				++it->second->visited1;
				if(activeBranches.size()>0 || activeBranches2>0){
					it->second->flagged=true;
				}
				if(I->form==SPLITTER || circles[it->second->id].form==SPLITTER) activeBranches[it->second->id]=true;
				else activeBranches2++;
			}
			//for in directions the other way around
			else{
				if(I->form==SPLITTER || circles[it->second->id].form==SPLITTER) activeBranches.erase(it->second->id);
				else activeBranches2=max((int)0,(int)activeBranches2-1); //shouldn't go lower than zero
				
				
				if(activeBranches.size()>0 || activeBranches2>0){
					it->second->flagged=true;
				}
			}
			
			it=increaseBranchInterator(it, *I);

			if(round==0 && it->second->visited1==1) ++round;
			else if(round==1 && it->second->visited1==2) done=true;

			
			
		}
		

		
		//The branches have been flagged to be deleted, here, we perform the actual deletion.
		it = I->intersectionBranches.begin();
		while(it != I->intersectionBranches.end()){
			erased=false;
			if(it->second->flagged){
				//dist0 = buriedness(hemisphere, *I, circles[it->second.id], it->second.PHI.rotation, buffer0, buffer1);
				//fprintf(stderr,"%f\n",dist0);
				it=I->intersectionBranches.erase(it);
				erased=true;
			}
			
			if(!erased) ++it;
				
		}
			
	}		
	
}






/**
 * This method removes invalid splitter intersectionpoints
 * */
void Tessellation::splitterSanityCheck(CircularInterfacesPerAtom &circles){
	int splitter=0;
	IntersectionBranches::iterator it;
	map<SplitterIntersection,bool,SplitterIntersectionComparator> splitterIntersections;
	CircularInterface *I;
	bool erased;
	SplitterIntersection s;
	
	//In a first pass, we traverse through the interfaces and register intersections with the splitter + id of the splitter
	for(unsigned int i=0; i<circles.size(); ++i){
		if(circles[i].form==SPLITTER){
			splitter=i;
		}
		else{
			I = &circles[i];
			for(it = I->intersectionBranches.begin(); it != I->intersectionBranches.end(); ++it){
				if(circles[it->second->id].form==SPLITTER){
					if(it->second->direction==IN) s.direction=OUT;
					else s.direction=IN;
					s.id=i;
					splitterIntersections[s]=true;
				} 
			}

		}
	}
	I = &circles[splitter];
	//In a second pass we traverse through the splitter's branches and check for validity
	it = I->intersectionBranches.begin();
	while(it != I->intersectionBranches.end()){
		erased=false;
		s.id=it->second->id;
		s.direction=it->second->direction;
		if(splitterIntersections.find(s)!=splitterIntersections.end()){
			//do nothing, the branch is valid
		}
		else{
			//invalid branch
			it=I->intersectionBranches.erase(it);
			erased=true;
		}
		if(!erased) ++it;
	}
	
	
}


/**
 * This method creates a duplicate of the entire intersectiongraph.
 * It involves a recalculation of all phi angles for the new graph, with respect to a new tessellation axis.
 * */
void Tessellation::copyIntersectionGraph(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, CircularInterfacesPerAtom &newCircles){
	SASASegment sasaSegment;
	
	
	CircularInterface *I, *J;
	vector<RhoContainer>::iterator it_j;
	
	
	IntersectionBranches::iterator it;
	PHIContainer PHII;
	
	RhoContainer rhoContainer;

	//clear intersectionbranches. They will be refreshed in the next step.
	newCircles=CircularInterfacesPerAtom(circles);
	for(unsigned int i=0; i < newCircles.size(); ++i){
		newCircles[i].intersectionBranches.clear();
		newCircles[i].hasDerivatives=false;
	}
	
	
	
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		
		//here, the intersection graph is refreshed with the new PHI values with respect to the new tessellation axis
		for(it = I->intersectionBranches.begin(); it != I->intersectionBranches.end(); ++it){
			J=&circles[it->second->id];
			rhoContainer=it->second->rho;
			PHII = calculatePHI(tessellationAxis, *I, *J, radius, rhoContainer);
			
			createIntersectionBranch(PHII, newCircles[i], newCircles[it->second->id], rhoContainer, *(it->second));
			
		}
		

	}

}






/**
 * Some interfaces do not intersect with any other interface. In order to handle them using the same pipeline, we instert artificial "fake"
 * intersection points anywhere on the interface.
 * */
void Tessellation::buildIntersectionGraphArtificialPointsPass(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, unsigned int &globalSegmentCounter){
	SASASegment sasaSegment;
	
	CircularInterface *I;
	
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		
		//if the interface has no intersections with other interfaces, it needs to be treated specially
		if(I->circularIntersections.size()==0 && I->form!=SPLITTER){
			insertArtificialIntersectionPoints(*I,tessellationAxis,hemisphere,sasa, globalSegmentCounter);			
		}
	}
}



void Tessellation::sortGaussBonnetPaths(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, string filename, MultiLayeredDepthBuffer &buffer0, MultiLayeredDepthBuffer &buffer1, bool useDepthBuffer, bool split, unsigned int &globalSegmentCounter, bool derivatives){
	SASASegment sasaSegment;
	
	CircularInterface *I;
	vector<RhoContainer>::iterator it_j;
	
	
	Interfaces interfacesJ, interfacesI;
	IntersectionBranches::iterator it, it2;
	PHIContainer PHIJ, PHII;
	
	SegmentGraph segmentGraph0;
	SegmentGraph segmentGraph1;
	pair<PartialSegmentID,SegmentList::iterator> p;
	SegmentGraph::iterator it_s0, it_s1;
	PartialSegmentID ps0,ps1,pempty;
	bool segmentSearchForward;
	SegmentList segmentList;
	SegmentList::iterator it_sl, it_sl2, it_sl2_start;
	SegmentInfo seginfo;
	vector<vector<SegmentList::iterator> > segmentPointerLists;
	vector<SegmentList::iterator> segmentPointers;
	vector<vector<SegmentList::iterator> >::iterator it_sp;
	int segmentCount;
	int id1;
	
	unsigned int k;
	Vector v0,v1;
	float s_direction;
	float arclength;
	float PHI0,PHI1;
	float kappa;
	
	//count segments;
	segmentCount=0;
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		for(it = I->intersectionBranches.begin(); it != I->intersectionBranches.end(); ++it){
			if(it->second->direction==IN){
				segmentCount++;
			}
		}
	}
	

	//build segmentgraph and assign distances
	segmentList.clear(),
	segmentList.reserve(segmentCount*2);
	
	pempty.i0=-1;
	pempty.i1=-1;
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		if(I->intersectionBranches.size()>0){
			for(it = I->intersectionBranches.begin(); it != I->intersectionBranches.end(); ++it){
				if(it->second->direction==IN){
					it2=increaseBranchInterator(it,*I);
					
					ps0.i0 = it->second->id;
					ps0.i1 = I->id;
					
					ps1.i0 = I->id;
					ps1.i1 = it2->second->id;
					
					
					
					//printf("EXPOSING %d %d %d\n",circles[it->second.id].index,circles[I->id].index,circles[it2->second.id].index);
					
					seginfo.id0 = ps0;
					seginfo.id1 = ps1;
					seginfo.backw=pempty;
					seginfo.forw=pempty;
					seginfo.hasForward=false;
					seginfo.hasBackward=false;
					seginfo.PHI0 = it->second->PHI;
					seginfo.PHI1 = it2->second->PHI;
					seginfo.visited=false;
					seginfo.source=it->second;
					if(hasDepthBuffer && useDepthBuffer) seginfo.dist=exposition(hemisphere, it, it2, *I);
					it_sl = segmentList.insert(segmentList.end(), seginfo);
					
					p.first=ps0;
					p.second=it_sl;
					segmentGraph0.insert(p);
					
					p.first=ps1;
					p.second=it_sl;
					segmentGraph1.insert(p);
					
					it_s0 = segmentGraph0.find(ps1);
					if(it_s0!=segmentGraph0.end()){
						it_sl2 = it_s0->second;
						it_sl->it_forw=it_sl2;
						it_sl->hasForward=true;
						it_sl2->it_backw=it_sl;
						it_sl2->hasBackward=true;
					}
					
					it_s1 = segmentGraph1.find(ps0);
					if(it_s1!=segmentGraph1.end()){
						it_sl2 = it_s1->second;
						it_sl->it_backw=it_sl2;
						it_sl->hasBackward=true;
						it_sl2->it_forw=it_sl;
						it_sl2->hasForward=true;
					}
					
					
					
				}
			}
		}
	}
	
	
	
	//create a vector of SASAs
	segmentPointerLists.clear();
		
		
	for(it_sl=segmentList.begin(); it_sl!=segmentList.end(); ++it_sl){
		segmentPointers.clear();
		segmentSearchForward=true;
		it_sl2 = it_sl;
		it_sl2_start=it_sl2;
		while(it_sl2!=segmentList.end() && !it_sl2->visited){
			segmentPointers.push_back(it_sl2);
			it_sl2->visited=true;
			if(segmentSearchForward){
				if(it_sl2->hasForward) it_sl2=it_sl2->it_forw;
				else {
					//it's not guaranteed that the sasa is a cycle. 
					//however, we get better results disregarding broken cycles
					segmentSearchForward=false;
					it_sl2=segmentList.end();
					segmentPointers.clear();
					
					
					it_sl2=it_sl2_start;
				}
			}
			
			if(!segmentSearchForward){
				if(it_sl2->hasBackward) it_sl2=it_sl2->it_backw;
				else {
					it_sl2=segmentList.end();

				}
			}
			
			
			
		}
		if(segmentPointers.size()>0) segmentPointerLists.push_back(segmentPointers);
	}			
	
	
	//add intersection-points and semi-segment points
	for(unsigned int i=0; i<segmentPointerLists.size(); ++i){
		for(unsigned int j=0; j<segmentPointerLists[i].size(); ++j){
				id1 = segmentPointerLists[i][j]->id0.i1;
			
				
				PHI0=segmentPointerLists[i][j]->PHI0.rotation;
				PHI1=segmentPointerLists[i][j]->PHI1.rotation;
				
				if(PHI1 >= PHI0) s_direction=1;
				else s_direction=-1;
				arclength = s_direction*PHI1 - s_direction*PHI0;
				
 				if(PHI1 < PHI0)
 					arclength = 2*M_PI - arclength;
				
				if(circles[id1].form!=CONVEX)
					arclength = 2*M_PI - arclength;
				
				arclength=arclength*sin(circles[id1].lambda.rotation);			
				
				kappa=calculateKappa(tessellationAxis, circles[id1].normal, circles[id1].form);
				v0=PHI2V(tessellationAxis,PHI0, circles[id1].psi.rotation, circles[id1].lambda.rotation, kappa, circles[id1].form);
				v1=PHI2V(tessellationAxis,PHI0+arclength/2.0, circles[id1].psi.rotation, circles[id1].lambda.rotation, kappa, circles[id1].form);
				
				k=j+1;
				if(k>=segmentPointerLists[i].size()) k=0;
				
				segmentPointerLists[i][j]->source->v0=v0;
				segmentPointerLists[i][j]->source->v1=v1;
				//distribute weights
				segmentPointerLists[i][j]->source->weight+=arclength;
				segmentPointerLists[i][j]->source->weight1=arclength;
				segmentPointerLists[i][k]->source->weight+=arclength;
				segmentPointerLists[i][k]->source->i=i;
				
				
				
				
		}
	}
	
	globalSegmentCounter=segmentPointerLists.size();

	
	

			
	//@TODO WE WILL DO THAT STUFF LATER
	
	//We compare here to the previous tessellation
	//if there are differences, we break calculation and redo without depthbuffer
	
	/*
	set<FullSegmentID,SegmentSetComparator> segmentSet;
	FullSegmentID segmentID;
	if(hasDepthBuffer && useDepthBuffer){
		//printf("NEW SEGMENTS: ");
	
		for(unsigned int i=0; i<segmentPointerLists.size(); ++i){
			for(unsigned int j=0; j<segmentPointerLists[i].size(); ++j){
				id0 = segmentPointerLists[i][j]->id0.i0;
				id1 = segmentPointerLists[i][j]->id0.i1;
				id2 = segmentPointerLists[i][j]->id1.i1;
				
				segmentID.i0 = circles[id0].index;
				segmentID.i1 = circles[id1].index;
				segmentID.i2 = circles[id2].index;
				
				//printf(" (%d %d %d) ",segmentID.i0,segmentID.i1,segmentID.i2);
				
				segmentSet.insert(segmentID);
			}
		}
		
		//printf("\n");
		
		//printf("OLD SEGMENTS: ");
		
		int ac=0;
		int bc=0;
		//now we check if the previous and current tessellation are equivalent
		for(unsigned int i=0; i<prevSasasForMolecule[l].size(); ++i){
			segmentID.i0=prevSasasForMolecule[l][i].index0;
			segmentID.i1=prevSasasForMolecule[l][i].index1;
			segmentID.i2=prevSasasForMolecule[l][i].index2;
			
			if(prevSasasForMolecule[l][i].hemisphere==hemisphere){
			//printf(" (%d %d %d) ",segmentID.i0,segmentID.i1,segmentID.i2);
			
			if(segmentSet.find(segmentID)==segmentSet.end()){
				//return false;		
				//printf("%d ",l);
				return false;
				ac++;
			}else bc++;
			}
		}
		//printf("\n");
		//printf("AC: %d %d\n", ac,bc);
//if(mismatch) return false;
	}
	if(hasDepthBuffer && !useDepthBuffer){
		//printf("REVISED SEGMENTS: ");
		
		for(unsigned int i=0; i<segmentPointerLists.size(); ++i){
			for(unsigned int j=0; j<segmentPointerLists[i].size(); ++j){
				id0 = segmentPointerLists[i][j]->id0.i0;
				id1 = segmentPointerLists[i][j]->id0.i1;
				id2 = segmentPointerLists[i][j]->id1.i1;
				
				segmentID.i0 = circles[id0].index;
				segmentID.i1 = circles[id1].index;
				segmentID.i2 = circles[id2].index;
				
				//printf(" (%d %d %d) ",segmentID.i0,segmentID.i1,segmentID.i2);
				
			}
		}
		//printf("\n");
		
	}
	
*/

	
}




bool Tessellation::buildIntersectionGraphCollectionPass(int l, float radius, TessellationAxis &tessellationAxis, CircularInterfacesPerAtom &circles, SASASegmentList &sasa, Hemisphere hemisphere, string filename, MultiLayeredDepthBuffer &buffer0, MultiLayeredDepthBuffer &buffer1, bool useDepthBuffer, bool split, unsigned int &globalSegmentCounter, bool derivatives){
	SASASegment sasaSegment;
	
	CircularInterface *I;
	vector<RhoContainer>::iterator it_j;
	
	
	Interfaces interfacesJ, interfacesI;
	IntersectionBranches::iterator it, it2;
	PHIContainer PHIJ, PHII;
	
	SegmentGraph segmentGraph0;
	SegmentGraph segmentGraph1;
	pair<PartialSegmentID,SegmentList::iterator> p;
	SegmentGraph::iterator it_s0, it_s1;
	PartialSegmentID ps0,ps1,pempty;
	bool segmentSearchForward;
	SegmentList segmentList;
	SegmentList::iterator it_sl, it_sl2, it_sl2_start;
	SegmentInfo seginfo;
	vector<vector<SegmentList::iterator> > segmentPointerLists;
	vector<SegmentList::iterator> segmentPointers;
	vector<vector<SegmentList::iterator> >::iterator it_sp;
	int segmentCount;
	int id0,id1,id2;
	
	//count segments;
	segmentCount=0;
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		for(it = I->intersectionBranches.begin(); it != I->intersectionBranches.end(); ++it){
			if(it->second->direction==IN){
				segmentCount++;
			}
		}
	}
	

	//build segmentgraph and assign distances
	segmentList.clear(),
	segmentList.reserve(segmentCount*2);
	
	pempty.i0=-1;
	pempty.i1=-1;
	for(unsigned int i=0; i < circles.size(); ++i){
		I = &circles[i];
		if(I->intersectionBranches.size()>0){
			if(hasDepthBuffer && useDepthBuffer) convertExposedVectors2PHIValues(tessellationAxis, hemisphere, *I);
			for(it = I->intersectionBranches.begin(); it != I->intersectionBranches.end(); ++it){
				if(it->second->direction==IN){
					it2=increaseBranchInterator(it,*I);
					
					ps0.i0 = it->second->id;
					ps0.i1 = I->id;
					
					ps1.i0 = I->id;
					ps1.i1 = it2->second->id;
					
					
					
					printf("EXPOSING %d %d %d\n",circles[it->second->id].index,circles[I->id].index,circles[it2->second->id].index);
					
					seginfo.id0 = ps0;
					seginfo.id1 = ps1;
					seginfo.backw=pempty;
					seginfo.forw=pempty;
					seginfo.hasForward=false;
					seginfo.hasBackward=false;
					seginfo.PHI0 = it->second->PHI;
					seginfo.PHI1 = it2->second->PHI;
					seginfo.visited=false;
					seginfo.v0=it->second->v0;
					seginfo.v1=it->second->v1;
					seginfo.weight=it->second->weight;
					seginfo.weight1=it->second->weight1;
					seginfo.i=it->second->i;
					if(hasDepthBuffer && useDepthBuffer) seginfo.dist=exposition(hemisphere, it, it2, *I);
					it_sl = segmentList.insert(segmentList.end(), seginfo);
					
					p.first=ps0;
					p.second=it_sl;
					segmentGraph0.insert(p);
					
					p.first=ps1;
					p.second=it_sl;
					segmentGraph1.insert(p);
					
					it_s0 = segmentGraph0.find(ps1);
					if(it_s0!=segmentGraph0.end()){
						it_sl2 = it_s0->second;
						it_sl->it_forw=it_sl2;
						it_sl->hasForward=true;
						it_sl2->it_backw=it_sl;
						it_sl2->hasBackward=true;
					}
					
					it_s1 = segmentGraph1.find(ps0);
					if(it_s1!=segmentGraph1.end()){
						it_sl2 = it_s1->second;
						it_sl->it_backw=it_sl2;
						it_sl->hasBackward=true;
						it_sl2->it_forw=it_sl;
						it_sl2->hasForward=true;
					}
					
					
					
				}
			}
		}
	}
	
	//create a vector of SASAs
	segmentPointerLists.clear();
	if(hasDepthBuffer && useDepthBuffer){
		for(it_sl=segmentList.begin(); it_sl!=segmentList.end(); ++it_sl){
			segmentPointers.clear();
			segmentSearchForward=true;
			it_sl2 = it_sl;
			it_sl2_start=it_sl2;
			while(it_sl2!=segmentList.end() && !it_sl2->visited){
				segmentPointers.push_back(it_sl2);
				it_sl2->visited=true;
				if(segmentSearchForward){
					if(it_sl2->hasForward) it_sl2=it_sl2->it_forw;
					else {
						//it's not guaranteed that the sasa is a cycle. 
						//however, we get better results disregarding broken cycles
						segmentSearchForward=false;
						//segmentPointers.clear();
						
						
						it_sl2=it_sl2_start;
					}
				}
				
				if(!segmentSearchForward){
					if(it_sl2->hasBackward) it_sl2=it_sl2->it_backw;
					else {
						it_sl2=segmentList.end();

					}
				}
				
				
				
			}
			if(segmentPointers.size()>0) segmentPointerLists.push_back(segmentPointers);
		}
		
		
		//filter SASAs
		it_sp=segmentPointerLists.begin();
		while(it_sp!=segmentPointerLists.end()){
			buffer1.startNewCycle();
			for(unsigned int j=0; j<it_sp->size(); ++j){
				buffer1.addProbe((*it_sp)[j]->dist);
				
			}
			if(!buffer1.isCycleExposed()){
				it_sp=segmentPointerLists.erase(it_sp);
			}
			else ++it_sp;
		}
		
	}
	else{
		
		//we do this here just for the sake of printing extra information that will be used to build sasea tables. This will not be called in normal production runs
		{//if(!split){
		
			for(it_sl=segmentList.begin(); it_sl!=segmentList.end(); ++it_sl){
				segmentPointers.clear();
				segmentSearchForward=true;
				it_sl2 = it_sl;
				it_sl2_start=it_sl2;
				while(it_sl2!=segmentList.end() && !it_sl2->visited){
					segmentPointers.push_back(it_sl2);
					it_sl2->visited=true;
					if(segmentSearchForward){
						if(it_sl2->hasForward) it_sl2=it_sl2->it_forw;
						else {
							//it's not guaranteed that the sasa is a cycle. 
							//however, we get better results disregarding broken cycles
							segmentSearchForward=false;
							
							it_sl2=it_sl2_start;
						}
					}
					
					if(!segmentSearchForward){
						if(it_sl2->hasBackward) it_sl2=it_sl2->it_backw;
						else {
							it_sl2=segmentList.end();

						}
					}
					
					
					
				}
				if(segmentPointers.size()>0) segmentPointerLists.push_back(segmentPointers);
			}			
		}
// 		else{
// 		
// 			
// 			//if there's no depth-buffer information, there is no real sense in figuring out which segments are joined.
// 			segmentPointers.clear();
// 			for(it_sl=segmentList.begin(); it_sl!=segmentList.end(); ++it_sl){
// 				segmentPointers.push_back(it_sl);
// 			}
// 			segmentPointerLists.push_back(segmentPointers);
// 		}
	}
	

	//We compare here to the previous tessellation
	//if there are differences, we break calculation and redo without depthbuffer
	set<FullSegmentID,SegmentSetComparator> segmentSet;
	FullSegmentID segmentID;
	if(hasDepthBuffer && useDepthBuffer){
		//printf("NEW SEGMENTS: ");
	
		for(unsigned int i=0; i<segmentPointerLists.size(); ++i){
			for(unsigned int j=0; j<segmentPointerLists[i].size(); ++j){
				id0 = segmentPointerLists[i][j]->id0.i0;
				id1 = segmentPointerLists[i][j]->id0.i1;
				id2 = segmentPointerLists[i][j]->id1.i1;
				
				segmentID.i0 = circles[id0].index;
				segmentID.i1 = circles[id1].index;
				segmentID.i2 = circles[id2].index;
				
				//printf(" (%d %d %d) ",segmentID.i0,segmentID.i1,segmentID.i2);
				
				segmentSet.insert(segmentID);
			}
		}
		
		//printf("\n");
		
		//printf("OLD SEGMENTS: ");
		
		int ac=0;
		int bc=0;
		//now we check if the previous and current tessellation are equivalent
		for(unsigned int i=0; i<prevSasasForMolecule[l].size(); ++i){
			segmentID.i0=prevSasasForMolecule[l][i].index0;
			segmentID.i1=prevSasasForMolecule[l][i].index1;
			segmentID.i2=prevSasasForMolecule[l][i].index2;
			
			if(prevSasasForMolecule[l][i].hemisphere==hemisphere){
			//printf(" (%d %d %d) ",segmentID.i0,segmentID.i1,segmentID.i2);
			
			if(segmentSet.find(segmentID)==segmentSet.end()){
				//return false;		
				//printf("%d ",l);
				return false;
				ac++;
			}else bc++;
			}
		}
		//printf("\n");
		//printf("AC: %d %d\n", ac,bc);
//if(mismatch) return false;
	}
	if(hasDepthBuffer && !useDepthBuffer){
		//printf("REVISED SEGMENTS: ");
		
		for(unsigned int i=0; i<segmentPointerLists.size(); ++i){
			for(unsigned int j=0; j<segmentPointerLists[i].size(); ++j){
				id0 = segmentPointerLists[i][j]->id0.i0;
				id1 = segmentPointerLists[i][j]->id0.i1;
				id2 = segmentPointerLists[i][j]->id1.i1;
				
				segmentID.i0 = circles[id0].index;
				segmentID.i1 = circles[id1].index;
				segmentID.i2 = circles[id2].index;
				
				//printf(" (%d %d %d) ",segmentID.i0,segmentID.i1,segmentID.i2);
				
			}
		}
		//printf("\n");
		
	}
	
		for(unsigned int i=0; i<segmentPointerLists.size(); ++i){
			for(unsigned int j=0; j<segmentPointerLists[i].size(); ++j){
				id0 = segmentPointerLists[i][j]->id0.i0;
				id1 = segmentPointerLists[i][j]->id0.i1;
				id2 = segmentPointerLists[i][j]->id1.i1;
				
				segmentID.i0 = circles[id0].index;
				segmentID.i1 = circles[id1].index;
				segmentID.i2 = circles[id2].index;
				
				printf(" (%d %d %d) ",segmentID.i0,segmentID.i1,segmentID.i2);
				
			}
		printf(".\n");
		}
		printf("#\n");
				
	for(unsigned int i=0; i<segmentPointerLists.size(); ++i){
		for(unsigned int j=0; j<segmentPointerLists[i].size(); ++j){
				id0 = segmentPointerLists[i][j]->id0.i0;
				id1 = segmentPointerLists[i][j]->id0.i1;
				id2 = segmentPointerLists[i][j]->id1.i1;
				
		
				sasaSegment.i = segmentPointerLists[i][j]->i;
				sasaSegment.i2 = j;
				
				sasaSegment.id0 = id0;
				sasaSegment.id1 = id1;
				sasaSegment.id2 = id2;
				
				
				sasaSegment.index0 = circles[id0].index;
				sasaSegment.index1 = circles[id1].index;
				sasaSegment.index2 = circles[id2].index;
				
				sasaSegment.PHI0=segmentPointerLists[i][j]->PHI0;
				sasaSegment.PHI1=segmentPointerLists[i][j]->PHI1;
				
				sasaSegment.tessellationAxis = tessellationAxis;
				sasaSegment.hemisphere = hemisphere;
				sasaSegment.depthBufferEstimatedArea=depthBufferEstimatedArea;
				sasaSegment.radius=radius;
				
				sasaSegment.normalForCircularInterface = circles[id1].normal;
				sasaSegment.form0 = circles[id0].form;
				sasaSegment.form1 = circles[id1].form;
				sasaSegment.form2 = circles[id2].form;
				
				
 				if(!circles[id0].hasDerivatives) 
					calculateProjectionAndDerivatives(tessellationAxis, circles[id0]);
 				if(!circles[id1].hasDerivatives) 
					calculateProjectionAndDerivatives(tessellationAxis, circles[id1]);
 				if(!circles[id2].hasDerivatives) 
					calculateProjectionAndDerivatives(tessellationAxis, circles[id2]);
				
				sasaSegment.lambda = circles[id1].lambda;
				sasaSegment.psi = circles[id1].psi;
				
				sasaSegment.rotation0.rotation = sasaSegment.PHI0.rotation;
				sasaSegment.rotation1.rotation = sasaSegment.PHI1.rotation;
				

				sasaSegment.v0=segmentPointerLists[i][j]->v0;
				sasaSegment.v1=segmentPointerLists[i][j]->v1;
				sasaSegment.weight=segmentPointerLists[i][j]->weight;
				sasaSegment.weight1=segmentPointerLists[i][j]->weight1;

				if(derivatives){
					
					sasaSegment.rotation0 = calculatePHIDerivatives(sasaSegment.PHI0, circles, tessellationAxis);
					sasaSegment.rotation1 = calculatePHIDerivatives(sasaSegment.PHI1, circles, tessellationAxis);
					
					
				}
				
				sasaSegment.artificial=false;
			
				
				sasa.push_back(sasaSegment);
				
			
			
			
				
				
		}
		globalSegmentCounter++;
		
		
	}
	


	//indicating tessellation is complete
	return true;


	
	
}






void Tessellation::outputGaussBonnetData(string filename, float radius, CircularInterfacesPerAtom &circles, SASAs &sasas, IntersectionGraph &intersectionGraph){
/*	FILE* file;
	
	file = fopen (filename.c_str(),"w");
	
	SASA sasa;
	SASANodeList::iterator it;
	CircularInterface circle;
	float area=0;
	int k;
	Vector origin(3), v(3);
	
	
	fprintf(file, "atom radius %f\n", radius);
	
	
	for(int i=0; i< circles.size(); ++i){
		circle = circles[i];
		if(circle.form==SPLITTER){
			fprintf(file, "circularregion %d radius %f vector %f %f %f form %d\n", circle.id, circle.a, 0.0001, 0.0, 0.0, 2);
		}
		else{
			fprintf(file, "circularregion %d radius %f vector %f %f %f form %d\n", circle.id, circle.a, circle.normal(0)*circle.g, circle.normal(1)*circle.g, circle.normal(2)*circle.g, circle.form);
			//fprintf(file, "intersector %d radius %f vector %f %f %f\n", i, circle.sphereRadius, circle.vector(0), circle.vector(1), circle.vector(2));
		}
	} 
	
	
	for(int i=0;i<sasas.size();++i){
		sasa = sasas[i];
		fprintf(file, "sasa %d size %d\n", i, sasa.sasa.size());
		k=0;
		for(int j=0; j<sasa.sasa.size(); ++j){
			fprintf(file, "intersectionpoint %d vector %f %f %f\n", j, sasa.sasa[j].vector(0), sasa.sasa[j].vector(1), sasa.sasa[j].vector(2));
		}
		
	}
	
	fclose(file);
	
	*/	
	

}


void Tessellation::outputTessellation(string filename){
	SASAs::iterator it_atom;
	SASASegmentList::iterator it;
	FILE* file;
	
	file = fopen (filename.c_str(),"w");
	
	int i;
	i=0;
	for(it_atom=sasasForMolecule.begin(); it_atom!=sasasForMolecule.end(); ++it_atom){
		for(it=it_atom->begin(); it!=it_atom->end(); ++it){
			fprintf(file, "%f %f %f %f %f %f %f %f %d %d %d\n",
				it->normalForCircularInterface(0),it->normalForCircularInterface(1),it->normalForCircularInterface(2),
				it->rotation0.rotation, it->rotation1.rotation,
				it->psi.rotation,
				it->lambda.rotation,
				it->radius,
				it->hemisphere,
				it->form1,
				i);
		}
		++i;
							
	}
	
	fclose(file);
}




float Tessellation::calculateKappa(TessellationAxis tessellationAxis, Vector &v, CircularInterfaceForm form){
	float kappa;
	Vector n2,n;
	float s,dot_n_aux;
	Vector ex(3),ey(3),ez(3);
	
	if(form==SPLITTER) return 0;
	
	n2=cross(v, tessellationAxis.v);
	n=normalise(cross(tessellationAxis.v,n2));
	s = sgn(dot(n,tessellationAxis.planeNormal));
	dot_n_aux=dot(n,tessellationAxis.auxiliary);
	kappa = s*acos(dot_n_aux);
	//if(s<0) kappa=2*M_PI-kappa;
	
	
	
	
// 	n2=cross(v, ex);
// 	n=normalise(cross(ex,n2));
// 	s = sgn(dot(n,ez));
// 	dot_n_aux=dot(n,ey);
// 	kappa = s*acos(dot_n_aux);
	
	return kappa;
}


fmat33 Tessellation::rotz(float theta){
	fmat33 m;
	m(0,0) = cos(theta);
	m(0,1) = -sin(theta);
	m(0,2) = 0;
	
	m(1,0) = sin(theta);
	m(1,1) = cos(theta);
	m(1,2) = 0;
	
	m(2,0) = 0;
	m(2,1) = 0;
	m(2,2) = 1;
	
	return m;
		
}



fmat33 Tessellation::roty(float theta){
	fmat33 m;
	m(0,0) = cos(theta);
	m(0,1) = 0;
	m(0,2) = sin(theta);
	
	m(1,0) = 0;
	m(1,1) = 1;
	m(1,2) = 0;
	
	m(2,0) = -sin(theta);
	m(2,1) = 0;
	m(2,2) = cos(theta);
	
	return m;
		
}



fmat33 Tessellation::rotx(float theta){
	fmat33 m;
	m(0,0) = 1;
	m(0,1) = 0;
	m(0,2) = 0;
	
	m(1,0) = 0;
	m(1,1) = cos(theta);
	m(1,2) = -sin(theta);
	
	m(2,0) = 0;
	m(2,1) = sin(theta);
	m(2,2) = cos(theta);
	
	return m;
		
}
/*
Vector Tessellation::PHI2V(TessellationAxis tessellationAxis, float PHI, float psi, float lambda, float kappa){
	fmat33 T;
	fmat33 r0, r1;
	Vector n(3);
	Vector n0(3);
	Vector n1(3);
	Vector n2(3);
	Vector v(3),v2(2);
	Vector ex(3);
	Vector p(3);
	float ux,uy,uz;
	float C,S,t;
	float g;
	
	ex(0) = 1;
	ex(1) = 0;
	ex(2) = 0;

	r0 = rotz(psi);
	r1 = rotz(psi-lambda);
	
	g = 1-cos(lambda);

	n = r0 * ex;
	n0=g*n;
	n1 = r1 * ex;
	
	v = n1-n0;
	
	ux=n(0);
	uy=n(1);
	uz=n(2);
	C=cos(PHI);
	S=sin(PHI);
	t=1-cos(PHI);
	
	
	T(0,0) = t*ux*ux + C;
	T(0,1) = t*ux*uy - S*uz;
	T(0,2) = t*ux*uz + S*uy;
	
	T(1,0) = t*ux*uy + S*uz;
	T(1,1) = t*uy*uy + C;
	T(1,2) = t*uy*uz - S*ux;
	
	T(2,0) = t*ux*uz - S*uy;
	T(2,1) = t*uy*uz + S*ux;
	T(2,2) = t*uz*uz + C;
	
	v2 = T * v;
	
	p= n0 + v2;
	
	p = roty(kappa) * p;
	
	//if it's on the backside, flip it
	if(tessellationAxis.hemisphere==BACKHEMISPHERE) p=-p;
	
	return p;


}
*/


Vector Tessellation::PHI2V(TessellationAxis tessellationAxis, float PHI, float psi, float lambda, float kappa, CircularInterfaceForm form){
	fmat33 T;
	fmat33 r0, r1;
	Vector n(3);
	Vector n0(3);
	Vector n1(3);
	Vector n2(3);
	Vector v(3),v2(2);
	Vector ex(3),ey(3);
	Vector p(3);
	Vector x(3);
	double s1;
	
	
	if(tessellationAxis.hemisphere==BACKHEMISPHERE) s1=-1;
	else s1=1;
	
	
	
	ex(0) = 1;
	ex(1) = 0;
	ex(2) = 0;
	
	
	
	ey(0) = 0;
	ey(1) = 1;
	ey(2) = 0;
	
// 	ex=tessellationAxis.v;
// 	ey=tessellationAxis.auxiliary;
	
	
/*	
	
	x=-ey*sin(lambda);
	
	x = rotx(-PHI)*x;
	
	v = ex*cos(lambda);
	x=x+v;
	x=rotz(psi) * x;
	x=rotx(kappa) * x;
	*/


	
	
	x=-ey*sin(lambda);
	
	x = rotx(-s1*PHI)*x;
	
	v = ex*cos(lambda);
	x=x+v;
	x=rotz(s1*psi) * x;
	x=rotx(s1*kappa) * x;
	
	


	return x;


}



void Tessellation::print(FILE* outputfile){
	SASASegment s;
	Vector ip(3);
	
	fprintf(outputfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","atom","id","i","i2","index0","index1","index2","weight0","weight1","ipx","ipy","ipz","hspx","hspy","hspz");
	for(unsigned int i=0;i<sasasForMolecule.size();++i){
		for(unsigned int j=0;j<sasasForMolecule[i].size();++j){
			{//if(sasasForMolecule[i][j].form1==CONVEX){
				s=sasasForMolecule[i][j];

				//we filter to avoid outputting splitter intersectionpoints
				//printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n",i,j,s.i,s.i2,s.index0, s.index1,s.index2,s.weight);
				if(s.weight>0){
					if(s.i<=-2){
						if(s.form1==CONVEX){
							ip=s.normalForCircularInterface;
						}
						else{
							ip=-s.normalForCircularInterface;
						}
					}
					else{
						ip = s.v0;
					}
						
					fprintf(outputfile,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i,j,s.i,s.i2,s.index0, s.index1,s.index2,s.weight,s.weight1,ip(0),ip(1),ip(2),s.v1(0),s.v1(1),s.v1(2));
				}
			}
		}
	}
}




Benchmark Tessellation::getBenchmark(){
	return benchmark;
}



