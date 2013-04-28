/* Copyright 2013, Nils J. D. Drechsel & Jordi Villà-Freixa
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef MULTI_LAYERED_DEPTH_BUFFER_H_
#define MULTI_LAYERED_DEPTH_BUFFER_H_

#include <string>
#include <vector>
#include <set>
#include "boost/multi_array.hpp"

#include <armadillo>

#include "depth3d.h"
#include "data1d.h"


using namespace std;
using namespace arma;

#define THRESHOLD_KAPPA 0.5




enum LineType{
	FRONT=0,
	BACK=1,
};



typedef struct{
	double g;
	Vector v;
	bool flip;
	double psi;
}
LimitingInterface;


typedef struct{
	int i;
	double kappa;

}
DepthBufferProjection;



typedef multimap<double, LineType> DepthBufferLine;
typedef vector<DepthBufferLine>DepthBuffer;
typedef vector<ScanlineMode>DepthBufferMode;


class MultiLayeredDepthBuffer{
	
public:



	MultiLayeredDepthBuffer();
	MultiLayeredDepthBuffer(Depth3D &data, Data1D &occludedDistribution, Data1D &exposedDistribution, int m);
	void addSphere(Vector &v, double lambda, bool invert, double &kappa, double &psi, int index);
	bool passesBuffer(double kappa, double psi, double lambda, bool invert, int index, vector<Vector> &exposedVectors);
	bool getSplitterExposedVectors(vector<Vector> &exposedVectors);

	double probe(Vector &projection, bool &isExtended);
	Vector project(bool splitter, bool hemflip, Vector n, double kappa, double PHI, double psi, double lambda, bool invert, Vector &projection0, Vector &projection1);
	void extend();
	double g2;
	void print();
	void startNewCycle();
	void addProbe(double x);
	bool isCycleExposed();
	double exposedProbability();
	bool retrieveLimitingInterfaces(vector<LimitingInterface> &limitList);
	bool isExposed(double x);


	

private:
	int len;
	Depth3D data;
	DepthBuffer dbuffer;
	DepthBufferMode dmode;
	vector<int> extensionDistance;
	vector<double> extensionAngle;
	int mode;
	double exposed;
	double occluded;
	Data1D exposedDistribution;
	Data1D occludedDistribution;
	double prior_occluded;
	double prior_exposed;
	bool test;
	vector<double> gTable;
	vector<DepthBufferProjection> projections[2];
	vector<DepthBufferProjection> segment;
	vector<vector<DepthBufferProjection> > segments;
	
	
	Vector normalise(Vector x);
	int sgn(double d);
	Vector orientationPlaneNormal,  orientationAxisAuxiliary,  orientationAxis;
	
	
	ScanlineMode insertIntoLineBuffer(DepthBufferLine &line, double front, double back);
	bool wouldChangeLineBuffer(DepthBufferLine &line, double front, double back, bool &frontChanges, bool &backChanges);
	Vector convertToCartesian(DepthBufferProjection &pr);
	DepthBufferLine::iterator increaseLineInterator(DepthBufferLine::iterator it, DepthBufferLine &line);
	DepthBufferLine::iterator decreaseLineInterator(DepthBufferLine::iterator it, DepthBufferLine &line);
	
	void getDepthBounds(int p, double d, DepthBufferLine::iterator &it0, DepthBufferLine::iterator &it1);

	bool isWithinNumericalLimits(double x, double l);

	bool scanBuffer(double kappa, int index, vector<Vector> &exposedVectors, DepthInformation &dat);
	
	

};

#endif //MULTI_LAYERED_DEPTH_BUFFER_H_
