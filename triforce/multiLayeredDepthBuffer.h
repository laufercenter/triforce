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
	float g;
	Vector v;
	bool flip;
	float psi;
}
LimitingInterface;


typedef struct{
	int i;
	float kappa;

}
DepthBufferCoordinate;





typedef multimap<float, LineType> DepthBufferLine;
typedef vector<DepthBufferLine>DepthBuffer;
typedef vector<ScanlineMode>DepthBufferMode;



class MultiLayeredDepthBuffer{
	
public:



	MultiLayeredDepthBuffer();
	MultiLayeredDepthBuffer(Depth3D &data, Data1D &occludedDistribution, Data1D &exposedDistribution, int m);
	void addSphere(Vector &v, float lambda, bool invert, float &kappa, float &psi);
	bool passesBuffer(float kappa, float psi, float lambda, bool invert, vector<Vector> &exposedVectors);
	bool getSplitterExposedVectors(vector<Vector> &exposedVectors);
	void print();
	void startNewCycle();
	void addProbe(float x);
	bool isCycleExposed();
	float exposedProbability();
	bool isExposed(float x);
	float exposedArea();


	

private:
	unsigned int len;
	Depth3D data;
	DepthBuffer dbuffer;
	DepthBufferMode dmode;
	int mode;
	float exposed;
	float occluded;
	Data1D exposedDistribution;
	Data1D occludedDistribution;
	float prior_occluded;
	float prior_exposed;
	vector<float> gTable;
	float C;
	vector<DepthBufferCoordinate> projections[2];
	vector<DepthBufferCoordinate> segment;
	vector<vector<DepthBufferCoordinate> > segments;
	
	
	Vector orientationPlaneNormal,  orientationAxisAuxiliary,  orientationAxis;
	
	
	Vector normalise(Vector x);
	int sgn(float d);
	ScanlineMode insertIntoLineBuffer(DepthBufferLine &line, float front, float back);
	bool wouldChangeLineBuffer(DepthBufferLine &line, float front, float back, bool &frontChanges, bool &backChanges);
	Vector convertToCartesian(DepthBufferCoordinate &pr);
	DepthBufferLine::iterator increaseLineInterator(DepthBufferLine::iterator it, DepthBufferLine &line);
	DepthBufferLine::iterator decreaseLineInterator(DepthBufferLine::iterator it, DepthBufferLine &line);
	void getDepthBounds(int p, float d, DepthBufferLine::iterator &it0, DepthBufferLine::iterator &it1);
	bool isWithinNumericalLimits(float x, float l);
	bool scanBuffer(float kappa, vector<Vector> &exposedVectors, DepthInformation &dat);
	
	

};

#endif //MULTI_LAYERED_DEPTH_BUFFER_H_
