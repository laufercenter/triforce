/* Copyright 2012, Nils J. D. Drechsel & Jordi Villà-Freixa
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

#ifndef DATAFILE_H_
#define DATAFILE_H_


#include <string>
#include <vector>
#include <map>
#include <armadillo>
#include "data3d.h"
#include "data1d.h"
#include "topology.h"
#include "molecule.h"


using namespace std;
using namespace arma;





enum DataFileType{
	MapCSV,
	Binary
};


enum DataFileMode{
	FloatDataFile,
	VectorialDataFile
};




class DataFile{
public:
	DataFile();
	DataFile(string name);
	
	Data1D* digest1DBinaryTable();
	Data3D<float>* digest3DBinaryTable();
	Data3D<Vector>* digest3DBinaryVectorialTable();
	
	Data3D<float>* digest6DBinaryTable();
	Topology *digestMapCSV();
	Topology* digestTOP(TopologyMode topm, Topology* t=NULL);
	Molecule *digestGRO(Topology &top, bool useHydrogens);
	Molecule *digestPDB(Topology &top, bool useHydrogens);
	Molecule *digestXYZR();
	Molecule *digestTRI();

		
		
	static float string2float(string s);
	static int string2int(string s);
	static vector<string>* split(string &s, char delimiter);
	static string string2UpperCase(string s);
	static string int2string(int d);
		
private:
	
	string name;
	DataFileType type;
	
	unsigned int parameter0Dim;
	unsigned int parameter1Dim;
	unsigned int parameter2Dim;
	unsigned int derivativeLevel, containsAuxiliaryData;
	unsigned int dataDim;
	Data3D<float> *tbl3Dfloat;
	Data3D<Vector> *tbl3Dvectorial;
	vector<unsigned int> dimensions;
	fstream *f;
	DataFileMode dfm;
	unsigned int numberDimensions;
	

	float charArray2Double(char* data);
	int32_t charArray2FixedSignedInt32(char *data);
	float fixedSignedInt322Double(int32_t x, unsigned short fraction);
	int fixedSignedInt322Int(int32_t x);	
	
	
	
	
	
	void readHeader3D();
	void readHeaderData3D();
	void readFloatData3D();
	void readVectorialData3D();
	void readGradients3D();
	void readAuxiliaryFloatData3D();
	void readAuxiliaryVectorialData3D();
	

};
	



#endif //DATAFILE_H_