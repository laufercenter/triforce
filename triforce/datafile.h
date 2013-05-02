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




class DataFile{
public:
	DataFile();
	DataFile(string name);
	
	Data1D* digest1DBinaryTable();
	Data3D* digest3DBinaryTable();
	Topology *digestMapCSV();
	Topology* digestTOP(TopologyMode topm);
	Molecule *digestGRO(Topology &top, bool useHydrogens);
	Molecule *digestPDB(Topology &top, bool useHydrogens);
	Molecule *digestXYZR();

		
private:
	
	string name;
	DataFileType type;
	

	double charArray2Double(char* data);
	int32_t charArray2FixedSignedInt32(char *data);
	double fixedSignedInt322Double(int32_t x, unsigned short fraction);
	int fixedSignedInt322Int(int32_t x);	
	
	
	
	double string2double(string s);
	int string2int(string s);
	vector<string>* split(string &s, char delimiter);
	string string2UpperCase(string s);
	string int2string(int d);

};
	



#endif //DATAFILE_H_