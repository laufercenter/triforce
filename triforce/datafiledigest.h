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

#ifndef DATAFILEDIGEST_H_
#define DATAFILEDIGEST_H_


#include <string>
#include <vector>
#include <map>
#include <armadillo>


using namespace std;
using namespace arma;




typedef vec Vector;
typedef mat Matrix;


enum DataFileType{
	SEAWaterFile,
	ParametersFile
};


class DataFileDigest{
public:	
		
	DataFileDigest(string n="", DataFileType t=ParametersFile);
		
private:
	
	string name;
	DataFileType type;
		
	double string2double(string s);
	vector<string>* split(string &s, char delimiter);
	string string2UpperCase(string s);
	map<string,vector<double> >* digestParametersFile();
	vector<Matrix> *digestSEAWaterFile();

		
	
	
};
	



#endif //DATAFILEDIGEST_H_