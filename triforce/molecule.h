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

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <string>
#include <vector>
#include <map>

#include <armadillo>


using namespace std;
using namespace arma;

typedef vec Vector;

typedef struct
{
	double* x;
	double* y;
	double* z;
}
CoordinatesPointers; 


enum ForceField {
	Amber99SBildn,
	Amber99SB,
	Amber99,
	Amber96
};

typedef struct{
	double mass;
	double epsilon;
	double sigma;
}Parameters;



class Molecule{
	
public:
	Molecule(ForceField forcefield=Amber99SBildn);
	void addAtom(int i, double* x, double* y, double* z, string type);
	void update();
	vector<Vector> &coordinates();
	
private:
	
	double string2double(string s);
	vector<string> *split(string &s, char delimiter);
	string string2UpperCase(string s);
	

	
	
	vector<CoordinatesPointers> coordinatesPointers;
	vector<Vector> atoms;
	vector<double> sigmas;
	vector<double> epsilons;
	ForceField forcefield;
	map<string,Parameters> dict;
	
	
};

#endif //MOLECULE_H_
