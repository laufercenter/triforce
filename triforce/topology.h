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

#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

#include <string>
#include <vector>
#include <map>

#include <armadillo>


using namespace std;
using namespace arma;


typedef struct{
	double mass;
	double epsilon;
	double sigma;
}Parameters;



typedef map<string,Parameters > MapVector;
typedef map<string,string> MapString;






class Topology{
	
public:
	Topology();
	void setCell(string ident, Parameters p);
	Parameters getCell(string ident);
	Parameters getAssociatedCell(string ident);
	void setMassValue(string ident, double v);
	void setEpsilonValue(string ident, double v);
	void setSigmaValue(string ident, double v);
	void setAssociation(string ident0, string ident1);
	string getAssociation(string ident);
	bool contains(string ident);

	void print();
	
	
private:
	
	MapVector data;
	MapString associations;
	
	
	
};

#endif //TOPOLOGY_H_
