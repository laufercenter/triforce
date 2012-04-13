#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <iostream>
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
	Molecule(ForceField forcefield);
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
