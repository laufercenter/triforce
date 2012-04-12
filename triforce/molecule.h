#if !defined(DEF_MOLECULE_H)
#define DEF_MOLECULE_H

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
DoubleCoordinates; 


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
	
private:
	
	void update();
	vector<Vector>* getCoordinates();
	double string2double(string s);
	vector<string> *split(string &s, char delimiter);
	string string2UpperCase(string s);
	

	
	
	vector<DoubleCoordinates> doubleCoordinatePointers;
	vector<Vector> coordinates;
	vector<double> sigmas;
	vector<double> epsilons;
	ForceField forcefield;
	map<string,Parameters> dict;
	
	
};

#endif //DEF_MOLECULE_H
