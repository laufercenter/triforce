#include "datafile.h"

#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include "boost/multi_array.hpp"
#include <boost/algorithm/string.hpp>
#include "topology.h"

#define INT32BYTEMASK 255
#define BINARY_DATA_BLOCK_SIZE 8


DataFile::DataFile(){
}

DataFile::DataFile(const char* name, DataFileType type){
	this->name=string(name);
	this->type=type;
}



Data3D* DataFile::digest3DBinaryTable(){
	char buffer0[4];
	char buffer1[1];
	char buffer2[8];
	char shortbuffer[2];
	char *buffer;
	int numberDimensions;
	int nrowsHeader;
	vector<double> tmp;
	Data3D *tbl;
	vector<int> dimensions;
	int totalCells;
	int maxdim;
	int d;
	int parameter0Dim;
	int parameter1Dim;
	int parameter2Dim;
	
	//printf("digesting...\n");
	
	fstream f(name.c_str(),ios::binary|ios::in);
	
	//first 4 bytes should be "NILS" : 78, 73, 76, 83
	f.read(buffer0,4);
	if(buffer0[0] != 78 || buffer0[1] != 73 || buffer0[2] != 76 || buffer0[3] != 83) exit(-1);
	
	//this byte gives number of dimensions (should be 3)
	f.read(buffer1,1);
	numberDimensions = static_cast<int>(buffer1[0]);
	
	//read headers
	//header for PHI
	f.read(buffer2,6);
	//zeroth byte is how many dimensions the PHI header has. should be 1
	//first byte is number of PHI entries
	parameter0Dim = static_cast<int>(buffer2[1]);
	//2nd byte is number of psi dimensions (should be 1)
	//3rd byte is number of psi entries
	//4th byte is number of lambda dimensions (should be 1)
	//5th byte is number of lambda entries
	parameter1Dim = static_cast<int>(buffer2[3]);
	parameter2Dim = static_cast<int>(buffer2[5]);
	

	
	//printf("dimensions: %d, rows in header: %d, maxdim: %d\n",numberDimensions,nrowsHeader,maxdim); 
	
	tbl = new Data3D(parameter0Dim, parameter1Dim, parameter2Dim);
	
	//read PHI header
	buffer=new char[parameter0Dim*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,parameter0Dim*BINARY_DATA_BLOCK_SIZE);
	for(int i=0; i<parameter0Dim; i++){
		double v = charArray2Double(buffer+i*BINARY_DATA_BLOCK_SIZE);
		tbl->setHeaderParameter0Cell(i,v);
	}
	delete buffer;
	//read psi header
	buffer=new char[parameter1Dim*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,parameter1Dim*BINARY_DATA_BLOCK_SIZE);
	for(int i=0; i<parameter1Dim; i++){
		double v = charArray2Double(buffer+i*BINARY_DATA_BLOCK_SIZE);
		tbl->setHeaderParameter1Cell(i,v);
	}
	delete buffer;
	//read lambda header
	buffer=new char[parameter2Dim*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,parameter2Dim*BINARY_DATA_BLOCK_SIZE);
	for(int i=0; i<parameter2Dim; i++){
		double v = charArray2Double(buffer+i*BINARY_DATA_BLOCK_SIZE);
		tbl->setHeaderParameter2Cell(i,v);
	}
	delete buffer;
	
	//read data
	totalCells = parameter0Dim*parameter1Dim*parameter2Dim;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE);
	for(int z=0; z<parameter2Dim; z++)
		for(int y=0; y<parameter1Dim; y++)
			for(int x=0; x<parameter0Dim; x++){
				double v = charArray2Double(buffer+((z*parameter1Dim+y)*parameter0Dim+x)*BINARY_DATA_BLOCK_SIZE);
				tbl->setDataCell(x,y,z,v);
			}
	
	delete buffer;
	
	
	
	//read gradients
	totalCells = parameter0Dim*parameter1Dim*parameter2Dim*3;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE);
	for(int z=0; z<parameter2Dim; z++)
		for(int y=0; y<parameter1Dim; y++)
			for(int x=0; x<parameter0Dim; x++)
				for(int i=0; i<3; i++){
					double v = charArray2Double(buffer+(((z*parameter1Dim+y)*parameter0Dim+x)*3+i)*BINARY_DATA_BLOCK_SIZE);
					tbl->setGradientCell(x,y,z,i,v);
				}
	
	delete buffer;
	

	//read hessians
	totalCells = parameter0Dim*parameter1Dim*parameter2Dim*3*3;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE);
	for(int z=0; z<parameter2Dim; z++)
		for(int y=0; y<parameter1Dim; y++)
			for(int x=0; x<parameter0Dim; x++)
				for(int j=0; j<3; j++)
					for(int i=0; i<3; i++){
						double v = charArray2Double(buffer+((((z*parameter1Dim+y)*parameter0Dim+x)*3+j)*3+i)*BINARY_DATA_BLOCK_SIZE);
						tbl->setHessianCell(x,y,z,j,i,v);
					}
	
	delete buffer;
	
	tbl->init();
	

	return tbl;
}	






/*
void DataFile::double2charArray(double x, char* data){
	int exponent;
	double significand;
	int32_t significandInt32, exponentInt32;
	
	significand=frexp(x,&exponent);
	
	significandInt32 = double2FixedUnsignedInt32(significand, 31);
	fixedSignedInt322CharArray(significandInt32, data);
	
	exponentInt32=int2FixedSignedInt32(exponent);
	fixedSignedInt322CharArray(exponentInt32, data+4);
	
}
*/

double DataFile::charArray2Double(char* data){
	int exponent;
	double significand;
	int32_t significandInt32, exponentInt32;
	int32_t t0,t1,t2,t3;
	int32_t t4,t5,t6,t7;
	
	t0 = static_cast<uint32_t>(data[0]) & 255;
	t1 = static_cast<uint32_t>(data[1]) & 255;
	t2 = static_cast<uint32_t>(data[2]) & 255;
	t3 = static_cast<uint32_t>(data[3]) & 255;
	t4 = static_cast<uint32_t>(data[4]) & 255;
	t5 = static_cast<uint32_t>(data[5]) & 255;
	t6 = static_cast<uint32_t>(data[6]) & 255;
	t7 = static_cast<uint32_t>(data[7]) & 255;
	
	
	//this is how we save nan in the tables
	if(	t0 == 255 &&
		t1 == 255 &&
		t2 == 255 &&
		t3 == 255 &&
		t4 == 255 &&
		t5 == 255 &&
		t6 == 255 &&
		t7 == 255){
			return numeric_limits<double>::quiet_NaN();
	}

	significandInt32 = charArray2FixedSignedInt32(data);
	
	
	significand = fixedSignedInt322Double(significandInt32, 30);
	
	exponentInt32 = charArray2FixedSignedInt32(data+4);
	exponent = fixedSignedInt322Int(exponentInt32);
	
	return ldexp(significand, exponent);
}

/*
void DataFile::fixedSignedInt322CharArray(int32_t x, char *data){
	data[0]=static_cast<char>(x>>24);
	data[1]=static_cast<char>((x>>16) & INT32BYTEMASK);
	data[2]=static_cast<char>((x>>8) & INT32BYTEMASK);
	data[3]=static_cast<char>(x & INT32BYTEMASK);
}
*/

int32_t DataFile::charArray2FixedSignedInt32(char *data){
	int32_t x=0;
	int32_t t0,t1,t2,t3;
	t0 = static_cast<uint32_t>(data[0]) & 255;
	t1 = static_cast<uint32_t>(data[1]) & 255;
	t2 = static_cast<uint32_t>(data[2]) & 255;
	t3 = static_cast<uint32_t>(data[3]) & 255;
	
	//printf("T: %u %u %u %u\n",t0,t1,t2,t3);
	
	
	x |= t0 << 24;
	x |= t1 << 16;
	x |= t2 << 8;
	x |= t3 << 0;
	
	return x;
	
}	

/*
int32_t DataFile::double2FixedSignedInt32(double x, unsigned short fraction){
	int32_t d;
	int32_t factor= 1 << fraction;
	d = static_cast<int32_t>(x*factor);
	return d;
	
	
}
*/

double DataFile::fixedSignedInt322Double(int32_t x, unsigned short fraction){
	double d;
	int32_t factor = 1 << fraction;

	d = (static_cast<double>(x))/factor;
	return d;

}

/*
int32_t DataFile::int2FixedSignedInt32(int x){
	int32_t d;
	d = static_cast<int32_t>(x);
	return d;
	
}
*/

int DataFile::fixedSignedInt322Int(int32_t x){
	int d;
	d = static_cast<int>(x);
	return d;
}










double DataFile::string2double(string s){
    istringstream strm;
    double v;
    double d;
    
    strm.str(s);
    strm >> d;
    return d;
}


vector<string>* DataFile::split(string &s, char delimiter) {
	stringstream ss(s);
	string x;
	vector<string> *content=new vector<string>();

	while(getline(ss, x, delimiter)) {
		boost::trim(x);
		if(x.size()>0)
		content->push_back(x);
	}
	return content;
}

string DataFile::string2UpperCase(string s){
	string str=s;
	transform(str.begin(), str.end(),str.begin(), ::toupper);
	return str;
}



Topology* DataFile::digestMapCSV(){
/*	ifstream *ifs;
	vector<string> *content;
	string line;
	vector<double> v;
	string ident;
	int i;
	Topology *data;
	
	ifs = new ifstream(name.c_str(),ifstream::in);

	data = new Topology();
	
	while(ifs->good()){
		std::getline(*ifs,line);
		content=split(line,' ');
		
		if(content->size()>0){
		
			//printf("SIZE: %d\n",content->size());
			//for(i=0;i<content->size();++i)
			//	printf("CONTENT[%d]: %s\n",i,((*content)[i]).c_str());
			
			ident=string2UpperCase((*content)[0]);
			
			v.clear();
			for(i=1;i<content->size();i++){
				v.push_back(string2double((*content)[i]));
			}
			
			data->setCell(ident,v);
		}
		
	}
	
	return data;
	*/
}


Topology* DataFile::digestTOP(){
	ifstream *ifs;
	vector<string> *content;
	string line;
	vector<double> v;
	string ident,ident0,ident1;
	int i;
	Topology *data;
	string block;
	Parameters p;
	
	ifs = new ifstream(name.c_str(),ifstream::in);

	data = new Topology();
	
	p.mass = 0;
	p.epsilon= 0;
	p.sigma = 0;
	
	while(ifs->good()){
		std::getline(*ifs,line);
		content=split(line,' ');
		
		if(content->size()>0){
			
			//are we in section atomtypes?
			if(content->size()>=3 && (*content)[0][0]=='[' && (*content)[2][0]==']')
				block=(*content)[1];
			else{
				//are we collecting atomtype data?
				if(block=="atomtypes"){
					ident=string2UpperCase((*content)[0]);
					//is it not a comment?
					if((*content)[0][0]!=';'){
						//right amount of parameters?
						if(content->size()>=7){
							//is it already in the map?
							if(data->contains(ident)){
								data->setEpsilonValue(ident,string2double((*content)[6]));
								data->setSigmaValue(ident,string2double((*content)[5]));
							}
							else{
								p.mass=-1;
								p.epsilon=string2double((*content)[6]);
								p.sigma=string2double((*content)[5])*10;
								printf("SIGMA: %s: %s: %f\n",ident.c_str(),(*content)[5].c_str(),p.sigma);
								data->setCell(ident,p);
							}
						}
					}
				}
				if(block=="atoms"){
					//is it not a comment?
					if(content[0][0]!=";"){
						//right amount of parameters?
						if(content->size()>=8){
							ident1=string2UpperCase((*content)[1]);
							ident0=string2UpperCase((*content)[4]);
							//is it already in the map?
							if(data->contains(ident)){
								data->setMassValue(ident1,string2double((*content)[7]));
							}
							else{
								p.mass=string2double((*content)[7]);
								p.epsilon=-1;
								p.sigma=-1;
								data->setCell(ident1,p);
							}
							data->setAssociation(ident0,ident1);
						}
					}
				}
			}
			

		}
		
	}
	
	return data;
}




Molecule *DataFile::digestGRO(Topology &top){
	ifstream *ifs;
	vector<string> *content;
	string line;
	int numbAtoms;
	int i;
	string block;
	Molecule *mol;
	double x,y,z;
	string ident;
	
	mol = new Molecule(top);
	
	ifs = new ifstream(name.c_str(),ifstream::in);

	
	//first line frame header
	std::getline(*ifs,line);
	//second line, number of atoms
	std::getline(*ifs,line);
	boost::trim(line);
	numbAtoms = string2double(line);
	i=0;
		
	while(ifs->good() && i<numbAtoms){
		std::getline(*ifs,line);
		
		content=split(line,' ');
		
		if(content->size()>0){
			//is it not a comment?
			if((*content)[0][0]!=';'){
				//right amount of parameters?
				if(content->size()>=6){
					ident=string2UpperCase((*content)[1]);
					x=string2double((*content)[3])*10;
					y=string2double((*content)[4])*10;
					z=string2double((*content)[5])*10;
					mol->addRealAtom(x,y,z,ident);
					i++;
				}
			}
			

		}
		
	}
	
	return mol;
}








