/**
	Copyright (c) 2009-2014 Nils J. D. Drechsel, Christopher J. Fennell, Ken A. Dill, Jordi Villà-Freixa
	Email: nils.drechsel@gmail.com
	License: MIT-license, which can be found in file LICENSE.txt as well as here: http://opensource.org/licenses/MIT
*/

#include "datafile.h"

#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>
#include "topology.h"

#define INT32BYTEMASK 255
#define BINARY_DATA_BLOCK_SIZE_DOUBLE 8
#define BINARY_DATA_BLOCK_SIZE_FLOAT 4


using namespace boost;


DataFile::DataFile(){
}

DataFile::DataFile(string name){
	this->name=name;
}





Data1D* DataFile::digest1DBinaryTable(){
	unsigned char buffer0[4];
	unsigned char buffer1[1];
	unsigned char buffer2[2];
	char *buffer;
	unsigned int numberDimensions;
	vector<float> tmp;
	Data1D *tbl;
	vector<unsigned int> dimensions;
	unsigned int totalCells;
	unsigned int parameter0Dim;
	bool containsAuxiliaryData;
	
	//printf("digesting...\n");
	
	fstream f(name.c_str(),ios::binary|ios::in);
	
	//first 4 bytes should be "NILS" : 78, 73, 76, 83
	f.read((char*)buffer0,4);
	if(buffer0[0] != 78 || buffer0[1] != 73 || buffer0[2] != 76 || buffer0[3] != 83){
		printf("invalid data file\n");
		exit(-1);
	}

	//level of derivatives
	f.read((char*)buffer1,1);
	//derivativeLevel = static_cast<unsigned int>(buffer1[0]);

	//contains auxiliary values
	f.read((char*)buffer1,1);
	containsAuxiliaryData = (static_cast<unsigned int>(buffer1[0])==1);
	
	//this byte gives number of dimensions (should be 1)
	f.read((char*)buffer1,1);
	numberDimensions = static_cast<unsigned int>(buffer1[0]);
	if(numberDimensions!=1) exit(-1);
	
	//read headers
	f.read((char*)buffer2,2);
	//zeroth byte is how many dimensions the header has. should be 1
	//first byte is number of header entries
	parameter0Dim = (unsigned int)buffer2[1];
	
	
	tbl = new Data1D(parameter0Dim,containsAuxiliaryData);
	
	//read header
	buffer=new char[parameter0Dim*BINARY_DATA_BLOCK_SIZE_DOUBLE];
	f.read(buffer,parameter0Dim*BINARY_DATA_BLOCK_SIZE_DOUBLE);
	for(unsigned int i=0; i<parameter0Dim; i++){
		float v = charArray2Double(buffer+i*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		tbl->setHeaderParameter0Cell(i,v);
	}
	delete buffer;
	
	//read data
	totalCells = parameter0Dim;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE];
	f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE);
	for(unsigned int x=0; x<parameter0Dim; x++){
		float v = charArray2Double(buffer+x*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		tbl->setDataCell(x,v);
	}
	
	delete buffer;
	
	
	if(containsAuxiliaryData){
		//read data
		totalCells = parameter0Dim;
		buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE];
		f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		for(unsigned int x=0; x<parameter0Dim; x++){
			float v = charArray2Double(buffer+x*BINARY_DATA_BLOCK_SIZE_DOUBLE);
			tbl->setAuxiliaryCell(x,v);

		}
		
		delete buffer;
	}

	
	tbl->init();
	

	return tbl;
}	


void DataFile::readHeader3D(){
	unsigned char buffer0[4];
	unsigned char buffer1[1];
	unsigned char buffer2[9];
	
	//printf("digesting...\n");
	
	f = new fstream(name.c_str(),ios::binary|ios::in);
	
	//first 4 bytes should be "NILS" : 78, 73, 76, 83
	f->read((char*)buffer0,4);
	if(buffer0[0] != 78 || buffer0[1] != 73 || buffer0[2] != 76 || buffer0[3] != 83) exit(-1);

	//level of derivatives
	f->read((char*)buffer1,1);
	derivativeLevel = static_cast<unsigned int>(buffer1[0]);

	//contains phi values
	f->read((char*)buffer1,1);
	containsAuxiliaryData = static_cast<unsigned int>(buffer1[0]);
	
	//this byte gives number of dimensions (should be 3)
	f->read((char*)buffer1,1);
	numberDimensions = static_cast<unsigned int>(buffer1[0]);
	if(numberDimensions!=3) exit(-1);
	
	//read headers
	//header for PHI
	f->read((char*)buffer2,7);
	//zeroth byte is how many dimensions the PHI header has. should be 1
	//first byte is number of PHI entries
	parameter0Dim = static_cast<unsigned int>(buffer2[1]);
	//2nd byte is number of psi dimensions (should be 1)
	//3rd byte is number of psi entries
	//4th byte is number of lambda dimensions (should be 1)
	//5th byte is number of lambda entries
	parameter1Dim = static_cast<unsigned int>(buffer2[3]);
	parameter2Dim = static_cast<unsigned int>(buffer2[5]);
	//6th byte is number of data dimensions
	dataDim = static_cast<unsigned int>(buffer2[6]);	
	
}


void DataFile::readHeaderData3D(){
	char *buffer;
	
	//read PHI header
	buffer=new char[parameter0Dim*BINARY_DATA_BLOCK_SIZE_DOUBLE];
	f->read(buffer,parameter0Dim*BINARY_DATA_BLOCK_SIZE_DOUBLE);
	for(unsigned int i=0; i<parameter0Dim; i++){
		float v = charArray2Double(buffer+i*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		switch(dfm){
			case FloatDataFile: tbl3Dfloat->setHeaderCell(0,i,v);break;
			case VectorialDataFile: tbl3Dvectorial->setHeaderCell(0,i,v);break;
		}
	}
	delete buffer;
	//read psi header
	buffer=new char[parameter1Dim*BINARY_DATA_BLOCK_SIZE_DOUBLE];
	f->read(buffer,parameter1Dim*BINARY_DATA_BLOCK_SIZE_DOUBLE);
	for(unsigned int i=0; i<parameter1Dim; i++){
		float v = charArray2Double(buffer+i*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		switch(dfm){
			case FloatDataFile: tbl3Dfloat->setHeaderCell(1,i,v);break;
			case VectorialDataFile: tbl3Dvectorial->setHeaderCell(1,i,v);break;
		}
	}
	delete buffer;
	//read lambda header
	buffer=new char[parameter2Dim*BINARY_DATA_BLOCK_SIZE_DOUBLE];
	f->read(buffer,parameter2Dim*BINARY_DATA_BLOCK_SIZE_DOUBLE);
	for(unsigned int i=0; i<parameter2Dim; i++){
		float v = charArray2Double(buffer+i*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		switch(dfm){
			case FloatDataFile: tbl3Dfloat->setHeaderCell(2,i,v);break;
			case VectorialDataFile: tbl3Dvectorial->setHeaderCell(2,i,v);break;
		}
	}
	delete buffer;
}


void DataFile::readFloatData3D(){
	char *buffer;
	unsigned int totalCells;
	VectorInt p(3);
	
	//read data
	totalCells = parameter0Dim*parameter1Dim*parameter2Dim;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE];
	f->read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE);
	for(unsigned int z=0; z<parameter2Dim; z++){
		p(2)=z;
		for(unsigned int y=0; y<parameter1Dim; y++){
			p(1)=y;
			for(unsigned int x=0; x<parameter0Dim; x++){
				p(0)=x;
				float v = charArray2Double(buffer+((z*parameter1Dim+y)*parameter0Dim+x)*BINARY_DATA_BLOCK_SIZE_DOUBLE);
				tbl3Dfloat->setDataCell(p,v);
			}
		}
	}
	
	delete buffer;
}


void DataFile::readVectorialData3D(){
	char *buffer;
	unsigned int totalCells;
	Vector v;
	VectorInt p(3);
	v=Vector(dataDim);
	
	//read data
	totalCells = parameter0Dim*parameter1Dim*parameter2Dim*dataDim;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE];
	f->read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE);
	for(unsigned int z=0; z<parameter2Dim; z++){
		p(2)=z;
		for(unsigned int y=0; y<parameter1Dim; y++){
			p(1)=y;
			for(unsigned int x=0; x<parameter0Dim; x++){
				p(0)=x;
				for(unsigned int i=0; i<dataDim; i++){
					v(i) = charArray2Double(buffer+(((z*parameter1Dim+y)*parameter0Dim+x)*dataDim+i)*BINARY_DATA_BLOCK_SIZE_DOUBLE);
				}
				tbl3Dvectorial->setDataCell(p,v);
			}
		}
	}
	
	delete buffer;
}


void DataFile::readGradients3D(){
	char *buffer;
	unsigned int totalCells;
	VectorInt p(3);
	
	if(derivativeLevel>=2){
		//read gradients
		totalCells = parameter0Dim*parameter1Dim*parameter2Dim*3;
		buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE];
		f->read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		for(unsigned int z=0; z<parameter2Dim; z++){
			p(2)=z;
			for(unsigned int y=0; y<parameter1Dim; y++){
				p(1)=y;
				for(unsigned int x=0; x<parameter0Dim; x++){
					p(0)=x;
					for(unsigned int i=0; i<3; i++){
						float v = charArray2Double(buffer+(((z*parameter1Dim+y)*parameter0Dim+x)*3+i)*BINARY_DATA_BLOCK_SIZE_DOUBLE);
						switch(dfm){
							case FloatDataFile: tbl3Dfloat->setGradientCell(p,i,v);break;
							case VectorialDataFile: tbl3Dvectorial->setGradientCell(p,i,v);break;
						}
					}
				}
			}
		}
		
		delete buffer;
	}
	
	if(derivativeLevel>=3){
		//read hessians
		totalCells = parameter0Dim*parameter1Dim*parameter2Dim*3*3;
		buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE];
		f->read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		for(unsigned int z=0; z<parameter2Dim; z++){
			p(2)=z;
			for(unsigned int y=0; y<parameter1Dim; y++){
				p(1)=y;
				for(unsigned int x=0; x<parameter0Dim; x++){
					p(0)=x;
					for(unsigned int j=0; j<3; j++){
						for(unsigned int i=0; i<3; i++){
							float v = charArray2Double(buffer+((((z*parameter1Dim+y)*parameter0Dim+x)*3+j)*3+i)*BINARY_DATA_BLOCK_SIZE_DOUBLE);
							switch(dfm){
								case FloatDataFile: tbl3Dfloat->setHessianCell(p,j,i,v);break;
								case VectorialDataFile: tbl3Dvectorial->setHessianCell(p,j,i,v);break;
							}
						}
					}
				}
			}
		}
		
		delete buffer;
	}
}



void DataFile::readAuxiliaryFloatData3D(){
	char *buffer;
	unsigned int totalCells;
	VectorInt p(3);
	
	if(containsAuxiliaryData == 1){
		totalCells = parameter0Dim*parameter1Dim*parameter2Dim;
		buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE];
		f->read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		for(unsigned int z=0; z<parameter2Dim; z++){
			p(2)=z;
			for(unsigned int y=0; y<parameter1Dim; y++){
				p(1)=y;
				for(unsigned int x=0; x<parameter0Dim; x++){
					p(0)=x;
					float v = charArray2Double(buffer+((z*parameter1Dim+y)*parameter0Dim+x)*BINARY_DATA_BLOCK_SIZE_DOUBLE);
					tbl3Dfloat->setAuxiliaryCell(p,v);
				}
			}
		}
		
		delete buffer;
	}	
	
}



void DataFile::readAuxiliaryVectorialData3D(){
	char *buffer;
	unsigned int totalCells;
	Vector v;
	VectorInt p(3);
	v=Vector(dataDim);
	
	if(containsAuxiliaryData == 1){
		
		//read data
		totalCells = parameter0Dim*parameter1Dim*parameter2Dim*dataDim;
		buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE];
		f->read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_DOUBLE);
		for(unsigned int z=0; z<parameter2Dim; z++){
			p(2)=z;
			for(unsigned int y=0; y<parameter1Dim; y++){
				p(1)=y;
				for(unsigned int x=0; x<parameter0Dim; x++){
					p(0)=x;
					for(unsigned int i=0; i<dataDim; i++){
						v(i) = charArray2Double(buffer+(((z*parameter1Dim+y)*parameter0Dim+x)*dataDim+i)*BINARY_DATA_BLOCK_SIZE_DOUBLE);
					}
					tbl3Dvectorial->setAuxiliaryCell(p,v);
				}
			}
		}
		delete buffer;
}	}
	


Data3D<float>* DataFile::digest3DBinaryTable(){
	vector<float> tmp;
	vector<unsigned int> dimensions;
	VectorInt p(3);
	
	dfm=FloatDataFile;
	readHeader3D();
	p(0) = parameter0Dim;
	p(1) = parameter1Dim;
	p(2) = parameter2Dim;
	tbl3Dfloat = new Data3D<float>(p, derivativeLevel, containsAuxiliaryData);
	readHeaderData3D();
	readFloatData3D();
	readGradients3D();
	readAuxiliaryFloatData3D();
	
	tbl3Dfloat->init();
	

	return tbl3Dfloat;
}	


Data3D<Vector>* DataFile::digest3DBinaryVectorialTable(){
	vector<float> tmp;
	VectorInt p(3);
	dfm=VectorialDataFile;
	readHeader3D();
	p(0) = parameter0Dim;
	p(1) = parameter1Dim;
	p(2) = parameter2Dim;
	tbl3Dvectorial = new Data3D<Vector>(p, derivativeLevel, containsAuxiliaryData);
	readHeaderData3D();
	readVectorialData3D();
	readGradients3D();
	readAuxiliaryVectorialData3D();
	
	tbl3Dvectorial->init();
	

	return tbl3Dvectorial;
}	



/*

Data6D* DataFile::digest6DBinaryTable(){
	unsigned char buffer0[4];
	unsigned char buffer1[1];
	unsigned char buffer2[8];
	char *buffer;
	unsigned int numberDimensions;
	vector<float> tmp;
	Data3D *tbl;
	vector<unsigned int> dimensions;
	unsigned int totalCells;
	unsigned int parameter0Dim;
	unsigned int parameter1Dim;
	unsigned int parameter2Dim;
	unsigned int parameter3Dim;
	unsigned int parameter4Dim;
	unsigned int parameter5Dim;
	unsigned int derivativeLevel, containsAuxiliaryData;
	
	//printf("digesting...\n");
	
	fstream f(name.c_str(),ios::binary|ios::in);
	
	//first 4 bytes should be "NILS" : 78, 73, 76, 83
	f.read((char*)buffer0,4);
	if(buffer0[0] != 78 || buffer0[1] != 73 || buffer0[2] != 76 || buffer0[3] != 83) exit(-1);

	//level of derivatives (should be 0)
	f.read((char*)buffer1,1);
	derivativeLevel = static_cast<unsigned int>(buffer1[0]);

	//contains phi values (should be 0)
	f.read((char*)buffer1,1);
	containsAuxiliaryData = static_cast<unsigned int>(buffer1[0]);
	
	//this byte gives number of dimensions (should be 6)
	f.read((char*)buffer1,1);
	numberDimensions = static_cast<unsigned int>(buffer1[0]);
	if(numberDimensions!=6) exit(-1);
	
	//read headers
	//header for PHI
	f.read((char*)buffer2,12);
	//zeroth byte is how many dimensions the PHI header has. should be 1
	//first byte is number of PHI entries
	parameter0Dim = static_cast<unsigned int>(buffer2[1]);
	//2nd byte is number of psi dimensions (should be 1)
	//3rd byte is number of psi entries
	//4th byte is number of lambda dimensions (should be 1)
	//5th byte is number of lambda entries
	parameter1Dim = static_cast<unsigned int>(buffer2[3]);
	parameter2Dim = static_cast<unsigned int>(buffer2[5]);
	parameter3Dim = static_cast<unsigned int>(buffer2[7]);
	parameter4Dim = static_cast<unsigned int>(buffer2[9]);
	parameter5Dim = static_cast<unsigned int>(buffer2[11]);
	

	
	//printf("dimensions: %d, rows in header: %d, maxdim: %d\n",numberDimensions,nrowsHeader,maxdim); 
	
	tbl = new Data6D(parameter0Dim, parameter1Dim, parameter2Dim, parameter3Dim, parameter4Dim, parameter5Dim);
	
	//read PHI header
	buffer=new char[parameter0Dim*BINARY_DATA_BLOCK_SIZE_FLOAT];
	f.read(buffer,parameter0Dim*BINARY_DATA_BLOCK_SIZE_FLOAT);
	for(unsigned int i=0; i<parameter0Dim; i++){
		float v = charArray2Float(buffer+i*BINARY_DATA_BLOCK_SIZE_FLOAT);
		tbl->setHeaderParameter0Cell(i,v);
	}
	delete buffer;
	//read psi header
	buffer=new char[parameter1Dim*BINARY_DATA_BLOCK_SIZE_FLOAT];
	f.read(buffer,parameter1Dim*BINARY_DATA_BLOCK_SIZE_FLOAT);
	for(unsigned int i=0; i<parameter1Dim; i++){
		float v = charArray2Float(buffer+i*BINARY_DATA_BLOCK_SIZE_FLOAT);
		tbl->setHeaderParameter1Cell(i,v);
	}
	delete buffer;
	//read lambda header
	buffer=new char[parameter2Dim*BINARY_DATA_BLOCK_SIZE_FLOAT];
	f.read(buffer,parameter2Dim*BINARY_DATA_BLOCK_SIZE_FLOAT);
	for(unsigned int i=0; i<parameter2Dim; i++){
		float v = charArray2Float(buffer+i*BINARY_DATA_BLOCK_SIZE_FLOAT);
		tbl->setHeaderParameter2Cell(i,v);
	}
	delete buffer;
	//read d header
	buffer=new char[parameter3Dim*BINARY_DATA_BLOCK_SIZE_FLOAT];
	f.read(buffer,parameter3Dim*BINARY_DATA_BLOCK_SIZE_FLOAT);
	for(unsigned int i=0; i<parameter3Dim; i++){
		float v = charArray2Float(buffer+i*BINARY_DATA_BLOCK_SIZE_FLOAT);
		tbl->setHeaderParameter2Cell(i,v);
	}
	delete buffer;
	//read eps header
	buffer=new char[parameter4Dim*BINARY_DATA_BLOCK_SIZE_FLOAT];
	f.read(buffer,parameter4Dim*BINARY_DATA_BLOCK_SIZE_FLOAT);
	for(unsigned int i=0; i<parameter4Dim; i++){
		float v = charArray2Float(buffer+i*BINARY_DATA_BLOCK_SIZE_FLOAT);
		tbl->setHeaderParameter2Cell(i,v);
	}
	delete buffer;
	//read sig header
	buffer=new char[parameter5Dim*BINARY_DATA_BLOCK_SIZE_FLOAT];
	f.read(buffer,parameter5Dim*BINARY_DATA_BLOCK_SIZE_FLOAT);
	for(unsigned int i=0; i<parameter5Dim; i++){
		float v = charArray2Float(buffer+i*BINARY_DATA_BLOCK_SIZE_FLOAT);
		tbl->setHeaderParameter2Cell(i,v);
	}
	delete buffer;
	
	//read data
	totalCells = parameter0Dim*parameter1Dim*parameter2Dim*parameter3Dim*parameter4Dim*parameter5Dim;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE_FLOAT];
	f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE_FLOAT);
	for(unsigned int w=0; w<parameter5Dim; w++)
		for(unsigned int v=0; v<parameter4Dim; v++)
			for(unsigned int u=0; u<parameter3Dim; u++)
				for(unsigned int z=0; z<parameter2Dim; z++)
					for(unsigned int y=0; y<parameter1Dim; y++)
						for(unsigned int x=0; x<parameter0Dim; x++){
							float w0 = w*parameter4Dim*parameter3Dim*parameter2Dim*parameter1Dim*parameter0Dim;
							float v0 = v*parameter3Dim*parameter2Dim*parameter1Dim*parameter0Dim;
							float u0 = u*parameter2Dim*parameter1Dim*parameter0Dim;
							float z0 = z*parameter1Dim*parameter0Dim;
							float y0 = y*parameter0Dim;
							float x0 = x;
							float val = charArray2Float(buffer+(w0+v0+u0+z0+y0+x0)*BINARY_DATA_BLOCK_SIZE_FLOAT);
							tbl->setDataCell(x,y,z,u,v,w,val);
						}
	
	delete buffer;
	
	
	tbl->init();
	

	return tbl;
}	
*/


/*
void DataFile::float2charArray(float x, char* data){
	unsigned int exponent;
	float significand;
	int32_t significandInt32, exponentInt32;
	
	significand=frexp(x,&exponent);
	
	significandInt32 = float2FixedUnsignedInt32(significand, 31);
	fixedSignedInt322CharArray(significandInt32, data);
	
	exponentInt32=int2FixedSignedInt32(exponent);
	fixedSignedInt322CharArray(exponentInt32, data+4);
	
}
*/

float DataFile::charArray2Double(char* data){
	int exponent;
	float significand;
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
			return numeric_limits<float>::quiet_NaN();
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
int32_t DataFile::float2FixedSignedInt32(float x, short fraction){
	int32_t d;
	int32_t factor= 1 << fraction;
	d = static_cast<int32_t>(x*factor);
	return d;
	
	
}
*/

float DataFile::fixedSignedInt322Double(int32_t x, unsigned short fraction){
	float d;
	int32_t factor = 1 << fraction;

	d = (static_cast<float>(x))/factor;
	return d;

}

/*
int32_t DataFile::int2FixedSignedInt32(unsigned int x){
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










float DataFile::string2float(string s){
    istringstream strm;
    float d;
    
    strm.str(s);
    strm >> d;
    return d;
}


int DataFile::string2int(string s){
    istringstream strm;
    int d;
    
    strm.str(s);
    strm >> d;
    return d;
}


string DataFile::int2string(int d){
    stringstream strm;
    string s;
    
    strm << d;
    
    return strm.str();
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
	vector<float> v;
	string ident;
	unsigned int i;
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
				v.push_back(string2float((*content)[i]));
			}
			
			data->setCell(ident,v);
		}
		
	}
	
	return data;
	*/

	return new Topology();
}


Topology* DataFile::digestTOP(TopologyMode topm, Topology* t){
	ifstream *ifs;
	vector<string> *content;
	string line;
	vector<float> v;
	string ident,ident0,ident1;
	Topology *data;
	string block;
	Parameters p;
	string prefix=string("");
	unsigned int chain;
	string chainstr;
	int residuenum;
	int prev_residuenum;
	

	ifs = new ifstream(name.c_str(),ifstream::in);

	if(t!=NULL) data=t;
	else data = new Topology(topm);
	
	p.mass = 0;
	p.epsilon= 0;
	p.sigma = 0;
	
	chain=0;
	chainstr="0";
	prev_residuenum=std::numeric_limits<int>::min();
	
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
								//data->setEpsilonValue(ident,string2float((*content)[6]));
								//data->setSigmaValue(ident,string2float((*content)[5]));
							}
							else{
								p.mass=-1;
								p.epsilon=string2float((*content)[6]);
								p.sigma=string2float((*content)[5])*10;
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
							
							residuenum = string2int((*content)[2]);
							if(residuenum<prev_residuenum){
								chain++;
								chainstr=int2string(chain);
							}
							prev_residuenum=residuenum;
							
							if(topm==GROMACS)
								ident0=chainstr+" "+(*content)[2]+" "+string2UpperCase((*content)[3])+" "+string2UpperCase((*content)[4]);
							else if(topm==GROMACS_GENERIC || topm==GROMACS_GENERIC2)
								ident0=string2UpperCase((*content)[3])+" "+string2UpperCase((*content)[4]);
							else if(topm==GROMACS_ELEMENTAL)
								ident0=string2UpperCase((*content)[4]);
							//is it already in the map?
							if(data->contains(ident1)){
								data->setMassValue(ident1,string2float((*content)[7]));
							}
							/*
							else{
								p.mass=string2float((*content)[7]);
								p.epsilon=-1;
								p.sigma=-1;
								data->setCell(ident1,p);
							}
							*/
							data->setAssociation(ident0,ident1);
						}
					}
					
				}
			}
			

		}
		
	}
	
	ifs->close();
	
	
	return data;
}




Molecule *DataFile::digestGRO(Topology &top, bool useHydrogens){
	ifstream *ifs;
	string line;
	unsigned int numbAtoms;
	unsigned int i;
	string block;
	Molecule *mol;
	float x,y,z;
	string ident;
	bool disregardingHydrogens;
	string residue;
	string atom;
	string xstr,ystr,zstr;
	unsigned int chain;
	string chainstr;
	string residuenumstr;
	unsigned int residuenum;
	unsigned int prev_residuenum;
	
	disregardingHydrogens=false;
	
	mol = new Molecule(top);
	
	ifs = new ifstream(name.c_str(),ifstream::in);

	
	//first line frame header
	std::getline(*ifs,line);
	//second line, number of atoms
	std::getline(*ifs,line);
	boost::trim(line);
	numbAtoms = string2int(line);
	i=0;
	
	chain=0;
	chainstr=string("0");
	prev_residuenum=std::numeric_limits<unsigned int>::min();

		
	while(ifs->good() && i<numbAtoms){
		std::getline(*ifs,line);
		
		if(line.length()>=44){
			//is it not a comment?
			if(line[0]!=';'){
				residuenumstr = line.substr(0,5);
				trim(residuenumstr);
				residuenum=string2int(residuenumstr);
				
				if(residuenum<prev_residuenum){
					chain++;
					chainstr = int2string(chain);
				}
				
				prev_residuenum = residuenum;
				
				residue = line.substr(5,5);
				trim(residue);
				atom = string2UpperCase(line.substr(10,5));
				trim(atom);
				xstr = line.substr(20,8);
				trim(xstr);
				x=string2float(xstr)*10;
				ystr = line.substr(28,8);
				trim(ystr);
				y=string2float(ystr)*10;
				zstr = line.substr(36,8);
				trim(zstr);
				z=string2float(zstr)*10;
				
				
				
				ident=chainstr + " " + residuenumstr + " " + residue + " " + atom;
				name=chainstr + " " + residuenumstr + " " + residue + " " + atom;
				try{
					if(useHydrogens || atom[0]!='H')
						mol->addInternallyStoredAtom(x,y,z,name,ident);
				}
				catch(AssociationException const &e){
					if(atom[0]!='H')
						printf("disregarded heavy atom %s for which no parameters could be found\n",ident.c_str());
					else
						if(!disregardingHydrogens){
							printf("disregarded some hydrogens for which no parameters could be found\n");
							disregardingHydrogens=true;
						}
					
					//do nothing
				}
				i++;
				
			}
			

		}
		
	}
	
	ifs->close();
	
	mol->init();
	
	
	return mol;
}





Molecule *DataFile::digestPDB(Topology &top, bool useHydrogens){
	ifstream *ifs;
	string line;
	unsigned int i;
	string block;
	Molecule *mol;
	float x,y,z;
	string ident;
	bool disregardingHydrogens;
	string residue;
	string atom,atom2;
	string xstr,ystr,zstr;
	unsigned int chain;
	string chainstr;
	string residuenumstr;
	int residuenum;
	int prev_residuenum;
	
	disregardingHydrogens=false;
	
	mol = new Molecule(top);
	
	ifs = new ifstream(name.c_str(),ifstream::in);
	
	

	
	chain=0;
	chainstr=string("0");
	prev_residuenum=std::numeric_limits<int>::min();

		
	while(ifs->good()){
		std::getline(*ifs,line);
		
		if(line.length()>=4){
			//is it not a comment?
			if(line.compare(0,4,"ATOM") == 0){
				residuenumstr = line.substr(22,4);
				trim(residuenumstr);
				residuenum=string2int(residuenumstr);
				
				if(residuenum<prev_residuenum){
					chain++;
					chainstr = int2string(chain);
				}
				
				prev_residuenum = residuenum;
				
				residue = line.substr(17,4);
				trim(residue);
				if(top.mode==GROMACS || top.mode==GROMACS_GENERIC || top.mode==GROMACS_GENERIC2){
					atom = string2UpperCase(line.substr(11,6));
				}
				else if(top.mode==GROMACS_ELEMENTAL){
					atom = string2UpperCase(line.substr(76,2));
				}
				atom2 = string2UpperCase(line.substr(76,2));
				trim(atom2);
				
				trim(atom);
				xstr = line.substr(30,8);
				trim(xstr);
				x=string2float(xstr);
				ystr = line.substr(39,8);
				trim(ystr);
				y=string2float(ystr);
				zstr = line.substr(46,8);
				trim(zstr);
				z=string2float(zstr);
				
				
				
				name=chainstr + " " + residuenumstr + " " + residue + " " + atom;
				
				if(top.mode==GROMACS || top.mode==GROMACS_GENERIC || top.mode==GROMACS_GENERIC2){
					ident=residue+" "+atom;
				}
				else if(top.mode==GROMACS_ELEMENTAL){
					ident=atom;
				}
				
				try{
					if(useHydrogens || atom[0]!='H')
						mol->addInternallyStoredAtom(x,y,z,name,ident);
				}
				catch(AssociationException const &e){
					if(top.mode==GROMACS_GENERIC2){
						try{
							ident=atom2;
							mol->addInternallyStoredAtom(x,y,z,name,ident);
						}
						catch(AssociationException const &e){
							if(atom[0]!='H')
								printf("disregarded heavy atom %s for which no parameters could be found\n",name.c_str());
							else
								if(!disregardingHydrogens){
									printf("disregarded some hydrogens for which no parameters could be found\n");
									disregardingHydrogens=true;
								}
						}
						
					}
					else{
				
						
						if(atom[0]!='H')
							printf("disregarded heavy atom %s for which no parameters could be found\n",name.c_str());
						else
							if(!disregardingHydrogens){
								printf("disregarded some hydrogens for which no parameters could be found\n");
								disregardingHydrogens=true;
							}
					}
					
					//do nothing
				}
				i++;
				
			}
			

		}
		
	}
	
	ifs->close();
	
	mol->init();
	
	return mol;
}




Molecule *DataFile::digestXYZR(){
	ifstream *ifs;
	string line;
	unsigned int i;
	Molecule *mol;
	float x,y,z,r;
	string ident;
	vector<string> *content;
	string element;
	
	
	mol = new Molecule();
	
	ifs = new ifstream(name.c_str(),ifstream::in);
	
	

	i = 0;
	while(ifs->good()){
		std::getline(*ifs,line);
		content=split(line,' ');
		if(content->size()>=3){
			element = (*content)[0];
			x = string2float((*content)[1]);
			y = string2float((*content)[2]);
			z = string2float((*content)[3]);
			r = string2float((*content)[4]);
			//ident = int2string(i);
			ident=element;
			mol->addInternallyStoredAtom(x, y, z, r, ident);
			++i;
		}
	}
	
	ifs->close();
	
	mol->init();
	
	return mol;
}


Molecule *DataFile::digestTRI(){
	ifstream *ifs;
	string line;
	unsigned int i;
	Molecule *mol;
	float x,y,z,r;
	string ident;
	vector<string> *content;
	string element;
	
	
	mol = new Molecule();
	
	ifs = new ifstream(name.c_str(),ifstream::in);
	
	

	i = 0;
	if(ifs->good()) std::getline(*ifs,line); //ignore header
	while(ifs->good()){
		std::getline(*ifs,line);
		content=split(line,'\t');
		if(content->size()>=12){
			element = (*content)[0];
			x = string2float((*content)[9]);
			y = string2float((*content)[10]);
			z = string2float((*content)[11]);
			r = string2float((*content)[6]);
			ident = (*content)[1];
			mol->addInternallyStoredAtom(x, y, z, r, ident);
			++i;
		}
	}
	
	ifs->close();
	
	mol->init();
	
	return mol;
}

