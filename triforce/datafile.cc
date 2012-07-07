#include "datafile.h"

#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include "boost/multi_array.hpp"
#include <boost/algorithm/string.hpp>

#define INT32BYTEMASK 255
#define BINARY_DATA_BLOCK_SIZE 8


DataFile::DataFile(){
}

DataFile::DataFile(const char* name, DataFileType type){
	this->name=string(name);
	this->type=type;
}



Data3D* DataFile::digestBinaryTable(){
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
	
	printf("digesting...\n");
	
	fstream f(name.c_str(),ios::binary|ios::in);
	
	//read first row
	f.read(shortbuffer,2);
	numberDimensions = static_cast<int>(shortbuffer[0]);
	
	if(numberDimensions != 3) throw exception();
	
	nrowsHeader = static_cast<int>(shortbuffer[1]);
	
	
	
	
	//read second row
	maxdim=0;
	buffer=new char[numberDimensions];
	f.read(buffer,numberDimensions);
	for(int i=0; i<numberDimensions; i++){
		d = static_cast<int>(buffer[i]);
		if(d>maxdim) maxdim=d;
		dimensions.push_back(d);
	}
	delete buffer;

	
	printf("dimensions: %d, rows in header: %d, maxdim: %d\n",numberDimensions,nrowsHeader,maxdim); 
	
	tbl = new Data3D(dimensions);
	
	//read header
	buffer=new char[nrowsHeader*maxdim*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,nrowsHeader*maxdim*BINARY_DATA_BLOCK_SIZE);
	for(int i=0; i<maxdim; i++){
		for(int j=0; j<nrowsHeader; j++){
			double v = charArray2Double((buffer+(i*nrowsHeader+j)*BINARY_DATA_BLOCK_SIZE));
			tbl->setHeaderCell(j,i,v);
		}
	}
	delete buffer;
	
	//read data
	totalCells = dimensions[0]*dimensions[1]*dimensions[2];
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE);
	for(int z=0; z<dimensions[2]; z++)
		for(int y=0; y<dimensions[1]; y++)
			for(int x=0; x<dimensions[0]; x++){
				double v = charArray2Double(buffer+((z*dimensions[1]+y)*dimensions[0]+x)*BINARY_DATA_BLOCK_SIZE);
				tbl->setDataCell(x,y,z,v);
			}
	
	delete buffer;
	
	
	
	//read gradients
	totalCells = dimensions[0]*dimensions[1]*dimensions[2]*3;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE);
	for(int z=0; z<dimensions[2]; z++)
		for(int y=0; y<dimensions[1]; y++)
			for(int x=0; x<dimensions[0]; x++)
				for(int i=0; i<3; i++){
					double v = charArray2Double(buffer+(((z*dimensions[1]+y)*dimensions[0]+x)*3+i)*BINARY_DATA_BLOCK_SIZE);
					tbl->setGradientCell(x,y,z,i,v);
				}
	
	delete buffer;
	

	//read hessians
	totalCells = dimensions[0]*dimensions[1]*dimensions[2]*3*3;
	buffer=new char[totalCells*BINARY_DATA_BLOCK_SIZE];
	f.read(buffer,totalCells*BINARY_DATA_BLOCK_SIZE);
	for(int z=0; z<dimensions[2]; z++)
		for(int y=0; y<dimensions[1]; y++)
			for(int x=0; x<dimensions[0]; x++)
				for(int j=0; j<3; j++)
					for(int i=0; i<3; i++){
						double v = charArray2Double(buffer+((((z*dimensions[1]+y)*dimensions[0]+x)*3+j)*3+i)*BINARY_DATA_BLOCK_SIZE);
						tbl->setHessianCell(x,y,z,j,i,v);
					}
	
	delete buffer;
	

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



DataMapCSV* DataFile::digestMapCSV(){
	ifstream *ifs;
	vector<string> *content;
	string line;
	vector<double> v;
	string ident;
	int i;
	DataMapCSV *data;
	
	ifs = new ifstream(name.c_str(),ifstream::in);

	data = new DataMapCSV();
	
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
}


DataMapCSV* DataFile::digestTOP(){
	ifstream *ifs;
	vector<string> *content;
	string line;
	vector<double> v;
	string ident,ident0,ident1;
	int i;
	DataMapCSV *data;
	string block;
	
	ifs = new ifstream(name.c_str(),ifstream::in);

	data = new DataMapCSV();
	
	v.push_back(0);
	v.push_back(0);
	v.push_back(0);
	
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
								data->setCellValue(ident,1,string2double((*content)[6]));
								data->setCellValue(ident,2,string2double((*content)[5]));
							}
							else{
								v[0]=-1;
								v[1]=string2double((*content)[6]);
								v[2]=string2double((*content)[5]);
								data->setCell(ident,v);
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
								data->setCellValue(ident1,0,string2double((*content)[7]));
							}
							else{
								v[0]=string2double((*content)[7]);
								v[1]=-1;
								v[2]=-1;
								data->setCell(ident1,v);
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




Molecule *DataFile::digestGRO(){
	ifstream *ifs;
	vector<string> *content;
	string line;
	DataMapCSV *data;
	int numbAtoms;
	int i;
	
	Molecule* mol = new Molecule();
	
	ifs = new ifstream(name.c_str(),ifstream::in);

	
	//first line frame header
	std::getline(*ifs,line);
	//second line, number of atoms
	std::getline(*ifs,line);
	boost::trim(line)
	numbAtoms = string2double(line);
	i=0;
		
	while(ifs->good() && i<numbAtoms){
		
		content=split(line,' ');
		
		if(content->size()>0){
			//is it not a comment?
			if((*content)[0][0]!=';'){
				//right amount of parameters?
				if(content->size()>=6){
					
				}
			}
			

		}
		
	}
	
	return data;
}









