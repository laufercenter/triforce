#include "datafiledigest.h"

#include <sstream>
#include <iostream>


DataFileDigest::DataFileDigest(string n, DataFileType t){
	name = n;
	type = t;
}


double DataFileDigest::string2double(string s){
    istringstream strm;
    double v;
    double d;
    
    strm.str(s);
    strm >> d;
    return d;
}


vector<string>* DataFileDigest::split(string &s, char delimiter) {
    stringstream ss(s);
    string x;
    vector<string> *content=new vector<string>();
    
    while(getline(ss, x, delimiter)) {
        content->push_back(x);
    }
    return content;
}

string DataFileDigest::string2UpperCase(string s){
	string str=s;
	transform(str.begin(), str.end(),str.begin(), ::toupper);
	return str;
}



map<string,vector<double> >* DataFileDigest::digestParametersFile(){
	ifstream *ifs;
	vector<string> *content;
	string line;
	vector<double> v;
	map<string,vector<double> > *dict = new map<string,vector<double> >();
	string ident;
	int i;
	
	ifs = new ifstream(name.c_str(),ifstream::in);
	
	while(ifs->good()){
		std::getline(*ifs,line);
		content=split(line,' ');
		
		ident=string2UpperCase((*content)[0]);
		
		v.clear();
		for(i=1;i<content->size();i++){
			v.push_back(string2double((*content)[i]));
		}
		
		(*dict)[ident]=v;
	}
	
	return dict;
}


vector<Matrix> *DataFileDigest::digestSEAWaterFile(){
	ifstream *ifs;
	vector<string> *content;
	string line;
	vector<double> v;
	vector<Matrix> *B = new vector<Matrix>();
	Matrix M;
	bool sizeInformationGiven=false;
	bool digestingBlock=false;
	vector<double> tmpRow;
	int rows,cols, r, c;
	int i;
	
	vector<vector<double> > tmpBlock;
	
	ifs = new ifstream(name.c_str(),ifstream::in);
	
	while(ifs){
		std::getline(*ifs,line);
		
		if(line.length()>0){
			if(line[0]!='#'){
				content=split(line,',');
				if(content->size()==3 && !sizeInformationGiven) sizeInformationGiven=true;
				else{
					digestingBlock=true;
					v.clear();
					for(i=0; i<content->size(); ++i){
						v.push_back(string2double((*content)[i]));
					}
					tmpBlock.push_back(v);
				}
			}
		}
		else{
			if(digestingBlock){
				digestingBlock=false;
				rows=tmpBlock.size();
				cols=tmpBlock[0].size();
				
				M = Matrix(rows,cols);
				
				for(r=0;r<rows;++r)
					for(c=0;c<cols;++c){
						M(r,c) = tmpBlock[r][c];
					}
				B->push_back(M);
			}
		}
	}
	
	return B;
}

