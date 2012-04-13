#include "datafiledigest.h"



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



map<string,vector<double> > &DataFileDigest::digestParametersFile(){
	string line;
	vector<double> v;
	map<string,vector<double> > dict;
	
	ifs = new ifstream(name);
	
	while(ifs){
		std::getline(*ifs,line);
		content=split(line,' ');
		
		ident=string2UpperCase((*content)[0]);
		
		v.clear();
		for(i=1;i<content.size();i++){
			v.push_back((*content)[1]);
		}
		
		dict[ident]=v;
	}
	
	return dict;
}


vector<Matrix> &DataFileDigest::digestSEAWaterFile(){
	string line;
	vector<double> v;
	vector<Matrix> B;
	bool sizeInformationGiven=false;
	bool digestingBlock=false;
	vector<double> tmpRow;
	
	vector<vector<double> > tmpBlock;
	
	ifs = new ifstream(name);
	
	while(ifs){
		std::getline(*ifs,line);
		
		if(line.length>0){
			if(line[0]!="#"){
				content=split(line,',');
				if(content.size()==3 && !sizeInformationGiven) sizeInformationGiven=true;
				else{
					digestingBlock=true;
					tmpBlock.push_back(content);
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
						M(r)(c) = tmpBlock[r][c];
					}
				B.push_back(M);
			}
		}
	}
	
	return B;
}

