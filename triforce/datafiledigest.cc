#include "datafiledigest.h"

#include <sstream>
#include <iostream>

#define INT32BYTEMASK 255

/*
 * 	int exponent;
	double significand;
	
	significand=frexp(x,&exponent);

 * 
 * */

/*
void DataFileDigest::double2charArray(double x, char* data){
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

double DataFileDigest::charArray2Double(char* data){
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
void DataFileDigest::fixedSignedInt322CharArray(int32_t x, char *data){
	data[0]=static_cast<char>(x>>24);
	data[1]=static_cast<char>((x>>16) & INT32BYTEMASK);
	data[2]=static_cast<char>((x>>8) & INT32BYTEMASK);
	data[3]=static_cast<char>(x & INT32BYTEMASK);
}
*/

int32_t DataFileDigest::charArray2FixedSignedInt32(char *data){
	int32_t x=0;
	x |= (static_cast<int32_t>(data[0])) << 24;
	x |= (static_cast<int32_t>(data[1])) << 16;
	x |= (static_cast<int32_t>(data[2])) << 8;
	x |= (static_cast<int32_t>(data[3])) << 0;
}	

/*
int32_t DataFileDigest::double2FixedSignedInt32(double x, unsigned short fraction){
	int32_t d;
	int32_t factor= 1 << fraction;
	d = static_cast<int32_t>(x*factor);
	return d;
	
	
}
*/

double DataFileDigest::fixedSignedInt322Double(int32_t x, unsigned short fraction){
	double d;
	int32_t factor = 1 << fraction;

	d = (static_cast<double>(x))/factor;
	return d;

}

/*
int32_t DataFileDigest::int2FixedSignedInt32(int x){
	int32_t d;
	d = static_cast<int32_t>(x);
	return d;
	
}
*/

int DataFileDigest::fixedSignedInt322Int(int32_t x){
	int d;
	d = static_cast<int>(x);
	return d;
}
