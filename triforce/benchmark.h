/* Copyright 2013, Nils J. D. Drechsel & Jordi Villà-Freixa
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

#ifndef BENCHMARK_H_
#define BENCHMARK_H_



#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <stdio.h>
#include <string>


using namespace std;


typedef map<string,float> EntityList;


class Benchmark{
	
public:
	Benchmark();
	Benchmark(string section);
	void start(string phase);
	void stop();
	void addQuantity(string quantity, float x);
	void print(FILE* outputfile);

private:
	
	string section;
	string phase;
	clock_t clock_start;
	EntityList times;
	EntityList stats;
	
	
	clock_t ms();
	
	
	
};

#endif //BENCHMARK_H_
