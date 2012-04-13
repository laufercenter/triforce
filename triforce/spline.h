/* Copyright 2011, Nils J. D. Drechsel & Jordi Villà-Freixa
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

#ifndef SPLINE_H_
#define SPLINE_H_

typedef vec Vector;


class Spline{
	
public:
	Spline(vector<Vector> geometry);
	double f(double x);
	

	
private:
	int logSearch(double x);
	
	Matrix B;
	vector<double> X;
	vector<Vector> G;
	
};


#endif //SPLINE_H_

