/* Copyright 2012, Nils J. D. Drechsel & Jordi Villà-Freixa
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

#ifndef INTEGRATORSTATISTICAL_H_
#define INTEGRATORSTATISTICAL_H_

#include <string>
#include <vector>
#include <map>
#include "boost/multi_array.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "tessellation.h"
#include "molecule.h"
#include "integrator.h"

#include <armadillo>


using namespace std;
using namespace arma;

class IntegratorStatistical: public Integrator{
	
public:
	IntegratorStatistical();
	IntegratorStatistical(int trials);
	float integrate(Molecule *molecule, Tessellation *tessellation);
	
	
private:
	Tessellation* tessellation;
	Molecule* molecule;
	int trials;
	
	float angle(Vector &a, Vector &b);
	float integrateAtomicSASA(SASAsForAtom sasasForAtom);
	
};

#endif //INTEGRATORSTATISTICAL_H_
