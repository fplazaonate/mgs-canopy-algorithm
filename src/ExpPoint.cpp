/*
 * ExpPoint.cpp
 *
 *  Created on: Sep 4, 2012
 *      Author: l-hhs
 */

#include "ExpPoint.h"
#include <math.h>
#include <string>
using namespace std;

double distance(ExpPoint &a, ExpPoint &b) {
	// TODO: return pierson correlation between vectors a and b
	return 1.0;
}

ExpPoint centroid(list<ExpPoint> &l) {
	vector<double> c;
	// TODO: return for each position, centroid
	//  of vectors in value
	return ExpPoint("canopy",c);
}
