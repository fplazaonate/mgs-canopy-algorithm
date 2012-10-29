/*
 * ExpPoint.h
 *
 *  Created on: Sep 4, 2012
 *      Author: l-hhs
 */

#ifndef EXPPOINT_H_
#define EXPPOINT_H_

#include "Point.h"
#include "Canopy.h"
#include <string>
#include <vector>
#include <iostream>
using namespace std;

class ExpPoint: public Point {
private:
	string name;
	vector<double> v;
public:
	ExpPoint(string name, vector<double> &v) {
		// TODO Assign vector
	}
	virtual ~ExpPoint() {};
	ostream& operator<<(std::ostream& ost, const ExpPoint& o)
	{
	        ost<<o.v;
	        return ost;
	}
	istream& operator>>(istream &ist,ExpPoint &o)
	{
	    ist>>o.v;
	    return ist;
	}
};

double distance(const ExpPoint &a, const ExpPoint &b);
ExpPoint centroid(list<ExpPoint> &l);

#endif /* EXPPOINT_H_ */
