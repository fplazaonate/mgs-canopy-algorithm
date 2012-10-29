/*
 * MyPoint.h
 *
 *  Created on: Sep 4, 2012
 *      Author: l-hhs
 */

#ifndef MYPOINT_H_
#define MYPOINT_H_

#include "Point.h"
#include <string>
#include <iostream>

class MyPoint: public Point {
private:
	double x;
	double y;
public:
	MyPoint(double x, double y) {
		this->x=x;
		this->y=y;
	}
	virtual ~MyPoint() {};
	ostream& operator<<(std::ostream& ost, const MyPoint& o)
	{
	        ost<<o.x<<','<<o.y;
	        return ost;
	}
	istream& operator>>(istream &ist,MyPoint &o)
	{
	    ist>>o.x>>','>>o.y;
	    return ist;
	}
};

double distance(const MyPoint &a, const MyPoint &b);
MyPoint centroid(list<MyPoint> &l);

#endif /* MYPOINT_H_ */
