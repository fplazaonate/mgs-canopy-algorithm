/*
 * Point.h
 *
 *  Created on: Sep 4, 2012
 *      Author: Hans-Henrik Stærfeldt
 *
 *  This is the parent class for all data point implementations
 */

#ifndef POINT_H_
#define POINT_H_
#include <string>
using namespace std;

class Point {
public:
	Point();
	virtual ~Point();
};

double distance(const Point &a, const Point &b);
Point centroid(list<Point> &l);

/*
 * 	Point Centroid(list <Point> points);
 */

#endif /* POINT_H_ */
