/*
 * MyPoint.cpp
 *
 *  Created on: Sep 4, 2012
 *      Author: l-hhs
 */

#include <MyPoint.hpp>
#include <math.h>
#include <string>
using namespace std;

double distance(MyPoint &a, MyPoint &b) {
	return sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) );
}

MyPoint centroid(list<MyPoint> &l) {
	MyPoint c;
	return c;
}
