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

MyPoint centroid(vector<MyPoint> &l) {

    double mean_x=0, mean_y=0;

    BOOST_FOREACH(MyPoint p){
        mean_x += p.x;
        mean_y += p.y;
    }
    mean_x /= l.size();
    mean_y /= l.size();

    //TODO: this is possibly dangerous as a point like this might already exist

    //TODO: this causes it to return a copy
	return MyPoint(-1, mean_x, mean_y);
}
