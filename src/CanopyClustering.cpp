//============================================================================
// Name        : CanopyClustering.cpp
// Author      : Hans Henrik Stærfeldt
// Version     :
// Copyright   : For CBS
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "MyPoint.h"
#include "Canopy.h"
#include "Clustering.h"
using namespace std;

int main(int argc, char *argv[]) {

	// Read data
	list<ExpPoint> points;

	// TODO: may need work :)
	cin >> points;

	// processing
	double t1 = 1.0;
	double t2 = 1.0;
	double t3 = 1.0;
	int citer = 1;
	int minsize = 0;

	Clustering<ExpPoint> *reducer = new Clustering<ExpPoint>(t1,t2,t3,citer,minsize);

	// result;
	list<Canopy<ExpPoint>> result = reducer->cluster(points);

	// TODO: may need work :)
	cout << result;

	return 0;
}
