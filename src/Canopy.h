/*
 * Canopy.h
 *
 *  Created on: Sep 4, 2012
 *      Author: Hans-Henrik Stærfeldt
 */

#ifndef CANOPY_H_
#define CANOPY_H_
#include <list>
using namespace std;

template <class P> class Canopy : public P {
private:
	list <*P> points;
public:
	Canopy(P &p) :P(p) {};
	void CalculateCentroid();
	void AddPoint(P *p) {
		points.push_back(p);
	}
	void RemovePoint(P *p){
		// TODO: use set?
	}
	~Canopy() {};
};

#endif /* CANOPY_H_ */
