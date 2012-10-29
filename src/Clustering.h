/*
 * Clustering.h
 *
 *  Created on: Oct 24, 2012
 *      Author: l-hhs
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_
#include <list>
#include <string>
#include "Canopy.h"
using namespace std;

// TODO: find fitting data structures to maintain fast access to searched
//   such as 'unmarked' and other critical things. Becareful to avoid
//   cloning objects, and use pointers where needed! Rely on initial points
//   list being stable.

template<class P> class Clustering {
	double t1;
	double t2;
	double t3;
	int     citer;
	int     minsize;
public:
	Clustering(double t1, double t2, double t3, int citer, int minsize) {
		// Options: distance to use for merging points into canopies
		this->t1=t1;

		// Options: distance to use for merging canopies into canopies
		this->t2=t2;

		// Options: distance to require for centroid to move (ends iterations)
		this->t3=t3;

		// Options: do not move centroid (above=inf)
		this->citer=citer;

		// Options: canopy minimal size
		this->minsize=minsize;
	}

	list<Canopy<P>> cluster(list<P> &l) {
		// TODO: perform cluster

		const int maxthreads=12;
		list<Canopy<P>> result[maxthreads];
		list<Canopy<P>> finalresult;
		list<Canopy<P>> searchresult;

		// Using OpenMP
		for (int i=0; i<maxthreads; i++) {
			// Make worker cluster 1/N of the points into clusters.
			// write results in locally initiated lists 'result[threadno]'
		}

		// Merge lists of canopies in result[*]
		// Loop through canopies, foreach canopy:
		//   Clear searchresult;
		//   Using OpenMP locate hits in each of 'result[threadno]'
		//     add results to searchresult using OpenMP semaphores!
		//     delete hits from local 'result'
		//   Merge results into one Canopy after search.
		//   Add new canopy to finalresult if big enough

		return finalresult;
	}

	void clusterworker() {
//     Iterate while free points
//
//	    The worker picks a free point P.
//
//	    The point is to be center of a temporary canopy Ct.
//	    Iterate, generating temporary canopies Ct
//	      The canopy Ct is compared with all assigned and free points in this thread
//	        to make a temporary list of points P1..i.
//	      The temporary list of points is used to make a centroid for a new canopy Ct+1.
//		  If distance between centroids of Ct and Ct+1<T1,
//		    Register the new canopy Cn=Ct+1, and assign the list of points to it.
//		  P is marked 'tested' (it might fall outside the resulting canopy).
//	    End iteration
//
//		Make single-point canopies of remaining points.

	}

	list<*Canopy<P>> findcanopies(double r,Canopy<P> c,list<Canopy<P>> &l){
		list<*Canopy<P>> result;
		// Locate in parallel
		return list;
	}

	virtual ~Clustering() {
	}
	;
};

#endif /* CLUSTERING_H_ */
