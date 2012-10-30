//============================================================================
// Name        : CanopyClustering.cpp
// Author      : Hans Henrik Stærfeldt
// Version     :
// Copyright   : For CBS
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <MyPoint.hpp>
#include <Canopy.hpp>
#include <Clustering.hpp>

#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

    vector<MyPoint> points;
	
    points = read_point_file("./points", 1);

    //Let the clustering begin!
    {
        //Prepare basic structures
        boost::unordered_set marked_points;
        boost::unordered_set unmarked_points(points);

        while(unmarked_points.size() > 0){
            point = take_point_by_random(unmarked_points);

            Canopy canopy = Canopy::find_canopy_until(double delta, *unmarked_points, *marked_points);



        }

        print(Canopies);
    }



	return 0;
}
