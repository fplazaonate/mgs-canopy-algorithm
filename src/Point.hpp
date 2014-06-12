/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
 *
 * This file is part of Metagenomics Canopy Clustering Implementation.
 *
 * Metagenomics Canopy Clustering Implementation is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Metagenomics Canopy Clustering Implementation is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef POINT
#define POINT

#include <string>
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;
class Point {
    public:
        double* sample_data;
        double* sample_data_pearson_precomputed;
        int num_data_samples;

        Point(const Point& p);
        Point(const char* line);
        virtual ~Point();
        
        std::string id;

        bool check_if_num_non_zero_samples_is_greater_than_x(int x);
        bool check_if_single_point_proportion_is_smaller_than(double x);
        bool check_if_top_three_point_proportion_is_smaller_than(double x);

        friend double* precompute_pearson_data(double* sample_data);
        friend std::size_t hash_value(const Point &p);
        friend std::ostream& operator<<(std::ostream& ost, const Point& ls);

        friend double get_distance_between_points(const Point* p1, const Point* p2);
        friend Point* get_centroid_of_points(const std::vector<Point*>& points);
        friend void verify_proper_point_input_or_die(const std::vector<Point*>& points);
};



#endif
