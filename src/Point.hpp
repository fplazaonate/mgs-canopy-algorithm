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

        friend double* precompute_pearson_data(double* sample_data);
        friend std::size_t hash_value(const Point &p);
        friend std::ostream& operator<<(std::ostream& ost, const Point& ls);

        friend double get_distance_between_points(const Point* p1, const Point* p2);
        friend Point* get_centroid_of_points(const std::vector<Point*>& points);
        friend void verify_proper_point_input_or_die(const std::vector<Point*>& points);
        friend void filter_out_input_points(std::vector<Point*>& points, int min_non_zero_data_samples);
};



#endif
