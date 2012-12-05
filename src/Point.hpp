
#ifndef POINT
#define POINT

#include <string>
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;
class Point {
    private:
        double* sample_data;
        double* sample_data_pearson_precomputed;
        int num_data_samples;




    public:
        Point(const Point& p);
        Point(const char* line);
        ~Point();
        
        bool operator==(const Point& other) const;

        std::string id;

        static double get_distance_between_points(const Point* p1, const Point* p2);
        static Point* get_centroid_of_points(const std::vector<Point*>& points);
        static void verify_proper_point_input_or_die(const std::vector<Point*>& points);
        static void filter_out_input_points(std::vector<Point*>& points);

        friend double* precompute_pearson_data(double* sample_data);
        friend double morten_pearsoncorr(int n, const vector<double>& v1, const vector<double>& v2);

        friend std::size_t hash_value(const Point &p);

        friend std::ostream& operator<<(std::ostream& ost, const Point& ls);

};



#endif
