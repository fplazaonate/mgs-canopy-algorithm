
#ifndef POINT
#define POINT

#include <string>
#include <vector>
#include <iostream>

#include <ap.h>

class Point {
    private:
        alglib::real_1d_array sample_data;

    public:
        Point();
        Point(std::string point_file_line);
        
        bool operator==(const Point& other) const;

        std::string id;

        static double get_distance_between_points(const Point& p1, const Point& p2);
        static Point get_centroid_of_points(const std::vector<Point>& points);
        static void verify_proper_point_input_or_die(const std::vector<Point>& points);
        static std::vector<Point> filter_out_input_points(const std::vector<Point>& points);

        friend std::size_t hash_value(const Point &p);

        friend std::ostream& operator<<(std::ostream& ost, const Point& ls);

};


#endif
