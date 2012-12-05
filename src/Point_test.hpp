#ifndef POINT_TEST
#define POINT_TEST

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include <Point.hpp>

using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE( test_centroid_calculation){

    vector<Point*> points;
    points.push_back(new Point("id_1 1 1 1 1 1 1 1 1"));
    points.push_back(new Point("id_2 2 2 2 2 2 2 2 2"));
    points.push_back(new Point("id_3 3 3 3 3 3 3 3 3"));
    points.push_back(new Point("id_4 4 4 4 4 4 4 4 4"));
    points.push_back(new Point("id_5 5 5 5 5 5 5 5 5"));
    points.push_back(new Point("id_6 6 6 6 6 6 6 6 6"));
    points.push_back(new Point("id_7 7 7 7 7 7 7 7 7"));
    points.push_back(new Point("id_8 8 8 8 8 8 8 8 8"));
    points.push_back(new Point("id_9 9 9 9 9 9 9 9 9"));
    points.push_back(new Point("id_10 10 10 10 10 10 10 10 10"));
    points.push_back(new Point("id_11 11 11 11 11 11 11 11 11"));

    Point::verify_proper_point_input_or_die(points);

    //The centroid should be a new point with only 6es
    for(int i = 0; i < points[0]->num_data_samples; i++){
        Point* p = Point::get_centroid_of_points(points);
        BOOST_CHECK_CLOSE(p->sample_data[i], 6, 0.001);
        delete p;
    }

    points.push_back(new Point("id_12 12 12 12 12 12 12 12 12"));

    //The centroid should be a new point with only 6.5es
    for(int i = 0; i < points[0]->num_data_samples; i++){
        Point* p = Point::get_centroid_of_points(points);
        BOOST_CHECK_CLOSE(p->sample_data[i],6.5,0.001);
        delete p;
    }

    BOOST_FOREACH(Point* p, points)
        delete p;

}


BOOST_AUTO_TEST_CASE( test_filtering_out_input_points ){

    vector<Point*> points;
    points.push_back(new Point("id_1 1 1 1 1 1 1 1 1"));
    points.push_back(new Point("id_2 0 0 0 0 0 2 2 2"));
    points.push_back(new Point("id_3 0 0 0 0 0 0 3 3"));
    points.push_back(new Point("id_4 4 4 4 4 4 4 4 4"));
    points.push_back(new Point("id_5 0 0 0 0 0 0 0 0"));
    points.push_back(new Point("id_6 0 0 0 0 0 0 0 6"));
    points.push_back(new Point("id_7 7 7 7 7 7 7 7 7"));
    points.push_back(new Point("id_8 0 0 0 0 0 0 0 0"));
    points.push_back(new Point("id_9 9 9 9 9 9 9 9 9"));
    points.push_back(new Point("id_10 0 0 0 0 10 10 10 10"));
    points.push_back(new Point("id_11 11 11 11 11 11 11 11 11"));

    Point::verify_proper_point_input_or_die(points);

    BOOST_CHECK_EQUAL(points.size(), 11);

    Point::filter_out_input_points(points);

    BOOST_CHECK_EQUAL(points.size(), 7);

    BOOST_CHECK(points[2] != NULL);
    BOOST_CHECK(points[4] != NULL);
    BOOST_CHECK(points[5] != NULL);
    BOOST_CHECK(points[7] != NULL);
    
    BOOST_FOREACH(Point* p, points)
        delete p;
}

#endif
