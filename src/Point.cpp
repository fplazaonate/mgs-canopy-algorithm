
#include <iostream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>

#include <statistics.h>

#include <Point.hpp>
#include <Log.hpp>

Point::Point(): id("-1"){
}

//Point::Point(const Point& p){
//    num_data_samples = p.num_data_samples;
//    for(int i= 0; i < num_data_samples; i++){
//        sample_data.push_back(p.sample_data[i]);
//    }
//}

Point::Point(std::string point_file_line){

    std::vector<std::string> sample_data_vector;
    boost::split(sample_data_vector, point_file_line, boost::is_any_of(" \t"), boost::token_compress_on);


    //Read ID
    id = sample_data_vector[0];
    _log(logDEBUG2)<< "\"" << id << "\""; 

    //Copy data from temp vector to array
    //sample_data.setlength(sample_data_vector.size()-1);
    num_data_samples = sample_data_vector.size() - 1;

    //_log(logDEBUG2) << sample_data.length() << "\t" << sample_data_vector.size();
    for(int i = 1; i < sample_data_vector.size(); i++){
    //    sample_data[i-1] = (double)atof(sample_data_vector[i].c_str());
        sample_data.push_back((double)atof(sample_data_vector[i].c_str()));
        _log(logDEBUG3) << "\"" << sample_data_vector[i] << "\"" << "\t" << atof(sample_data_vector[i].c_str()) << "\t" << sample_data[i-1];
    }
}

Point::~Point(){
//    if(sample_data)
//        delete sample_data;
}

double morten_pearsoncorr(int n, const vector<double>& v1, const vector<double>& v2){
	int	i;
	double x0=0,y0=0;
	double t, nx, ny;
	double c;

    for(int i = 0; i < n; i++){
        x0+=v1[i];
        y0+=v2[i];
    }
    x0/=n;
    y0/=n;

	t = nx = ny = 0.0;

	for ( i=0;i<n;i++ ) {
		t += ( v1[i] - x0 ) * ( v2[i] - y0 );
		nx += ( v1[i] - x0 ) * ( v1[i] - x0 );
		ny += ( v2[i] - y0 ) * ( v2[i] - y0 );
	}

	if ( nx * ny == 0.0 )
		c = 0.0;
	else
		c = t/sqrt(nx*ny);

	return( c );
}

double Point::get_distance_between_points(const Point& p1, const Point& p2){

    int len = p1.num_data_samples;
    double dist = morten_pearsoncorr(len, p1.sample_data, p2.sample_data);
    //double dist = pearsoncorr2(p1.sample_data, p2.sample_data);

    if(log_level >= logDEBUG3){
        _log(logDEBUG3) << "<<<<<<DISTANCE<<<<<<";
        _log(logDEBUG3) << "point: " << p1.id;
        for(int i=0; i < p1.num_data_samples; i++){
            _log(logDEBUG3) << "\t"<<p1.sample_data[i];
        }
        _log(logDEBUG3) << "point: " << p2.id;
        for(int i=0; i < p2.num_data_samples; i++){
            _log(logDEBUG3) << "\t"<<p2.sample_data[i];
        }
        _log(logDEBUG3) << "distance: " << dist;
    }

    return dist; 
}

void Point::verify_proper_point_input_or_die(const std::vector<Point>& points){
    
    //Verify all points have the same number of samples
    //int num_samples = points[0].sample_data.length();
    //BOOST_FOREACH(const Point point, points){
    //    _log(logDEBUG) <<  point;
    //    assert(point.sample_data.length() == num_samples);
    //}
    int num_samples = points[0].num_data_samples;
    BOOST_FOREACH(const Point point, points){
        _log(logDEBUG) <<  point;
        assert(point.num_data_samples == num_samples);
    }

    _log(logINFO) << "Finished reading point input file";
    _log(logINFO) << "Observed number of samples per point: " << num_samples;
    _log(logINFO) << "Number of points read: " << points.size();
        
    //TODO check if samples vary

}

Point Point::get_centroid_of_points(const std::vector<Point>& points){

    //TODO: median should be estimated using boost/accumulators/statistics/median.hpp

    Point centroid(points[0]);
    
    assert(points.size());

    //std::cout << points[0] << std::endl;
    
    //int num_samples = points[0].sample_data.length();

    //centroid.sample_data.setlength(num_samples);

    //for(int sample = 0; sample < num_samples; sample++){
    int num_samples = points[0].num_data_samples;

    for(int sample = 0; sample < num_samples; sample++){

        std::vector<double> point_samples;

        BOOST_FOREACH(const Point& p, points){

            //TODO: this is slow as hell
            point_samples.push_back(p.sample_data[sample]);

        }

        std::sort(point_samples.begin(), point_samples.end());

        double median = point_samples[int(point_samples.size()/2)];

        centroid.sample_data[sample] = median;
    }
    
    return centroid;
}


std::size_t hash_value(const Point& p){
    boost::hash<std::string> hasher;
    return hasher(p.id);
}

bool Point::operator==(const Point& other) const {
    if(id == other.id){
        return true;
    } else {
        return false;
    }
}

std::ostream& operator<<(std::ostream& ost, const Point& p)
{
        ost << "============================" << std::endl;
        ost << "Point: " << p.id << std::endl;
        //for(int i=0; i < p.sample_data.length(); i++){
        for(int i=0; i < p.num_data_samples; i++){
            ost << p.sample_data[i] << "\t" ;
        }
        ost << std::endl;
        ost << "============================" << std::endl;
        
        return ost;
}

