
#include <iostream>
#include <sstream>
#include <algorithm>
#include <assert.h>

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>

#include <statistics.h>

#include <Point.hpp>

Point::Point(): id("-1"){}

Point::Point(std::string point_file_line){

    std::stringstream line_stream(point_file_line);

    //Read ID
    line_stream >> id;

    //Read Sample Values into temporary vector
    std::vector<double> sample_data_vector;
    double sample_value;
    while(line_stream >> sample_value){
        sample_data_vector.push_back(sample_value);
    }

    //Copy data from temp vector to array
    sample_data.setlength(sample_data_vector.size());
    for(int i = 0; i < sample_data_vector.size(); i++)
        sample_data[i] = sample_data_vector[i];

    
}

double Point::get_distance_between_points(const Point& p1, const Point& p2){
    return alglib::pearsoncorr2(p1.sample_data, p2.sample_data);
}

void Point::verify_proper_point_input_or_die(const std::vector<Point>& points){
    
    //Verify all points have the same number of samples
    int num_samples = points[0].sample_data.length();
    BOOST_FOREACH(const Point point, points)
        assert(point.sample_data.length() == num_samples);

    std::cout << "Observed number of samples per point: " << num_samples << std::endl;
        
    //TODO check if samples vary

}

Point Point::get_centroid_of_points(const std::vector<Point>& points){

    //TODO: median should be estimated using boost/accumulators/statistics/median.hpp

    Point centroid;
    
    int num_samples = points[0].sample_data.length();

    centroid.sample_data.setlength(num_samples);

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
        for(int i=0; i < p.sample_data.length(); i++){
            ost << p.sample_data[i] << "\t" ;
        }
        ost << std::endl;
        ost << "============================" << std::endl;
        
        return ost;
}
