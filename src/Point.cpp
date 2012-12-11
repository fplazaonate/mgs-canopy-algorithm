
#include <iostream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>

#include <Point.hpp>
#include <Log.hpp>
#include <Stats.hpp>

using namespace std;

Point::Point(const char* line){
    //Copy line to private buffer - strtok will modify it
    char* private_line = new char[strlen(line) + 1];
    strcpy(private_line,line);
    _log(logDEBUG2)<< "Point constructor, got: \"" << line << "\""; 

    //Read gene id - first word in the line
    char* word = strtok(private_line, "\t ");
    id = string(word);
    _log(logDEBUG2)<< "Point constructor, point id: \"" << id << "\""; 

    //Fill vector with data samples
    std::vector<double> sample_data_vector;
    sample_data_vector.reserve(700);

    word = strtok(NULL, "\t ");
    while( word != NULL ){
        sample_data_vector.push_back((double)atof(word));
        word = strtok(NULL, "\t ");
    }

    //Get number of samples for this point
    num_data_samples = sample_data_vector.size();
    _log(logDEBUG2)<< "Point constructor, num data samples: \"" << num_data_samples << "\""; 

    //Allocate and copy samples into array
    sample_data = new double[num_data_samples];
    sample_data_pearson_precomputed = new double[num_data_samples]; 
    for(int i = 0; i < sample_data_vector.size(); i++){
        sample_data[i] = sample_data_vector[i];
    }

    precompute_pearson_data(num_data_samples, sample_data, sample_data_pearson_precomputed);

    delete private_line;
}

Point::Point(const Point& p){
    id = p.id;
    num_data_samples = p.num_data_samples;

    sample_data = new double[num_data_samples];
    for(int i=0; i < num_data_samples;i++){
        sample_data[i] = p.sample_data[i];
    }

    sample_data_pearson_precomputed = new double[num_data_samples];
    for(int i=0; i < num_data_samples;i++){
        sample_data_pearson_precomputed[i] = p.sample_data_pearson_precomputed[i];
    }
}


Point::~Point(){
    delete sample_data;
    delete sample_data_pearson_precomputed;
}

bool Point::check_if_num_non_zero_samples_is_greater_than_x(int x){

    int num_non_zero_medians = 0;
    for(int i=0; i < num_data_samples; i++){
        if(sample_data[i] > 0.0000001){
            num_non_zero_medians++;
            if(num_non_zero_medians >= x)
                return true;
        }
    }
    return false;
}

bool Point::check_if_single_point_proportion_is_smaller_than(double x){
    double sum_data_samples = 0;
    double max_data_sample = 0;

    for(int i=0; i < num_data_samples; i++){
        if(max_data_sample < sample_data[i])
            max_data_sample = sample_data[i];

        sum_data_samples += sample_data[i];
    }

    return (max_data_sample / sum_data_samples) < x;
}

void verify_proper_point_input_or_die(const std::vector<Point*>& points){
    
    //Verify all points have the same number of samples
    int num_samples = points[0]->num_data_samples;
    BOOST_FOREACH(const Point* point, points){
        _log(logDEBUG) <<  *point;
        assert(point->num_data_samples == num_samples);
    }

    _log(logINFO) << "Finished reading point input file";
    _log(logINFO) << "Observed number of samples per point: " << num_samples;
    _log(logINFO) << "Number of points read: " << points.size();
        
    //TODO check if samples vary

}

double get_distance_between_points(const Point* p1, const Point* p2){

    int len = p1->num_data_samples;
    double dist = 1 - fabs(pearsoncorr_from_precomputed(len, p1->sample_data_pearson_precomputed, p2->sample_data_pearson_precomputed));

    //if(log_level >= logDEBUG3){
    //    _log(logDEBUG3) << "<<<<<<DISTANCE<<<<<<";
    //    _log(logDEBUG3) << "point: " << p1->id;
    //    for(int i=0; i < p1->num_data_samples; i++){
    //        _log(logDEBUG3) << "\t"<<p1->sample_data[i];
    //    }
    //    _log(logDEBUG3) << "point: " << p2->id;
    //    for(int i=0; i < p2->num_data_samples; i++){
    //        _log(logDEBUG3) << "\t"<<p2->sample_data[i];
    //    }
    //    _log(logDEBUG3) << "distance: " << dist;
    //}

    return dist; 
}

Point* get_centroid_of_points(const std::vector<Point*>& points){

    //TODO: median should be estimated using boost/accumulators/statistics/median.hpp

    Point* centroid = new Point(*points[0]);
    //TODO: Could be done better
    centroid->id = "!GENERATED!";
    
    assert(points.size());

    int num_samples = points[0]->num_data_samples;

    _log(logDEBUG4) << "num samples: " << num_samples;

    for(int i = 0; i < num_samples; i++){

        std::vector<double> point_samples;

        BOOST_FOREACH(const Point* p, points){

            //TODO: this is slow as hell
            point_samples.push_back(p->sample_data[i]);

        }

        std::sort(point_samples.begin(), point_samples.end());

        double median = -1;

        int mid = floor((point_samples.size() - 1)/2);
        if(!(point_samples.size()%2)){
            median = (point_samples[mid] + point_samples[mid+1])/2.0; 
        } else {
            median = point_samples[mid];
        }

        assert(median != -1);

        centroid->sample_data[i] = median;
    }

    precompute_pearson_data(centroid->num_data_samples, centroid->sample_data, centroid->sample_data_pearson_precomputed);
    
    return centroid;
}

void filter_out_input_points(std::vector<Point*>& points){

    //TODO: make a parameter
    int min_non_zero_data_samples = 3;

    int num_points = points.size();
    int num_data_samples = points[0]->num_data_samples;

    vector<int> indexes_to_drop;

    for(int i = 0; i < num_points; i++){
        int num_non_zero_samples = 0;
        for(int j = 0; j < num_data_samples; j++){
            if( points[i]->sample_data[j] > 0.0000001 ){
                num_non_zero_samples++;
                if(num_non_zero_samples >= min_non_zero_data_samples)
                    break;
            }
        }
        if(num_non_zero_samples < min_non_zero_data_samples)
            indexes_to_drop.push_back(i);
    }

    std::sort(indexes_to_drop.begin(), indexes_to_drop.end());
    std::reverse(indexes_to_drop.begin(), indexes_to_drop.end());

    for(int i = 0; i < indexes_to_drop.size(); i++){
        delete points[indexes_to_drop[i]];
        points.erase(points.begin() + indexes_to_drop[i]);
    }



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
        for(int i=0; i < p.num_data_samples; i++){
            ost << p.sample_data[i] << "\t" ;
        }
        ost << std::endl;
        ost << "============================" << std::endl;
        
        return ost;
}

