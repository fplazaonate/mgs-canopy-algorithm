#include <stdio.h>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <boost/foreach.hpp>

#include <Point.hpp>
#include <CanopyClustering.hpp>
#include <Log.hpp>

using namespace std;
using namespace boost::iostreams;
using namespace boost::program_options;



void verify_input_correctness(const options_description& command_line_desc, const variables_map& command_line_variable_map){
    if (command_line_variable_map.count("help")) {
        cout << command_line_desc << "\n";
        exit(1);
    }

    if (command_line_variable_map.count("point_input_file") != 1 ){
        cout << "Incorrect specification of point input file path!" << endl;
        exit(1);
    }
}

int main(int argc, const char* argv[])
{
    //log_level = logDEBUG4;
    log_level = logINFO;
    //
    //Parse input options
    //
    options_description command_line_desc("Allowed options");

    command_line_desc.add_options()
        ("help", "write help message")
        ("point_input_file", value<string>(), "Point input file");

    positional_options_description command_line_positional_desc;

    command_line_positional_desc.add("point_input_file",-1);

    variables_map command_line_variable_map;
    store(command_line_parser(argc,argv).options(command_line_desc).positional(command_line_positional_desc).run(), command_line_variable_map);
    notify(command_line_variable_map);

    //Verify input correctness
    verify_input_correctness(command_line_desc, command_line_variable_map);


    //
    //Parse point description file
    //
    vector<Point> points;
    
    stream<mapped_file_source> point_input_file;
    point_input_file.open(command_line_variable_map["point_input_file"].as<string>());

    string line;
    getline(point_input_file, line);//Ignore first line
    while(getline(point_input_file, line)){
        points.push_back(Point(line.c_str()));
    }

    Point::verify_proper_point_input_or_die(points);

    _log(logINFO) << "Finished reading point input file";
    _log(logINFO) << "Number of points read: " << points.size();

    //cout << points[0] << endl;
    //cout << points[1] << endl;
    //cout << points[2] << endl;
    //cout << points[3] << endl;

    
    //
    //Run Canopy Clustering
    //
    std::vector<Canopy> canopies;
    canopies = CanopyClusteringAlg::single_core_run_clustering_on(points);

    //int i = 0;
    cout << "####################Results####################" << endl;
    cout << "####################Results####################" << endl;
    cout << "####################Results####################" << endl;
    BOOST_FOREACH(const Canopy& c, canopies){
        //cout << "Canopy: " << ++i << "\t\t";
        //BOOST_FOREACH(const Point& p, c.neighbours){
        //    cout << p.id << "\t";
        //}
        //cout << endl;
        cout << c;
    }



    return 0;
}
