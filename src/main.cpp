#include <stdio.h>
#include <string>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <fcntl.h>

#include <boost/program_options.hpp>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <boost/foreach.hpp>

#include <Point.hpp>
#include <CanopyClustering.hpp>
#include <Log.hpp>

using namespace std;
//using namespace boost::iostreams;
using namespace boost::program_options;


void verify_input_correctness(const options_description& command_line_desc, const variables_map& command_line_variable_map){
    if (command_line_variable_map.count("help")) {
        cout << command_line_desc << "\n";
        exit(1);
    }
         
    double max_canopy_dist = command_line_variable_map["max_canopy_dist"].as<double>();
    double max_close_dist = command_line_variable_map["max_close_dist"].as<double>();
    double max_merge_dist = command_line_variable_map["max_merge_dist"].as<double>();
    double max_step_dist = command_line_variable_map["max_step_dist"].as<double>();

    if( max_canopy_dist < 0 || max_canopy_dist > 1 || 
        max_close_dist < 0 || max_close_dist > 1 ||
        max_merge_dist < 0 || max_merge_dist > 1 ||
        max_step_dist < 0 || max_step_dist > 1 
      ){
        cout << "Distance values must be a number within range <0,1>" << endl;
        exit(1);
    }

    if (command_line_variable_map.count("point_input_file") != 1 ){
        cout << "Incorrect specification of point input file path!" << endl;
        exit(1);
    }
}

int main(int argc, char* argv[])
{
    //
    //Initialization
    //
    
    //Set logging level
    log_level = logINFO;

    //Prepare variables for command line input
    double max_canopy_dist;
    double max_close_dist;
    double max_merge_dist;
    double max_step_dist;
    string point_input_file;

    //Define and read command line options
    options_description command_line_desc("Allowed options");

    command_line_desc.add_options()
        ("help", "write help message")
        ("max_canopy_dist", value<double>(&max_canopy_dist)->default_value(0.1), "Max distance between a canopy center and a point in which the point belongs to the canopy")
        ("max_close_dist", value<double>(&max_close_dist)->default_value(0.4), "Max distance between a canopy center and a point in which the point will be considered close to the canopy")
        ("max_merge_dist", value<double>(&max_merge_dist)->default_value(0.03), "Max distance between two canopy centers in which the canopies should be merged")
        ("max_step_dist", value<double>(&max_step_dist)->default_value(0.1), "Max distance between canopy center and canopy centroid in which the centroid will be used as an origin for a new canpy")
        ("point_input_file", value<string>(&point_input_file), "Point input file");

    positional_options_description command_line_positional_desc;

    command_line_positional_desc.add("point_input_file",-1);

    variables_map command_line_variable_map;
    store(command_line_parser(argc,argv).options(command_line_desc).positional(command_line_positional_desc).run(), command_line_variable_map);
    notify(command_line_variable_map);

    //Verify command line input parameters
    verify_input_correctness(command_line_desc, command_line_variable_map);


    //
    //Parse point description file
    //
    vector<Point*> points;

    int point_file;
    char* point_file_mmap;
    struct stat statbuf;

    /* open the input file */
    point_file = open(point_input_file.c_str(), O_RDONLY);

    /* find size of input file */
    fstat(point_file,&statbuf);

    /* mmap the input file */
    point_file_mmap = (char*)mmap (0, statbuf.st_size, PROT_WRITE, MAP_PRIVATE, point_file, 0);

    _log(logINFO) << "File loaded into memory, generating points";


    char* line_start_ptr = point_file_mmap;
    char* line_end_ptr = point_file_mmap;
    char* mmap_end_ptr = point_file_mmap + statbuf.st_size;
    while(line_start_ptr < mmap_end_ptr){
        line_end_ptr = line_start_ptr;
        while(*line_end_ptr != '\n' && *line_end_ptr != '\r' && line_end_ptr < mmap_end_ptr){
            line_end_ptr++;
        }
        if(line_end_ptr != line_start_ptr && line_start_ptr != point_file_mmap){//Check if the line is not empty
            *line_end_ptr = '\0';
            //cout << line_start_ptr << endl;
            points.push_back(new Point(line_start_ptr));
        }
        line_start_ptr = ++line_end_ptr;
    }

    _log(logINFO) << "Points read, dropping file from memory";

    /* drop the file from memory*/
    munmap(point_file_mmap, statbuf.st_size);
    
    _log(logINFO) << "Running basic validation of points";

    verify_proper_point_input_or_die(points);

    _log(logINFO) << "Filtering points";

    filter_out_input_points(points);


    _log(logINFO) << "Finished input points processing";
    
    _log(logINFO) << "Number of points read: " << points.size();
    
    //
    //Run Canopy Clustering
    //
    std::vector<Canopy*> canopies;

    canopies = CanopyClusteringAlg::multi_core_run_clustering_on(points, max_canopy_dist, max_close_dist, max_merge_dist, max_step_dist);
    
    _log(logINFO) << "Finished clustering, number of canopies:" << canopies.size();

    CanopyClusteringAlg::filter_clusters_by_zero_medians(4, canopies);

    _log(logINFO) << "Finished filtering by medians, number of canopies:" << canopies.size();

    CanopyClusteringAlg::filter_clusters_by_single_point_skew(0.9, canopies);

    _log(logINFO) << "Finished filtering by single data point proportion, number of canopies:" << canopies.size();

    //Filter out canopies

    //cout << "####################Results####################" << endl;
    //cout << "####################Results####################" << endl;
    //cout << "####################Results####################" << endl;
    BOOST_FOREACH(Canopy* c, canopies){
        //cout << *c;
        cout << "canopy: ";
        BOOST_FOREACH(Point* p, c->neighbours){
            cout << p->id << ", ";
        }
        cout << endl;
    }

    BOOST_FOREACH(Canopy* c, canopies)
        delete c;

    BOOST_FOREACH(Point* point, points)
        delete point;




    return 0;
}
