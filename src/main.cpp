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


void verify_input_correctness(const options_description& all_options_desc, const variables_map& command_line_variable_map){
    if (command_line_variable_map.count("help")) {
        cout << all_options_desc<< "\n";
        exit(1);
    }
         
    int num_threads = command_line_variable_map["num_threads"].as<int>();
    double max_canopy_dist = command_line_variable_map["max_canopy_dist"].as<double>();
    double max_close_dist = command_line_variable_map["max_close_dist"].as<double>();
    double max_merge_dist = command_line_variable_map["max_merge_dist"].as<double>();
    double max_step_dist = command_line_variable_map["max_step_dist"].as<double>();
    string verbosity = command_line_variable_map["verbosity"].as<string>();

    if( max_canopy_dist < 0 || max_canopy_dist > 1 || 
        max_close_dist < 0 || max_close_dist > 1 ||
        max_merge_dist < 0 || max_merge_dist > 1 ||
        max_step_dist < 0 || max_step_dist > 1 
      ){
        cout << "Distance values must be a number within range <0,1>" << endl;
        exit(1);
    }

    if(num_threads < 1){
        cout << "Number of threads is the number of logical cpus you can utilize for the job, as shown by /proc/cpuinfo, it must be at least one" << endl;

    }

    if(num_threads > 999){
        cout << "if you really have more than 999 threads available... well congratulations - you need to recompile this program" << endl; 
    }

    if(verbosity != "error" && verbosity != "progress" && verbosity != "warn" && verbosity != "info" && verbosity != "debug" && verbosity != "debug1"){
        cout << "Verbosity must be one of: error, progress, warn, info, debug, debug1!" << endl;
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
    
    //Set initial logging level
    log_level = logINFO;

    //Prepare variables for command line input
    int num_threads;
    double max_canopy_dist;
    double max_close_dist;
    double max_merge_dist;
    double max_step_dist;
    string verbosity_option;
    string point_input_file;
    int min_num_non_zero_medians;
    double max_single_data_point_proportion;
    double stop_proportion_of_points;
    int stop_num_single_point_clusters;

    //Define and read command line options
    options_description all_options_desc("Allowed options");
    options_description general_options_desc("General");
    options_description algorithm_param_options_desc("Algorithm Parameters");
    options_description filter_options_desc("Output filter parameters");
    options_description early_stop_options_desc("Early stopping");
    options_description misc_options_desc("Miscellaneous");


    general_options_desc.add_options()
        ("point_input_file", value<string>(&point_input_file), "Point input file")
        ("num_threads", value<int>(&num_threads)->default_value(4), "IMPORTANT! Number of cpu threads to use.")
        ("verbosity,v", value<string>(&verbosity_option)->default_value("progress"), "Control how much information should be printed to the scree. Available levels according to their verbosity: error, progress, warn, info, debug, debug1.");

    algorithm_param_options_desc.add_options()
        ("max_canopy_dist", value<double>(&max_canopy_dist)->default_value(0.1), "Max distance between a canopy center and a point in which the point belongs to the canopy")
        ("max_close_dist", value<double>(&max_close_dist)->default_value(0.4), "Max distance between a canopy center and a point in which the point will be considered close to the canopy")
        ("max_merge_dist", value<double>(&max_merge_dist)->default_value(0.03), "Max distance between two canopy centers in which the canopies should be merged")
        ("max_step_dist", value<double>(&max_step_dist)->default_value(0.1), "Max distance between canopy center and canopy centroid in which the centroid will be used as an origin for a new canpy");

    filter_options_desc.add_options()
        ("filter_zero_medians", value<int>(&min_num_non_zero_medians)->default_value(4), "Return only those canopies that have at least N non-zero medians. Setting it to 0 will disable the filter.")
        ("filter_single_point", value<double>(&max_single_data_point_proportion)->default_value(0.9), "Don't return canopies containing a single median which divided by sum of all its medians is greater than X. Setting it to 1 disables the filter.");

    early_stop_options_desc.add_options()
        ("stop_on_proportion_of_points_clustered", value<double>(&stop_proportion_of_points)->default_value(0.5), "Stop clustering when X*total_number_of_points were clustered")
        ("stop_on_num_single_point_clusters", value<int>(&stop_num_single_point_clusters)->default_value(1000), "Stop clustering when X consecutive clusters had only one point in them");

    misc_options_desc.add_options()
        ("help", "write help message");

    all_options_desc.add(general_options_desc).add(algorithm_param_options_desc).add(filter_options_desc).add(early_stop_options_desc).add(misc_options_desc);

    positional_options_description command_line_positional_desc;
    command_line_positional_desc.add("point_input_file",-1);

    variables_map command_line_variable_map;
    store(command_line_parser(argc,argv).options(all_options_desc).positional(command_line_positional_desc).run(), command_line_variable_map);
    notify(command_line_variable_map);

    //Verify command line input parameters
    verify_input_correctness(all_options_desc, command_line_variable_map);


    //
    //Set user chosen logging level
    //
    if(verbosity_option == "error"){
        log_level = logERR;
    }else if(verbosity_option == "progress"){
        log_level = logPROGRESS;
    }else if(verbosity_option == "warn"){
        log_level = logWARN;
    }else if(verbosity_option == "info"){
        log_level = logINFO;
    }else if(verbosity_option == "debug"){
        log_level = logDEBUG;
    }else if(verbosity_option == "debug1"){
        log_level = logDEBUG1;
    }


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

    canopies = CanopyClusteringAlg::multi_core_run_clustering_on(points, num_threads, max_canopy_dist, max_close_dist, max_merge_dist, max_step_dist, stop_proportion_of_points, stop_num_single_point_clusters);
    
    _log(logINFO) << "Finished clustering, number of canopies:" << canopies.size();

    //
    //Filter out canopies
    //
    if(min_num_non_zero_medians){
        CanopyClusteringAlg::filter_clusters_by_zero_medians(min_num_non_zero_medians, canopies);
        _log(logINFO) << "Finished filtering by medians, number of canopies:" << canopies.size();
    }


    if(max_single_data_point_proportion < 0.99999){
        CanopyClusteringAlg::filter_clusters_by_single_point_skew(max_single_data_point_proportion, canopies);
        _log(logINFO) << "Finished filtering by single data point proportion, number of canopies:" << canopies.size();
    }


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
