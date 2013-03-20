#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <fcntl.h>

#include <boost/program_options.hpp>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <boost/foreach.hpp>
#include <boost/assign/std/vector.hpp> 

#include <Point.hpp>
#include <CanopyClustering.hpp>
#include <Log.hpp>
#include <program_options_misc.hpp>

using namespace std;
using namespace boost::program_options;
using namespace boost::assign;


int main(int argc, char* argv[])
{
    //
    //Initialization
    //
    
    //Set initial logging level
    log_level = logINFO;

    //Preapre Time Profile
    TimeProfile time_profile;
    time_profile.start_timer("Total");

    //Prepare variables for command line input
    string point_input_file;
    string output_file_path;
    string output_centers_file_path;
    string output_cluster_prefix;
    int num_threads;
    double max_canopy_dist;
    double max_close_dist;
    double max_merge_dist;
    double min_step_dist;
    string verbosity_option;
    int min_non_zero_data_samples;
    int min_num_non_zero_medians;
    double max_single_data_point_proportion;
    double stop_proportion_of_points;
    int stop_num_single_point_clusters;
    string canopy_size_stats_fp;
    bool show_progress_bar;
    bool print_time_statistics;


    //Define and read command line options
    options_description all_options_desc("Allowed options");
    options_description general_options_desc("General");
    options_description algorithm_param_options_desc("Algorithm Parameters");
    options_description filter_in_options_desc("Input filter parameters");
    options_description filter_out_options_desc("Output filter parameters");
    options_description early_stop_options_desc("Early stopping");
    options_description misc_options_desc("Miscellaneous");


    general_options_desc.add_options()
        ("point_input_file,i", value<string>(&point_input_file), "Point input file")
        ("output_file_path,o", value<string>(&output_file_path), "Provide path to file to which clusters will be written")
        ("output_centers_file_path,c", value<string>(&output_centers_file_path), "Provide path to file to which cluster centers will be written")
        ("output_clusters_prefix,p", value<string>(&output_cluster_prefix)->default_value("MGU"), "Provide path to file to which cluster centers will be written")
        ("num_threads,n", value<int>(&num_threads)->default_value(4), "IMPORTANT! Number of cpu threads to use.")
        ("verbosity,v", value<string>(&verbosity_option)->default_value("info"), "Control how much information should be printed to the scree. Available levels according to their verbosity: error, progress, warn, info, debug, debug1.");

    algorithm_param_options_desc.add_options()
        ("max_canopy_dist", value<double>(&max_canopy_dist)->default_value(0.1), "Max distance between a canopy center and a point in which the point belongs to the canopy")
        ("max_close_dist", value<double>(&max_close_dist)->default_value(0.4), "Max distance between a canopy center and a point in which the point will be considered close to the canopy")
        ("max_merge_dist", value<double>(&max_merge_dist)->default_value(0.03), "Max distance between two canopy centers in which the canopies should be merged")
        ("min_step_dist", value<double>(&min_step_dist)->default_value(0.1), "Min distance between canopy center and canopy centroid in which the centroid will be used as an origin for a new canpy");

    filter_in_options_desc.add_options()
        ("filter_min_non_zero_data_points", value<int>(&min_non_zero_data_samples)->default_value(3), "Use in the analysis only those points that have at least N non zero data points. Setting it to 0 will disable the filter");

    filter_out_options_desc.add_options()
        ("filter_zero_medians", value<int>(&min_num_non_zero_medians)->default_value(4), "Return only those canopies that have at least N non-zero medians. Setting it to 0 will disable the filter.")
        ("filter_single_point", value<double>(&max_single_data_point_proportion)->default_value(0.9), "Don't return canopies containing a single median which divided by sum of all its medians is greater than X. Setting it to 1 disables the filter.");

    early_stop_options_desc.add_options()
        ("stop_on_proportion_of_points_clustered", value<double>(&stop_proportion_of_points)->default_value(0.5), "Stop clustering when X*total_number_of_points were clustered")
        ("stop_on_num_single_point_clusters", value<int>(&stop_num_single_point_clusters)->default_value(1000), "Stop clustering when X consecutive clusters had only one point in them");

    misc_options_desc.add_options()
        ("save_canopy_size_statistics", value<string>(&canopy_size_stats_fp)->default_value(""), "Provide path to file to which statistics of canopy size vs number of canopies created will be saved")
        ("print_time_statistics,t", bool_switch(&print_time_statistics), "Print wall clock time profiles of various analysis parts. This is not aggressive and won't increase compuatation time.")
        ("show_progress_bar,b", bool_switch(&show_progress_bar), "Show progress bar, nice if output is printed to console, don't use if you are redirecting to a file. Verbosity must be set to at least PROGRESS for it to have an effect.") 
        ("help", "write help message");

    all_options_desc.add(general_options_desc).add(algorithm_param_options_desc).add(filter_in_options_desc).add(filter_out_options_desc).add(early_stop_options_desc).add(misc_options_desc);

    positional_options_description command_line_positional_desc;
    command_line_positional_desc.add("point_input_file",1);
    command_line_positional_desc.add("output_file_path",2);
    //command_line_positional_desc.add("output_centers_file_path",3);

    variables_map command_line_variable_map;
    store(command_line_parser(argc,argv).options(all_options_desc).positional(command_line_positional_desc).run(), command_line_variable_map);
    notify(command_line_variable_map);

    //
    //Verify command line input parameters
    //
    //verify_input_correctness(all_options_desc, command_line_variable_map);
    if (command_line_variable_map.count("help")) {
        cout << "Usage: cc.bin [options] POINTS_INPUT_FILE CLUSTERS_OUTPUT_FILE" << endl << endl;;
        cout << all_options_desc<< "\n";
        exit(1);
    }

    check_if_file_is_readable("point_input_file",point_input_file);
    check_if_file_is_writable("output_file_path",output_file_path);
    check_if_file_is_writable("output_centers_file_path",output_file_path);
    vector<string> valid_verbosities;
    valid_verbosities += "error", "progress", "warn", "info", "debug", "debug1", "debug2";
    check_if_one_of("verbosity_option",verbosity_option, valid_verbosities);
    check_if_within_bounds("num_threads",num_threads,1,999);//Not exactly future proof, but let's put foolproofness first
    check_if_within_bounds("max_canopy_dist",max_canopy_dist,0.0,1.0);
    check_if_within_bounds("max_close_dist",max_close_dist,0.0,1.0);
    check_if_within_bounds("max_merge_dist",max_merge_dist,0.0,1.0);
    check_if_within_bounds("min_step_dist",min_step_dist,0.0,1.0);

    check_if_within_bounds("min_non_zero_data_samples",min_non_zero_data_samples,0,10000);
    check_if_within_bounds("min_num_non_zero_medians",min_num_non_zero_medians,0,10000);
    check_if_within_bounds("max_single_data_point_proportion",max_single_data_point_proportion,0.0,1.0);
    check_if_within_bounds("stop_proportion_of_points",stop_proportion_of_points,0.0,1.0);
    check_if_within_bounds("stop_num_single_point_clusters",stop_num_single_point_clusters,0,10000000);
    if(canopy_size_stats_fp != "")
        check_if_file_is_writable("canopy_size_stats_fp",canopy_size_stats_fp);

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
    time_profile.start_timer("File loading");

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

    time_profile.stop_timer("File loading");
    time_profile.start_timer("Reading points");

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

    time_profile.stop_timer("Reading points");

    _log(logINFO) << "Points read, dropping file from memory";

    /* drop the file from memory*/
    munmap(point_file_mmap, statbuf.st_size);
    
    _log(logINFO) << "Running basic validation of points";

    time_profile.start_timer("Point validation");
    verify_proper_point_input_or_die(points);
    time_profile.stop_timer("Point validation");


    time_profile.start_timer("Input point filtering");
    if(min_non_zero_data_samples > 0){
        _log(logINFO) << "Filtering points";
        filter_out_input_points(points, min_non_zero_data_samples);
    }

    time_profile.stop_timer("Input point filtering");

    _log(logINFO) << "Finished input points processing";
    
    _log(logINFO) << "Number of points read: " << points.size();
    
    //
    //Run Canopy Clustering
    //
    std::vector<Canopy*> canopies;

    canopies = CanopyClusteringAlg::multi_core_run_clustering_on(points, num_threads, max_canopy_dist, max_close_dist, max_merge_dist, min_step_dist, stop_proportion_of_points, stop_num_single_point_clusters, canopy_size_stats_fp, show_progress_bar, time_profile);

    _log(logINFO) << "Finished clustering, number of canopies:" << canopies.size();

    //
    //Filter out canopies
    //

    if(min_num_non_zero_medians){
        time_profile.start_timer("Filtering canopies by medians");
        CanopyClusteringAlg::filter_clusters_by_zero_medians(min_num_non_zero_medians, canopies);
        _log(logINFO) << "Finished filtering by medians, number of canopies:" << canopies.size();
        time_profile.stop_timer("Filtering canopies by medians");
    }


    if(max_single_data_point_proportion < 0.99999){ //It's due to a double comparison
        time_profile.start_timer("Filtering canopies by single point bias");
        CanopyClusteringAlg::filter_clusters_by_single_point_skew(max_single_data_point_proportion, canopies);
        _log(logINFO) << "Finished filtering by single data point proportion, number of canopies:" << canopies.size();
        time_profile.stop_timer("Filtering canopies by single point bias");
    }


    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "####################Writing Results####################" ;
    ofstream output_file;

    sort(canopies.begin(), canopies.end(), compare_canopy_ptrs);

    int num_digits = ceil(log10(canopies.size()));
    cout << std::setfill('0');


    int i =0;
    output_file.open(output_file_path.c_str(), ios::out | ios::trunc);
    BOOST_FOREACH(Canopy* c, canopies){
        BOOST_FOREACH(Point* p, c->neighbours){
            output_file << output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << i++ << "\t";
            output_file << p->id << endl;
        }
    }
    output_file.close();

    i=0;
    output_file.open(output_centers_file_path.c_str(), ios::out | ios::trunc);
    BOOST_FOREACH(Canopy* c, canopies){
        output_file << output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << i++ << "\t";
        
        for(int j=0; j < c->center->num_data_samples; j++){
            output_file << c->center->sample_data[j] << "\t" ;
        }

        output_file << endl;
    }
    output_file.close();

    BOOST_FOREACH(Canopy* c, canopies)
        delete c;

    BOOST_FOREACH(Point* point, points)
        if(point)//Some points were centers of canopies and were deleted already
            delete point;


    time_profile.stop_timer("Total");
    //Write output statistics
    if(print_time_statistics){
        cout << time_profile << endl;
    }


    return 0;
}
