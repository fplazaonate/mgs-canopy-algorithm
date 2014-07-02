/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
 *
 * This file is part of Metagenomics Canopy Clustering Implementation.
 *
 * Metagenomics Canopy Clustering Implementation is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Metagenomics Canopy Clustering Implementation is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */
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
#include <signal_handlers.hpp>

#include <omp.h>

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
    string points_filtered_out_top_three_prop_file_path;
    string points_filtered_out_at_least_non_zero_file_path;
    string output_file;
    string output_centers_file;
    string output_cluster_prefix;
    int num_threads;
    double max_canopy_dist;
    double max_close_dist;
    double max_merge_dist;
    double min_step_dist;
    string verbosity_option;
    int min_non_zero_data_samples;
    double max_top_three_data_point_proportion;
    int min_num_non_zero_medians;
    double max_single_data_point_proportion;
    double stop_proportion_of_points;
    string canopy_size_stats_file;
    string not_processed_points_file;
    bool show_progress_bar;
    bool print_time_statistics;
    bool die_on_kill;
    int max_num_canopy_walks;


    //Define and read command line options
    options_description all_options_desc("Allowed options");
    options_description general_options_desc("General");
    options_description algorithm_param_options_desc("Algorithm Parameters");
    options_description filter_in_options_desc("Input filter parameters");
    options_description filter_out_options_desc("Output filter parameters");
    options_description early_stop_options_desc("Early stopping");
    options_description misc_options_desc("Miscellaneous");


    general_options_desc.add_options()
        ("input_file_path,i", value<string>(&point_input_file), "Path to the input file")
        ("output_clusters_file_path,o", value<string>(&output_file), "Path to file to which clusters will be written")
        ("output_cluster_profiles_file,c", value<string>(&output_centers_file), "Path to file to which cluster profiles will be written")
        ("cluster_name_prefix,p", value<string>(&output_cluster_prefix)->default_value("CAG"), "Prefix prepended to output cluster names")
        ("num_threads,n", value<int>(&num_threads)->default_value(4), "IMPORTANT! Number of cpu threads to use.")
        ("verbosity,v", value<string>(&verbosity_option)->default_value("info"), "Control how much information should be printed to the screen. Available levels according to their verbosity: error, progress, warn, info, debug, debug1.");

    algorithm_param_options_desc.add_options()
        ("max_canopy_dist", value<double>(&max_canopy_dist)->default_value(0.1), "Max pearson correlation difference between a canopy center and a point included to the canopy")
        ("max_close_dist", value<double>(&max_close_dist)->default_value(0.4), "Max pearson correlation difference between a canopy center and a point in which the point will be considered close to the canopy. As a heuristc, only points within this distance will be considered as potential neighbours during the canopy walk.")
        ("max_merge_dist", value<double>(&max_merge_dist)->default_value(0.05), "Max pearson correlation difference between two canopy centers in which the canopies should be merged. Please note, that the final canopy profiles are calculated after the merge step and consequently some final canopies might have profiles that are closer then max_merge_dist specifies.")
        ("min_step_dist", value<double>(&min_step_dist)->default_value(0.01), "Min pearson correlation difference between canopy center and canopy centroid in which the centroid will be used as an origin for a new canpy (canopy walk). This is a stop criterion for canopy walk.")
        ("max_num_canopy_walks", value<int>(&max_num_canopy_walks)->default_value(3), "Max number of times the canopy will walk. This is a stop criterion for canopy walk.");

    filter_in_options_desc.add_options()
        ("filter_min_obs", value<int>(&min_non_zero_data_samples)->default_value(3), "Discard those points which have fewer than N non-zero data points (observations). Setting it to 0 will disable the filter.")
        ("filter_max_dominant_obs", value<double>(&max_top_three_data_point_proportion)->default_value(0.9), "Discard those points for which top 3 data points constitute more than X fraction of the total signal. Setting it to 1 will disable the filter")
        ("filtered_out_points_min_obs_file", value<string>(&points_filtered_out_at_least_non_zero_file_path)->default_value(""), "File to which write out those files that didn't match the filter_min_obs filter")
        ("filtered_out_points_max_dominant_obs_file", value<string>(&points_filtered_out_top_three_prop_file_path)->default_value(""), "File to which write out those files that didn't match the filter_max_dominant_obs filter.");

    filter_out_options_desc.add_options()
        ("filter_zero_medians", value<int>(&min_num_non_zero_medians)->default_value(3), "Return only those canopies that have at least N non-zero cluster profile observations. Setting it to 0 will disable the filter.")
        ("filter_single_point", value<double>(&max_single_data_point_proportion)->default_value(0.9), "Don't return canopies containing a single profile observation which constitutes to more than X fraction of the total profile. Setting it to 1 disables the filter.");

    early_stop_options_desc.add_options()
        ("stop_fraction", value<double>(&stop_proportion_of_points)->default_value(1.0), "Stop clustering after X fraction of all points have been clustered. Setting it to 1 will disable this stop criterion.");

    misc_options_desc.add_options()
        ("die_on_kill", bool_switch(&die_on_kill), "If set, after receiving a KILL signal, the program will die and no results will be produced. By default clustering will stop but clusters will be merged and partial results will be printed as usual.")
        ("not_processed_points_file", value<string>(&not_processed_points_file)->default_value(""), "Path to file to which unprocessed origins will be dumped at KILL signal")
        ("print_time_statistics,t", bool_switch(&print_time_statistics), "Print wall clock time profiles of various analysis parts. This is not aggressive and won't increase compuatation time.")
        ("show_progress_bar,b", bool_switch(&show_progress_bar), "Show progress bar, nice if output is printed to console, don't use if you are redirecting to a file. Verbosity must be set to at least PROGRESS for it to have an effect.") 
        ("canopy_size_stats_file", value<string>(&canopy_size_stats_file)->default_value(""), "If set, to this file current progress after each processed origin will be dumped in format <index> <num_left_origins> <this_canopy_size> <total_num_thread_collisions>")
        ("help", "write help message");

    all_options_desc.add(general_options_desc).add(algorithm_param_options_desc).add(filter_in_options_desc).add(filter_out_options_desc).add(early_stop_options_desc).add(misc_options_desc);

    positional_options_description command_line_positional_desc;
    command_line_positional_desc.add("point_input_file",1);
    command_line_positional_desc.add("output_file",1);
    command_line_positional_desc.add("output_centers_file",1);

    variables_map command_line_variable_map;
    store(command_line_parser(argc,argv).options(all_options_desc).positional(command_line_positional_desc).run(), command_line_variable_map);
    notify(command_line_variable_map);

    //
    //Verify command line input parameters
    //
    //verify_input_correctness(all_options_desc, command_line_variable_map);
    if (command_line_variable_map.count("help") || argc < 3) {
        cout << "Usage: cc.bin [options] POINTS_INPUT_FILE CLUSTERS_OUTPUT_FILE" << endl << endl;;
        cout << all_options_desc<< "\n";
        exit(1);
    }

    check_if_file_is_readable("point_input_file",point_input_file);
    check_if_file_is_writable("output_file",output_file);
    check_if_file_is_writable("output_centers_file",output_centers_file);
    check_if_file_is_writable("points_filtered_out_top_three_prop_file_path",points_filtered_out_top_three_prop_file_path);
    check_if_file_is_writable("points_filtered_out_at_least_non_zero_file_path",points_filtered_out_at_least_non_zero_file_path);
    vector<string> valid_verbosities;
    valid_verbosities += "error", "progress", "warn", "info", "debug", "debug1", "debug2", "debug3";
    check_if_one_of("verbosity_option",verbosity_option, valid_verbosities);
    check_if_within_bounds("num_threads",num_threads,1,999);//Not exactly future proof, but let's put foolproofness first
    check_if_within_bounds("max_canopy_dist",max_canopy_dist,0.0,1.0);
    check_if_within_bounds("max_close_dist",max_close_dist,0.0,1.0);
    check_if_within_bounds("max_merge_dist",max_merge_dist,0.0,1.0);
    check_if_within_bounds("min_step_dist",min_step_dist,0.0,1.0);
    check_if_within_bounds("max_num_canopy_walks",max_num_canopy_walks,0,100);

    check_if_within_bounds("min_non_zero_data_samples",min_non_zero_data_samples,0,10000);
    check_if_within_bounds("max_top_three_data_point_proportion",max_top_three_data_point_proportion,0.0,1.0);
    check_if_within_bounds("min_num_non_zero_medians",min_num_non_zero_medians,0,10000);
    check_if_within_bounds("max_single_data_point_proportion",max_single_data_point_proportion,0.0,1.0);
    check_if_within_bounds("stop_proportion_of_points",stop_proportion_of_points,0.0,1.0);
    if(canopy_size_stats_file != "")
        check_if_file_is_writable("canopy_size_stats_file",canopy_size_stats_file);
    if(not_processed_points_file!= "")
        check_if_file_is_writable("not_processed_points_file",not_processed_points_file);
       
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
    }else if(verbosity_option == "debug2"){
        log_level = logDEBUG2;
    }else if(verbosity_option == "debug3"){
        log_level = logDEBUG3;
    }

    _log(logINFO) << "";
    _log(logINFO) << "Files:";
    _log(logINFO) << "point_input_file:\t " << point_input_file;
    _log(logINFO) << "output_centers_file:\t " << output_centers_file;
    _log(logINFO) << "canopy_size_stats_file:\t " << canopy_size_stats_file;
    _log(logINFO) << "not_processed_points_file:\t " << not_processed_points_file;
    _log(logINFO) << "";
    

    //Set signal handler
    if(die_on_kill) 
        signal(SIGINT, signal_callback_die_handler);
    else    
        signal(SIGINT, signal_callback_gentle_handler);

    //Set number of threads
    _log(logINFO) << "";
    _log(logINFO) << "General:";
    _log(logINFO) << "num_threads:\t " << num_threads;
    _log(logINFO) << "";

    omp_set_num_threads(num_threads);

    //
    //Parse point description file
    //
    time_profile.start_timer("File loading");

    vector<Point*> points;
    vector<Point*> filtered_points;

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
        if(line_end_ptr != line_start_ptr){//Check if the line is not empty
            *line_end_ptr = '\0';
            //cout << line_start_ptr << endl;
            points.push_back(new Point(line_start_ptr));
        }
        line_start_ptr = ++line_end_ptr;
        die_if_true(terminate_called);
    }

    time_profile.stop_timer("Reading points");

    _log(logINFO) << "Points read, dropping file from memory";
    _log(logINFO) << "";

    /* drop the file from memory*/
    munmap(point_file_mmap, statbuf.st_size);
    
    _log(logINFO) << "Running basic validation of points";
    _log(logINFO) << "max_top_three_data_point_proportion:\t " << max_top_three_data_point_proportion;
    _log(logINFO) << "min_non_zero_data_samples:\t " << min_non_zero_data_samples;
    _log(logINFO) << "";

    time_profile.start_timer("Point validation");
    verify_proper_point_input_or_die(points);
    time_profile.stop_timer("Point validation");

    vector<Point*> points_filtered_out_due_to_three_point_proportion_filter;
    vector<Point*> points_filtered_out_due_to_num_non_zero_samples_filter;

    time_profile.start_timer("Input point filtering");

#pragma omp parallel for shared(points_filtered_out_due_to_three_point_proportion_filter, points_filtered_out_due_to_num_non_zero_samples_filter, filtered_points)
    for(int i = 0; i < points.size(); i++){
        //Both filters are set
        if((min_non_zero_data_samples > 0) && (max_top_three_data_point_proportion < 0.9999)){
            bool point_is_valid = true;
                
            if( ! points[i]->check_if_num_non_zero_samples_is_greater_than_x(min_non_zero_data_samples))
            {
#pragma omp critical
                points_filtered_out_due_to_num_non_zero_samples_filter.push_back(points[i]);
                point_is_valid = false;
            }

            if( ! points[i]->check_if_top_three_point_proportion_is_smaller_than(max_top_three_data_point_proportion))
            {
#pragma omp critical
                points_filtered_out_due_to_three_point_proportion_filter.push_back(points[i]);
                point_is_valid = false;
            }

            if(point_is_valid){
#pragma omp critical
                filtered_points.push_back(points[i]);
            }
        } else if (min_non_zero_data_samples > 0){ 
            if(points[i]->check_if_num_non_zero_samples_is_greater_than_x(min_non_zero_data_samples)){
#pragma omp critical
                filtered_points.push_back(points[i]);
            }
            else 
            {
#pragma omp critical
                points_filtered_out_due_to_num_non_zero_samples_filter.push_back(points[i]);
            }
        } else if (max_top_three_data_point_proportion < 0.9999){ 
            if(points[i]->check_if_top_three_point_proportion_is_smaller_than(max_top_three_data_point_proportion)){ 
#pragma omp critical
                filtered_points.push_back(points[i]);
            }
            else 
            {
#pragma omp critical
                points_filtered_out_due_to_three_point_proportion_filter.push_back(points[i]);
            }
        }
    }
    if(points_filtered_out_at_least_non_zero_file_path != ""){
        ofstream filtered_point_file;
        filtered_point_file.open(points_filtered_out_at_least_non_zero_file_path.c_str(), ios::out | ios::trunc);
        for(int i = 0; i < points_filtered_out_due_to_num_non_zero_samples_filter.size(); i++){
            filtered_point_file << points[i]->id << "\n";
        }
        filtered_point_file.close();
    }

    if(points_filtered_out_top_three_prop_file_path != ""){
        ofstream filtered_point_file;
        filtered_point_file.open(points_filtered_out_top_three_prop_file_path.c_str(), ios::out | ios::trunc);
        for(int i = 0; i < points_filtered_out_due_to_three_point_proportion_filter.size(); i++){
            filtered_point_file << points[i]->id << "\n";
        }
        filtered_point_file.close();
    }

    time_profile.stop_timer("Input point filtering");
    _log(logINFO) << "Number of points filtered out due to three point sample values proportion filter: " << points_filtered_out_due_to_three_point_proportion_filter.size();
    _log(logINFO) << "Number of points filtered out due to non zero samples number filter: " << points_filtered_out_due_to_num_non_zero_samples_filter.size();
    _log(logINFO) << "Number of points filtered out: " << points.size() - filtered_points.size(); 

    _log(logINFO) << "Finished input points processing";
    
    _log(logINFO) << "Number of points after filtering: " << filtered_points.size();
    
    die_if_true(terminate_called);
    die_if_true(filtered_points.size() < 1);
    //
    //Run Canopy Clustering
    //
    std::vector<Canopy*> canopies;

    canopies = CanopyClusteringAlg::multi_core_run_clustering_on(filtered_points, num_threads, max_canopy_dist, max_close_dist, max_merge_dist, min_step_dist, max_num_canopy_walks, stop_proportion_of_points, canopy_size_stats_file, not_processed_points_file, show_progress_bar, time_profile);

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

    {
        time_profile.start_timer("Filtering canopies by size");
        CanopyClusteringAlg::filter_clusters_by_size(canopies);
        _log(logINFO) << "Finished filtering by size(number of neighbours must be bigger than 1), number of canopies:" << canopies.size();
        time_profile.stop_timer("Filtering canopies by size");
    }

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "####################Writing Results####################" ;
    ofstream out_file;

    sort(canopies.begin(), canopies.end(), compare_canopy_ptrs_by_canopy_size);

    int num_digits = ceil(log10(canopies.size()));
    cout << std::setfill('0');


    int i =1;
    out_file.open(output_file.c_str(), ios::out | ios::trunc);
    BOOST_FOREACH(Canopy* c, canopies){
        BOOST_FOREACH(Point* p, c->neighbours){
            out_file << output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << i << "\t";
            out_file << p->id << "\n";
        }
        i++;
    }
    out_file.close();

    i=1;
    out_file.open(output_centers_file.c_str(), ios::out | ios::trunc);
    BOOST_FOREACH(Canopy* c, canopies){
        out_file << output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << i << "\t";
        
        for(int j=0; j < c->center->num_data_samples; j++){
            out_file << c->center->sample_data[j] << "\t" ;
        }

        i++;
        out_file << "\n";
    }
    out_file.close();

    //
    //Clean up
    //
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
