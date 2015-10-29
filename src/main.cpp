/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
 * Copyright (C) 2015 Enterome
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


#include <boost/foreach.hpp>

#include <Options.hh>

#include <Point.hpp>
#include <CanopyClustering.hpp>
#include <Log.hpp>
#include <signal_handlers.hpp>

#include <Stats.hpp>
#include <omp.h>

using namespace std;



int main(int argc, char* argv[])
{
    //Set initial logging level
	Logger::set_log_level("info");

    //Prepare Time Profile
    TimeProfile time_profile;
    time_profile.start_timer("Total");


    //Parse command line options
	const Options& options = Options::parse(argc, argv);
	_log(logINFO) << options;


    //Set user chosen logging level
	Logger::set_log_level(options.verbosity_option);

    //Set signal handler
    if(options.die_on_kill) 
        signal(SIGINT, signal_callback_die_handler);
    else    
        signal(SIGINT, signal_callback_gentle_handler);

    //Set number of threads
    omp_set_num_threads(options.num_threads);

    //
    //Parse point description file
    //

    vector<Point*> points;
    vector<Point*> filtered_points;

	std::ifstream point_file(options.point_input_file.c_str());    

    time_profile.start_timer("Reading points");

    _log(logINFO) << "Reading point input file.";

	// Should be sufficient for the biggest count matrices we have
	char line[100000];

	while (point_file.getline(line, 100000)) {

		char* line_ptr = line;

		// Skip blank spaces at the beginning of the line
		while (*line_ptr != '\0' && std::isspace(*line_ptr))
		{
			line_ptr++;
		}

		// Line is empty
		if (*line_ptr == '\0')
			continue;

		points.push_back(new Point(line_ptr));

		die_if_true(terminate_called);
	}

    _log(logINFO) << "Finished reading point input file.";
    _log(logINFO) << "";

    time_profile.stop_timer("Reading points");

	if (points.empty()) {
        _log(logERR) << "Input has no data.";
		exit(1);
	}

    _log(logINFO) << "Running basic validation of points";
    _log(logINFO) << "";

    time_profile.start_timer("Point validation");
    verify_proper_point_input_or_die(points);
    time_profile.stop_timer("Point validation");

    vector<Point*> points_filtered_out;
    vector<Point*> points_filtered_out_due_to_three_point_proportion_filter;
    vector<Point*> points_filtered_out_due_to_num_non_zero_samples_filter;

    time_profile.start_timer("Input point filtering");

#pragma omp parallel for shared(points_filtered_out_due_to_three_point_proportion_filter, points_filtered_out_due_to_num_non_zero_samples_filter, filtered_points)
    for(size_t i = 0; i < points.size(); i++){
        //Both filters are set
        if((options.min_non_zero_data_samples > 0) && (options.max_top_three_data_point_proportion < 0.9999)){
            bool point_is_valid = true;
                
            if( ! points[i]->check_if_num_non_zero_samples_is_greater_than_x(options.min_non_zero_data_samples))
            {
#pragma omp critical
                points_filtered_out_due_to_num_non_zero_samples_filter.push_back(points[i]);
                point_is_valid = false;
            }

            if( ! points[i]->check_if_top_three_point_proportion_is_smaller_than(options.max_top_three_data_point_proportion))
            {
#pragma omp critical
                points_filtered_out_due_to_three_point_proportion_filter.push_back(points[i]);
                point_is_valid = false;
            }

            if(point_is_valid){
#pragma omp critical
                filtered_points.push_back(points[i]);
            } else {
#pragma omp critical
				points_filtered_out.push_back(points[i]);
			}
        } else if (options.min_non_zero_data_samples > 0){ 
            if(points[i]->check_if_num_non_zero_samples_is_greater_than_x(options.min_non_zero_data_samples)){
#pragma omp critical
                filtered_points.push_back(points[i]);
            }
            else 
#pragma omp critical
            {
                points_filtered_out_due_to_num_non_zero_samples_filter.push_back(points[i]);
				points_filtered_out.push_back(points[i]);
            }
        } else if (options.max_top_three_data_point_proportion < 0.9999){ 
            if(points[i]->check_if_top_three_point_proportion_is_smaller_than(options.max_top_three_data_point_proportion)){ 
#pragma omp critical
                filtered_points.push_back(points[i]);
            }
            else 
#pragma omp critical
            {
                points_filtered_out_due_to_three_point_proportion_filter.push_back(points[i]);
				points_filtered_out.push_back(points[i]);
            }
        }
    }
    if(options.points_filtered_out_at_least_non_zero_file_path != ""){
        ofstream filtered_point_file;
        filtered_point_file.open(options.points_filtered_out_at_least_non_zero_file_path.c_str(), ios::out | ios::trunc);
        for(size_t i = 0; i < points_filtered_out_due_to_num_non_zero_samples_filter.size(); i++){
            filtered_point_file << points[i]->id << "\n";
        }
        filtered_point_file.close();
    }

    if(options.points_filtered_out_top_three_prop_file_path != ""){
        ofstream filtered_point_file;
        filtered_point_file.open(options.points_filtered_out_top_three_prop_file_path.c_str(), ios::out | ios::trunc);
        for(size_t i = 0; i < points_filtered_out_due_to_three_point_proportion_filter.size(); i++){
            filtered_point_file << points[i]->id << "\n";
        }
        filtered_point_file.close();
    }
	
	BOOST_FOREACH(Point* p, points_filtered_out)
		delete p;

    time_profile.stop_timer("Input point filtering");
    _log(logINFO) << "Number of points filtered out due to three point sample values proportion filter: " << points_filtered_out_due_to_three_point_proportion_filter.size();
    _log(logINFO) << "Number of points filtered out due to non zero samples number filter: " << points_filtered_out_due_to_num_non_zero_samples_filter.size();
    _log(logINFO) << "Number of points filtered out: " << points_filtered_out.size(); 

    _log(logINFO) << "Finished input points processing";
    
    _log(logINFO) << "Number of points after filtering: " << filtered_points.size();
    _log(logINFO) << "";
    
    die_if_true(terminate_called);
    die_if_true(filtered_points.size() < 1);


    time_profile.start_timer("Precomputing pearson correlations.");
    _log(logINFO) << "Precomputing pearson correlations.";

	for (size_t curr_point = 0; curr_point < filtered_points.size(); curr_point++)
	{
		Point* point = filtered_points[curr_point];
		point->sample_data_pearson_precomputed = new double[point->num_data_samples];
	}

	#pragma omp parallel for
	for (size_t curr_point = 0; curr_point < filtered_points.size(); curr_point++)
	{
		Point* point = filtered_points[curr_point];
		precompute_pearson_data(point->num_data_samples, point->sample_data, point->sample_data_pearson_precomputed);
	}
    _log(logINFO) << "Done.";
    _log(logINFO) << "";
    time_profile.stop_timer("Precomputing pearson correlations.");

    //
    //Run Canopy Clustering
    //
    std::vector<Canopy*> canopies;

    canopies = CanopyClusteringAlg::multi_core_run_clustering_on(filtered_points, options.max_canopy_dist, options.max_close_dist, options.max_merge_dist, options.min_step_dist, options.max_num_canopy_walks, options.stop_proportion_of_points, options.canopy_size_stats_file, options.not_processed_points_file, options.show_progress_bar, time_profile);

    _log(logINFO) << "Finished clustering, number of canopies:" << canopies.size();

    //
    //Filter out canopies
    //

    if(options.min_num_non_zero_medians){
        time_profile.start_timer("Filtering canopies by medians");
        CanopyClusteringAlg::filter_clusters_by_zero_medians(options.min_num_non_zero_medians, canopies);
        _log(logINFO) << "Finished filtering by medians, number of canopies:" << canopies.size();
        time_profile.stop_timer("Filtering canopies by medians");
    }


    if(options.max_single_data_point_proportion < 0.99999){ //It's due to a double comparison
        time_profile.start_timer("Filtering canopies by single point bias");
        CanopyClusteringAlg::filter_clusters_by_single_point_skew(options.max_single_data_point_proportion, canopies);
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

    int num_digits = ceil(log10(canopies.size()+1));
    cout << std::setfill('0');


    int i =1;
    out_file.open(options.output_file.c_str(), ios::out | ios::trunc);
    BOOST_FOREACH(Canopy* c, canopies){
		std::ostringstream cluster_name;
		cluster_name << options.output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << i;

        BOOST_FOREACH(Point* p, c->neighbours){
            out_file << cluster_name.str() << "\t" << p->id << "\n";
        }
        i++;
    }
    out_file.close();

    i=1;
    out_file.open(options.output_centers_file.c_str(), ios::out | ios::trunc);
    BOOST_FOREACH(Canopy* c, canopies){
        out_file << options.output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << i << "\t";
        
        for(size_t j=0; j < c->center->num_data_samples; j++){
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

    BOOST_FOREACH(Point* point, filtered_points)
        if(point)//Some points were centers of canopies and were deleted already
            delete point;


    time_profile.stop_timer("Total");
    //Write output statistics
    if(options.print_time_statistics){
        cout << time_profile << endl;
    }


    return 0;
}
