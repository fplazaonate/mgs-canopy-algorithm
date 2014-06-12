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
#ifndef CANOPY_CLUSTERING
#define CANOPY_CLUSTERING
#include <vector>
#include <boost/unordered_set.hpp>

#include <Point.hpp>
#include <Canopy.hpp>
#include <TimeProfile.hpp>

using namespace std;

/*
 * Static class representing the Canopy Clustering Algorithm
 */

class CanopyClusteringAlg{
    public:
        /**
         * Run the canopy clustering algorithm
         *
         * Parameters:
         * points - list of references to points to be clustered
         * num_threads - number of threads to be used in the calculation
         * max_canopy_dist, max_close_dist, max_merge_dist, max_num_canopy_walks, stop_proportion_of_points, show_progress_bar - see program parameters description 
         * canopy_size_stats_fp - absolute file path to the canopy size statistics file
         * not_procesed_point_fp - absolute file path to the file containing not processed points if early stopped
         * time_profile - TimeProfile object instance for gathering statistics on time it took for each of the analysis steps
         */
        static std::vector<Canopy*> multi_core_run_clustering_on(vector<Point*>& points, int num_threads, double max_canopy_dist, double max_close_dist, double max_merge_dist, double min_step_dist, int max_num_canopy_walks, double stop_proportion_of_points, string canopy_size_stats_fp, string not_processed_points_fp,  bool show_progress_bar, TimeProfile& time_profile);

        /**
         * Create canopy given an origin point
         *
         * Parameters:
         * origin - point being the canopy origin
         * points - list of all points
         * close_points - list of close_points
         * min_neighbour_correlation - minimum distance in correlation space between origin and a tested point for the point to be considered inside the canopy
         * min_close_correlation - minimum distance in correlation space between origin and a tested point for the point to be considered close to the canopy
         * sets_close_points - flag describing if the current execution of this function should set the close_points
         */
        static Canopy* create_canopy(Point* origin, vector<Point*>& points, vector<Point*>& close_points, double min_neighbour_correlation, double min_close_correlation, bool sets_close_points);
        
        /**
         * Execute the create_canopy function iteratively until a stable canopy is reached. 
         *
         * Parameters:
         * origin - point which is the canopy origin
         * points - list of all points
         * close_points - list of points which are within "close" distance to canopy
         * max_canopy_dist, max_close_dist, min_step_dist, max_num_canopy_walks  - see program parameters
         * num_canopy_jumps - number of times the create_canopy function was executed
         */
        static Canopy* canopy_walk(Point* origin, vector<Point*>& points, vector<Point*>& close_points, double max_canopy_dist, double max_close_dist, double min_step_dist, double max_num_canopy_walks, int& num_canopy_jumps);

        static void filter_clusters_by_zero_medians(int min_num_non_zero_medians, vector<Canopy*>& canopies_to_filter);
        static void filter_clusters_by_single_point_skew(double max_single_data_point_proportion, vector<Canopy*>& canopies_to_filter);
        static void filter_clusters_by_size(std::vector<Canopy*>& canopies_to_filter);
};

#endif 
