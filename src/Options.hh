#ifndef OPTIONS_HH
#define OPTIONS_HH

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

#include <string>
#include <iostream>
#include <vector>

struct Options
{
    std::string point_input_file;
    std::string points_filtered_out_top_three_prop_file_path;
    std::string points_filtered_out_at_least_non_zero_file_path;
    std::string output_file;
    std::string output_centers_file;
    std::string output_cluster_prefix;
    int num_threads;
    double max_canopy_dist;
    double max_close_dist;
    double max_merge_dist;
    double min_step_dist;
    std::string verbosity_option;
    int min_non_zero_data_samples;
    double max_top_three_data_point_proportion;
    int min_num_non_zero_medians;
    double max_single_data_point_proportion;
    double stop_proportion_of_points;
    std::string canopy_size_stats_file;
    std::string not_processed_points_file;
    bool show_progress_bar;
    bool print_time_statistics;
    bool die_on_kill;
    int max_num_canopy_walks;

	static Options parse(int argc, char* argv[]);
	static bool check_if_file_is_readable(const std::string& option_name, const std::string& path);
	static bool check_if_file_is_writable(const std::string& option_name, const std::string& path);
	static bool check_if_one_of(const std::string& option_name, const std::string& value, const std::vector<std::string>& valid_options);
	template <typename T>
	static bool check_if_within_bounds(const std::string& val_name, T value, T lower, T higher);
};

std::ostream& operator<<(std::ostream& os, const Options& parameters);

#endif // OPTIONS_HH
