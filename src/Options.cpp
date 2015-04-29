#include "Options.hh"
#include <fstream>
#include <limits>
#include <omp.h>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <boost/assign/list_of.hpp> 
namespace po = boost::program_options;
#include <boost/algorithm/string/join.hpp>
#include <Log.hpp>


const std::vector<std::string> Options::valid_verbosities = 
boost::assign::list_of("error")("progress")("warn")("info")("debug")("debug1")("debug2")("debug3");

Options Options::parse(int argc, char* argv[])
{
	Options options;

	//Define and read command line options
	po::options_description all_options_desc("Allowed options");
	po::options_description general_options_desc("General");
	po::options_description algorithm_param_options_desc("Algorithm Options");
	po::options_description filter_in_options_desc("Input filter options");
	po::options_description filter_out_options_desc("Output filter options");
	po::options_description early_stop_options_desc("Early stopping");
	po::options_description misc_options_desc("Miscellaneous");


	general_options_desc.add_options()
		("input_file_path,i", po::value<std::string>(&options.point_input_file), "Path to the input file")
		("output_clusters_file_path,o", po::value<std::string>(&options.output_file), "Path to file to which clusters will be written")
		("output_cluster_profiles_file,c", po::value<std::string>(&options.output_centers_file), "Path to file to which cluster profiles will be written")
		("cluster_name_prefix,p", po::value<std::string>(&options.output_cluster_prefix)->default_value("CAG"), "Prefix prepended to output cluster names")
		("num_threads,n", po::value<int>(&options.num_threads)->default_value(4), "IMPORTANT! Number of cpu threads to use.")
		("verbosity,v", po::value<std::string>(&options.verbosity_option)->default_value("info"), "Control how much information should be printed to the screen. Available levels according to their verbosity: error, progress, warn, info, debug, debug1.");

	algorithm_param_options_desc.add_options()
		("max_canopy_dist", po::value<double>(&options.max_canopy_dist)->default_value(0.1,"0.1"), "Max pearson correlation difference between a canopy center and a point included to the canopy")
		("max_close_dist", po::value<double>(&options.max_close_dist)->default_value(0.4,"0.4"), "Max pearson correlation difference between a canopy center and a point in which the point will be considered close to the canopy. As a heuristc, only points within this distance will be considered as potential neighbours during the canopy walk.")
		("max_merge_dist", po::value<double>(&options.max_merge_dist)->default_value(0.05,"0.05"), "Max pearson correlation difference between two canopy centers in which the canopies should be merged. Please note, that the final canopy profiles are calculated after the merge step and consequently some final canopies might have profiles that are closer then max_merge_dist specifies.")
		("min_step_dist", po::value<double>(&options.min_step_dist)->default_value(0.01,"0.01"), "Min pearson correlation difference between canopy center and canopy centroid in which the centroid will be used as an origin for a new canpy (canopy walk). This is a stop criterion for canopy walk.")
		("max_num_canopy_walks", po::value<int>(&options.max_num_canopy_walks)->default_value(3), "Max number of times the canopy will walk. This is a stop criterion for canopy walk.");

	filter_in_options_desc.add_options()
		("filter_min_obs", po::value<int>(&options.min_non_zero_data_samples)->default_value(3), "Discard those points which have fewer than N non-zero data points (observations). Setting it to 0 will disable the filter.")
		("filter_max_dominant_obs", po::value<double>(&options.max_top_three_data_point_proportion)->default_value(0.9,"0.9"), "Discard those points for which top 3 data points constitute more than X fraction of the total signal. Setting it to 1 will disable the filter")
		("filtered_out_points_min_obs_file", po::value<std::string>(&options.points_filtered_out_at_least_non_zero_file_path)->default_value(""), "File to which write out those files that didn't match the filter_min_obs filter")
		("filtered_out_points_max_dominant_obs_file", po::value<std::string>(&options.points_filtered_out_top_three_prop_file_path)->default_value(""), "File to which write out those files that didn't match the filter_max_dominant_obs filter.");

	filter_out_options_desc.add_options()
		("filter_zero_medians", po::value<int>(&options.min_num_non_zero_medians)->default_value(3), "Return only those canopies that have at least N non-zero cluster profile observations. Setting it to 0 will disable the filter.")
		("filter_single_point", po::value<double>(&options.max_single_data_point_proportion)->default_value(0.9,"0.9"), "Don't return canopies containing a single profile observation which constitutes to more than X fraction of the total profile. Setting it to 1 disables the filter.");

	early_stop_options_desc.add_options()
		("stop_fraction", po::value<double>(&options.stop_proportion_of_points)->default_value(1.0,"1.0"), "Stop clustering after X fraction of all points have been clustered. Setting it to 1 will disable this stop criterion.");

	misc_options_desc.add_options()
		("die_on_kill", po::bool_switch(&options.die_on_kill), "If set, after receiving a KILL signal, the program will die and no results will be produced. By default clustering will stop but clusters will be merged and partial results will be printed as usual.")
		("not_processed_points_file", po::value<std::string>(&options.not_processed_points_file)->default_value(""), "Path to file to which unprocessed origins will be dumped at KILL signal")
		("print_time_statistics,t", po::bool_switch(&options.print_time_statistics), "Print wall clock time profiles of various analysis parts. This is not aggressive and won't increase compuatation time.")
		("show_progress_bar,b", po::bool_switch(&options.show_progress_bar), "Show progress bar, nice if output is printed to console, don't use if you are redirecting to a file. Verbosity must be set to at least PROGRESS for it to have an effect.") 
		("canopy_size_stats_file", po::value<std::string>(&options.canopy_size_stats_file)->default_value(""), "If set, to this file current progress after each processed origin will be dumped in format <index> <num_left_origins> <this_canopy_size> <total_num_thread_collisions>")
		("help", "write help message");

	all_options_desc.add(general_options_desc).add(algorithm_param_options_desc).add(filter_in_options_desc).add(filter_out_options_desc).add(early_stop_options_desc).add(misc_options_desc);

	po::variables_map command_line_variable_map;
	po::store(po::command_line_parser(argc,argv).options(all_options_desc).run(), command_line_variable_map);
	po::notify(command_line_variable_map);

	//
	//Verify command line input options
	//
	//verify_input_correctness(all_options_desc, command_line_variable_map);
	if (command_line_variable_map.count("help") || argc < 3) {
		std::cout << all_options_desc<< "\n";
		exit(1);
	}

	check_if_file_is_readable("point_input_file", options.point_input_file);
	check_if_file_is_writable("output_file", options.output_file);
	check_if_file_is_writable("output_centers_file", options.output_centers_file);
	if (options.points_filtered_out_top_three_prop_file_path != "")
		check_if_file_is_writable("points_filtered_out_top_three_prop_file_path", options.points_filtered_out_top_three_prop_file_path);
	if (options.points_filtered_out_at_least_non_zero_file_path != "")
		check_if_file_is_writable("points_filtered_out_at_least_non_zero_file_path", options.points_filtered_out_at_least_non_zero_file_path);
	check_if_one_of("verbosity_option", options.verbosity_option, valid_verbosities);
	check_if_within_bounds("num_threads", options.num_threads,1,omp_get_max_threads());
	check_if_within_bounds("max_canopy_dist", options.max_canopy_dist,0.0,1.0);
	check_if_within_bounds("max_close_dist", options.max_close_dist,0.0,1.0);
	check_if_within_bounds("max_merge_dist", options.max_merge_dist,0.0,1.0);
	check_if_within_bounds("min_step_dist", options.min_step_dist,0.0,1.0);
	check_if_within_bounds("max_num_canopy_walks", options.max_num_canopy_walks,0,100);

	check_if_within_bounds("min_non_zero_data_samples", options.min_non_zero_data_samples,0,10000);
	check_if_within_bounds("max_top_three_data_point_proportion", options.max_top_three_data_point_proportion,0.0,1.0);
	check_if_within_bounds("min_num_non_zero_medians", options.min_num_non_zero_medians,0,10000);
	check_if_within_bounds("max_single_data_point_proportion", options.max_single_data_point_proportion,0.0,1.0);
	check_if_within_bounds("stop_proportion_of_points", options.stop_proportion_of_points,0.0,1.0);
	if(options.canopy_size_stats_file != "")
		check_if_file_is_writable("canopy_size_stats_file", options.canopy_size_stats_file);
	if(options.not_processed_points_file!= "")
		check_if_file_is_writable("not_processed_points_file", options.not_processed_points_file);

	return options;
}


template <typename T>
bool Options::check_if_within_bounds(const std::string& option_name, T value, T lower, T higher){
	if( value >= lower - std::numeric_limits<T>::epsilon() && value <= higher + std::numeric_limits<T>::epsilon()){
		return true;
	}else{ 
		_log(logERR) << "Option: \"" << option_name << "\" must be a value within range: <" << lower << ";" << higher << ">";
		exit(1);
	}
}

bool Options::check_if_file_is_readable(const std::string& option_name, const std::string& path){
	std::ifstream file;
	file.open(path.c_str());

	if (file.good()) {
		file.close();
		return true;
	} else {
		_log(logERR) << "Option: \"" << option_name << "\" must be accessible and readable.";
		exit(1);
	}
}

bool Options::check_if_file_is_writable(const std::string& option_name, const std::string& path){
	std::ofstream file;
	file.open(path.c_str());

	if (file.good()) {
		file.close();
		return true;
	} else {
		_log(logERR) << "Option: \"" << option_name << "\" must be accessible and writable.";
		exit(1);
	}
}

bool Options::check_if_one_of(const std::string& option_name, const std::string& value, const std::vector<std::string>& valid_options){
	BOOST_FOREACH(const std::string& valid_opt, valid_options){
		if( value == valid_opt){
			return true;
		}
	}
	_log(logERR) << "Option: \"" << option_name << "\" must be one of:" << boost::algorithm::join(valid_options, ", ");
	exit(1);
}


std::ostream& operator<<(std::ostream& os, const Options& options)
{
	os << "---------------------\n";
	os << "Options summary:\n\n";

	os << "General options:\n";
	os << "--input_file_path = " << options.point_input_file << '\n';
	os << "--output_clusters_file_path = " << options.output_file << '\n';
	os << "--output_cluster_profiles_file = " << options.output_centers_file << '\n';
	os << "--cluster_name_prefix = " << options.output_cluster_prefix << '\n';
	os << "--num_threads = " << options.num_threads << '\n';
	os << "--verbosity = " << options.verbosity_option << '\n';
	os << '\n';

	os << "Algorithm options:\n";
	os << "--max_canopy_dist = " << options.max_canopy_dist << '\n';
	os << "--max_close_dist = " << options.max_close_dist << '\n';
	os << "--max_merge_dist = " << options.max_merge_dist << '\n';
	os << "--min_step_dist = " << options.min_step_dist << '\n';
	os << "--max_num_canopy_walks = " << options.max_num_canopy_walks << '\n';
	os << '\n';

	os << "Input filtering options:\n";
	os << "--filter_min_obs = " << options.min_non_zero_data_samples << '\n';
	os << "--filter_max_dominant_obs = " << options.max_top_three_data_point_proportion << '\n';
	if (options.points_filtered_out_at_least_non_zero_file_path != "")
		os << "--filtered_out_points_min_obs_file = " << options.points_filtered_out_at_least_non_zero_file_path << '\n';
	if (options.points_filtered_out_top_three_prop_file_path != "")
		os << "--filtered_out_points_max_dominant_obs_file = " << options.points_filtered_out_top_three_prop_file_path << '\n';
	os << '\n';

	os << "Output filtering options:\n";
	os << "--filter_zero_medians = " << options.min_num_non_zero_medians << '\n';
	os << "--filter_single_point = " << options.max_single_data_point_proportion << '\n';
	os << '\n';

	os << "Early stop options:\n";
	os << "--stop_fraction = " << options.stop_proportion_of_points << '\n';
	os << '\n';

	os << "Miscellaneous:\n";
	os << "--die_on_kill = " << std::boolalpha << options.die_on_kill << '\n';
	if (options.not_processed_points_file != "")
		os << "--not_processed_points_file = " << options.not_processed_points_file << '\n';
	os << "--print_time_statistics = " << std::boolalpha << options.print_time_statistics << '\n';
	os << "--show_progress_bar = " << std::boolalpha << options.show_progress_bar << '\n';
	if (options.canopy_size_stats_file != "")
		os << "--canopy_size_stats_file = " << options.canopy_size_stats_file << '\n';
	os << "---------------------\n";

	return os;
}
