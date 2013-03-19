#ifndef CANOPY_CLUSTERING
#define CANOPY_CLUSTERING
#include <vector>
#include <boost/unordered_set.hpp>

#include <Point.hpp>
#include <Canopy.hpp>
#include <TimeProfile.hpp>

using namespace std;

class CanopyClusteringAlg{
    public:
        static std::vector<Canopy*> single_core_run_clustering_on(vector<Point*>& points);
        static std::vector<Canopy*> multi_core_run_clustering_on(vector<Point*>& points, int num_threads, double max_canopy_dist, double max_close_dist, double max_merge_dist, double max_step_dist, double stop_proportion_of_points, int stop_num_single_point_clusters, string canopy_size_stats_fp, bool show_progress_bar, TimeProfile& time_profile);
        static Canopy* create_canopy(Point* origin, vector<Point*>& points, vector<Point*>& close_points, double min_neighbour_correlation, double min_close_correlation, bool sets_close_points);

        static void filter_clusters_by_zero_medians(int min_num_non_zero_medians, vector<Canopy*>& canopies_to_filter);
        static void filter_clusters_by_single_point_skew(double max_single_data_point_proportion, vector<Canopy*>& canopies_to_filter);
};

#endif 
