#include <vector>
#include <boost/unordered_set.hpp>

#include <Point.hpp>
#include <Canopy.hpp>

using namespace std;

class CanopyClusteringAlg{
    public:
        static std::vector<Canopy*> single_core_run_clustering_on(std::vector<Point*>& points);
        static std::vector<Canopy*> multi_core_run_clustering_on(std::vector<Point*>& points);
        static Canopy* create_canopy(Point* origin, vector<Point*>& points, vector<Point*>& close_points, double min_neighbour_correlation, double min_close_correlation, bool sets_close_points);

        static void filter_clusters_by_zero_medians(int min_num_non_zero_medians, std::vector<Canopy*>& canopies_to_filter);
        static void filter_clusters_by_single_point_skew(double max_single_data_point_proportion, std::vector<Canopy*>& canopies_to_filter);
};
