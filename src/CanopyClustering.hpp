#include <vector>
#include <boost/unordered_set.hpp>

#include <Point.hpp>

struct Canopy{
    Point center;
    std::vector<Point> neighbours;
};

class CanopyClusteringAlg{
    public:
        static std::vector<Canopy> single_core_run_clustering_on(std::vector<Point>& points);
        static Canopy create_canopy(const Point& center, boost::unordered_set<Point>& marked_points, const std::vector<Point>& points);
};
