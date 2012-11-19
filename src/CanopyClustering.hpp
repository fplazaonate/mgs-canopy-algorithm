#include <vector>
#include <boost/unordered_set.hpp>

#include <Point.hpp>

struct Canopy{
    Point* center;
    Point* origin;
    std::vector<Point*> neighbours;
};

std::ostream& operator<<(std::ostream& ost, const Canopy& c);

class CanopyClusteringAlg{
    public:
        static std::vector<Canopy> single_core_run_clustering_on(std::vector<Point*>& points);
        static Canopy create_canopy(Point* origin, boost::unordered_set<Point*>& marked_points, std::vector<Point*>& points, double min_correlation);
};
