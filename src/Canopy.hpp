#ifndef CANOPY 
#define CANOPY

#include <string>

#include <Point.hpp>

using namespace std;

class Canopy {
    public:
        Canopy();
        Canopy(Point* origin, Point* center, std::vector<Point*> neighbours);
        

        Point* origin;
        Point* center;
        std::vector<Point*> neighbours;

        friend std::ostream& operator<<(std::ostream& ost, const Canopy& c);
        friend bool compare_canopy_ptrs(const Canopy* a, const Canopy* b);
};

bool compare_canopy_ptrs(const Canopy* a, const Canopy* b);


#endif
