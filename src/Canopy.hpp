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
};


#endif
