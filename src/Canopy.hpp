#ifndef CANOPY 
#define CANOPY

#include <string>

#include <Point.hpp>

using namespace std;

class Canopy {
    public:
        Canopy(std::vector<Point*> neighbours);
        virtual ~Canopy();

        void find_and_set_center();
        
        Point* center;
        std::vector<Point*> neighbours;

        friend std::ostream& operator<<(std::ostream& ost, const Canopy& c);
        friend bool compare_canopy_ptrs(const Canopy* a, const Canopy* b);
};

bool compare_canopy_ptrs(const Canopy* a, const Canopy* b);


#endif
