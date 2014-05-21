#ifndef CANOPY 
#define CANOPY

#include <string>

#include <Point.hpp>

using namespace std;

/**
 * Represents a canopy
 */
class Canopy {
    public:
        //Constructor - assigns the neighour points and creates new point for center (representing canopy profile) 
        Canopy(std::vector<Point*> neighbours);

        //Destructor - deletes only the center point - not the neighours
        virtual ~Canopy();

        //Set's the profile for the center point
        void find_and_set_center();
        
        //Center point representing the canopy profile
        Point* center;

        //List of points belonging to the canopy
        std::vector<Point*> neighbours;

        //Debugging printout of the canopy
        friend std::ostream& operator<<(std::ostream& ost, const Canopy& c);
        
        //Comparison operator of Canopy objects by the number of their neighbours
        friend bool compare_canopy_ptrs_by_canopy_size(const Canopy* a, const Canopy* b);
};

bool compare_canopy_ptrs_by_canopy_size(const Canopy* a, const Canopy* b);


#endif
