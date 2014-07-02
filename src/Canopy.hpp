/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
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
        //Constructor - copies the given point as center, neighbours are empty (used only in rare cases) 
        Canopy(Point* center_to_copy);

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
