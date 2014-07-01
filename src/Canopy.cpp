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


#include <boost/foreach.hpp>

#include <Canopy.hpp>

Canopy::Canopy(Point* center_to_copy){
    center = new Point(*center_to_copy);
    center->id = "!GENERATED!";
}

Canopy::Canopy(std::vector<Point*> neighbours): neighbours(neighbours){
    find_and_set_center();
}

Canopy::~Canopy(){
    delete center;
}

void Canopy::find_and_set_center(){

    center = get_centroid_of_points(neighbours);

}

std::ostream& operator<<(std::ostream& ost, const Canopy& c)
{
    ost << ">>>>>>>>>>Canopy>>>>>>>>" << std::endl;
    ost << "Center:" << std::endl;
    if(c.center != NULL)
        ost << *c.center;
    else
        ost << "===NONE===" << endl;
    ost << "Neighbours: " << c.neighbours.size() << std::endl;
    //BOOST_FOREACH(const Point* p, c.neighbours)
    //    ost << p->id << "\t";
    //ost << std::endl;
    ost << ">>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;


}

bool compare_canopy_ptrs_by_canopy_size(const Canopy* a, const Canopy* b){
    return (a->neighbours.size() > b->neighbours.size());
}


