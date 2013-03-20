
#include <boost/foreach.hpp>

#include <Canopy.hpp>

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
    ost << "Neighbours:" << std::endl;
    BOOST_FOREACH(const Point* p, c.neighbours)
        ost << p->id << "\t";
    ost << std::endl;
    ost << ">>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;


}

bool compare_canopy_ptrs(const Canopy* a, const Canopy* b){
    return (a->neighbours.size() > b->neighbours.size());
}


