
#include <boost/foreach.hpp>

#include <Canopy.hpp>

Canopy::Canopy(){
    center = NULL;
    origin = NULL;
}

Canopy::Canopy(Point* origin, Point* center, std::vector<Point*>& neighbours): origin(origin), center(center), neighbours(neighbours){}

std::ostream& operator<<(std::ostream& ost, const Canopy& c)
{
    ost << ">>>>>>>>>>Canopy>>>>>>>>" << std::endl;
    ost << "Origin:" << std::endl;
    if(c.origin != NULL)
        ost << *c.origin;
    else
        ost << "===NONE===" << endl;
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
