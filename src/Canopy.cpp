
#include <boost/foreach.hpp>

#include <Canopy.hpp>

Canopy::Canopy(Point* origin, Point* center, std::vector<Point*>& neighbours): origin(origin), center(center), neighbours(neighbours){}

std::ostream& operator<<(std::ostream& ost, const Canopy& c)
{
    ost << ">>>>>>>>>>Canopy>>>>>>>>" << std::endl;
    ost << "Origin:" << std::endl;
    ost << *c.origin;
    ost << "Center:" << std::endl;
    ost << *c.center;
    ost << "Neighbours:" << std::endl;
    BOOST_FOREACH(const Point* p, c.neighbours)
        ost << p->id << "\t";
    ost << std::endl;
    ost << ">>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;


}
