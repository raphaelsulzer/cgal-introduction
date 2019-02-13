#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_ply_points.h>
#include <utility>
#include <vector>
#include <fstream>
// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::cpp11::array<unsigned char, 3> Color;
// Point with normal, color and intensity
typedef CGAL::cpp11::tuple<Point, Vector, Color, int> PNCI;
typedef CGAL::Nth_of_tuple_property_map<0, PNCI> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNCI> Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNCI> Color_map;
typedef CGAL::Nth_of_tuple_property_map<3, PNCI> Intensity_map;


std::list<Point> readPLYfun(const char* fname)
{
    // Reads a .ply point set file with normal vectors and colors
    std::list<Point> points; // store points
    std::ifstream in(fname);
    CGAL::read_ply_points(in, std::back_inserter (points));

    // Display points read
    for (auto pt = points.begin(); pt != points.end(); pt++)
    {
        const Point& p = *pt;
        std::cout << p << std::endl;
    }
    
    return points;
}
