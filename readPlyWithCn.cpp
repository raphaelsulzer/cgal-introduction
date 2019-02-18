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

std::vector<PNCI> readPlyWithCnFun(const char* fname)
{

    std::vector<PNCI> points; // store points
    std::ifstream in(fname);

    CGAL::read_ply_points_with_properties
      (in,
       std::back_inserter (points),
       CGAL::make_ply_point_reader (Point_map()),
       std::make_pair (Intensity_map(),
                       CGAL::PLY_property<int>("intensity")),
       std::make_tuple (Color_map(),
                        CGAL::Construct_array(),
                        CGAL::PLY_property<unsigned char>("red"),
                        CGAL::PLY_property<unsigned char>("green"),
                        CGAL::PLY_property<unsigned char>("blue")),
       CGAL::make_ply_normal_reader (Normal_map()));

    // Display points read
    for (std::size_t i = 0; i < points.size (); ++ i)
    {
      const Point& p = get<0>(points[i]);
      const Vector& n = get<1>(points[i]);
      const Color& c = get<2>(points[i]);
      int I = get<3>(points[i]);
      std::cerr << "Point (" << p << ") with normal (" << n
                << "), color (" << int(c[0]) << " " << int(c[1]) << " " << int(c[2])
                << ") and intensity " << I << std::endl;
    }

    return points;
}
