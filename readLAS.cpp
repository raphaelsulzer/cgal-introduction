//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/property_map.h>
//#include <CGAL/IO/read_las_points.h>
//#include <utility>
//#include <vector>
//#include <fstream>
//// types
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef Kernel::FT FT;
//typedef Kernel::Point_3 Point;
//typedef CGAL::cpp11::array<unsigned short, 4> Color;
//typedef std::pair<Point, Color> PointWithColor;
//int main(int argc, char*argv[])
//{
//    const char* fname = (argc>1) ? argv[1] : "/Users/Raphael/Library/Mobile Documents/com~apple~CloudDocs/Studium/PhD/Paris/data/bouwpub.ply";
//    // Reads a .las point set file with normal vectors and colors
//    std::vector<Point> points; // store points
//    std::ifstream in(fname, std::ios_base::binary);
//    if (!in ||
//        !CGAL::read_las_points
//        (in,
//         std::back_inserter (points)))
//    {
//        std::cerr << "Error: cannot read file " << fname << std::endl;
//        return EXIT_FAILURE;
//    }
//    for (std::size_t i = 0; i < points.size(); ++ i)
//        std::cout << points[i] << std::endl;
//
//    return EXIT_SUCCESS;
//}
