#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

#endif // CGAL_TYPEDEFS_H

///////// FILE I/O /////////
#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Euclidean_distance.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;

// vertext base for point + info (=vector, color, intensity)
typedef Kernel::Vector_3                                            Vector;
typedef CGAL::cpp11::array<unsigned char, 3>                        Color;
//typedef std::tuple<Vector, Color>                                   VColb;
// vertext base for point + info (=vector)
typedef CGAL::Triangulation_vertex_base_with_info_3<Vector, Kernel> VNb;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, Kernel> VCb;


typedef CGAL::Delaunay_triangulation_cell_base_3<Kernel>            Cb;         // cell base
typedef CGAL::Triangulation_data_structure_3<VNb, Cb>               Tds;        // triangulation data structure
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>                 Delaunay;   // delaunay triangulation based on triangulation data structure
typedef Delaunay::Point                                             Point;
typedef Delaunay::Edge                                              Edge;
typedef Delaunay::Facet                                             Facet;
typedef Delaunay::Cell_handle                                       Cell_handle;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
// map cell of the Dt, to an index, the outside score, the inside score, and the final label
typedef std::map<Cell_handle, std::tuple<int, float, float, int>>   Cell_map;
typedef std::map<Vertex_handle, int>                                Vertex_map;
typedef std::map<Vertex_handle, std::pair<Point,double>>            VNC_map;

// for reading PLY file
typedef CGAL::cpp11::tuple<Point, int> PC;
typedef CGAL::Nth_of_tuple_property_map<0, PC> PointC_map;
typedef CGAL::Nth_of_tuple_property_map<1, PC> Camera_map;

typedef CGAL::cpp11::tuple<Point, Vector> PN;
typedef CGAL::Nth_of_tuple_property_map<0, PN> PointN_map;
typedef CGAL::Nth_of_tuple_property_map<1, PN> Normal_map;


///////// ray tracing /////////
typedef Kernel::Ray_3                                               Ray;
typedef Kernel::Triangle_3                                          Triangle;
typedef Kernel::Intersect_3                                         Intersect;
typedef Kernel::Segment_3                                           Segment;

