#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

struct cell_info{
    int idx;
    float outside_score;
    float inside_score;
    int final_label;
};

struct vertex_info{
    int idx;
    float outside_score;
    float inside_score;
    int final_label;
};


///////// FILE I/O /////////
#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
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
typedef std::tuple<int, double, Vector, Color>                      IdxSigNormCol;
typedef CGAL::Triangulation_vertex_base_with_info_3<IdxSigNormCol, Kernel>     VB;
//// vertex base for point + info (=vector or camera)
//typedef CGAL::Triangulation_vertex_base_with_info_3<Vector, Kernel>         VerNormB;
//typedef CGAL::Triangulation_vertex_base_with_info_3<int, Kernel>            VerCamB;



typedef CGAL::Triangulation_cell_base_with_info_3<cell_info, Kernel>           CB;         // cell base
typedef CGAL::Triangulation_data_structure_3<VB, CB>                Tds;        // triangulation data structure
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>                 Delaunay;   // delaunay triangulation based on triangulation data structure
typedef Delaunay::Point                                             Point;
typedef Delaunay::Edge                                              Edge;
typedef Delaunay::Facet                                             Facet;
typedef Delaunay::Cell_handle                                       Cell_handle;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
// map cell of the Dt, to an index, the outside score, the inside score, and the final label
typedef std::map<Cell_handle, std::tuple<int, float, float, int>>   Cell_map;   // not really needed anymore, just for old code
typedef std::map<Vertex_handle, int>                                Vertex_map;

///////// read PLY /////////
typedef CGAL::cpp11::tuple<Point, Vector, Color> PNC;

///////// ray tracing /////////
typedef Kernel::Ray_3                                               Ray;
typedef Kernel::Triangle_3                                          Triangle;
typedef Kernel::Intersect_3                                         Intersect;
typedef Kernel::Segment_3                                           Segment;

///////// Polyhedron mesh /////////
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/tags.h>

typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>      Polyhedron;


#endif // CGAL_TYPEDEFS_H



