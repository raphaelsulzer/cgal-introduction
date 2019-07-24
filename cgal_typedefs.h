#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H


///////// FILE I/O /////////
#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
///////// tetTracing /////////
#include <unordered_set>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Euclidean_distance.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel           EPECK;

// vertext base for point + info (=vector, color, intensity)
typedef EPICK::Vector_3                                            Vector;
typedef EPICK::Point_3                                             Point;
typedef CGAL::cpp11::array<unsigned char, 3>                       Color;
typedef CGAL::cpp11::array<float, 3>                               Sensor;

struct vertex_info{
    int idx = 0;
    double sigma = 0.0;
    Color color;
    Vector normal;
    Vector sensor_vec;
    Point sensor_pos;
    std::vector<int> sensor_tet;
};
typedef CGAL::Triangulation_vertex_base_with_info_3<vertex_info, EPICK>    VB;

struct cell_info{
    int idx = 0;
    float outside_score = 0.0;
    float inside_score = 0.0;
    int final_label = 0;
};
typedef CGAL::Triangulation_cell_base_with_info_3<cell_info, EPICK>        CB;         // cell base

// Delaunay triangulation data structure
typedef CGAL::Triangulation_data_structure_3<VB, CB>                Tds;        // triangulation data structure
typedef CGAL::Delaunay_triangulation_3<EPICK, Tds>                 Delaunay;   // delaunay triangulation based on triangulation data structure
//typedef Delaunay::Point                                             Point;
typedef Delaunay::Edge                                              Edge;
typedef Delaunay::Facet                                             Facet;
typedef Delaunay::Cell_handle                                       Cell_handle;
typedef Delaunay::Vertex_handle                                     Vertex_handle;

typedef boost::tuple<Point, vertex_info>                            point_info;

// map cell of the Dt, to an index, the outside score, the inside score, and the final label
// not really needed anymore, just for old code
typedef std::map<Cell_handle, std::tuple<int, float, float, int>>   Cell_map;
typedef std::map<Vertex_handle, int>                                Vertex_map;
typedef std::map<Vertex_handle, std::pair<Point,double>>            VPS_map;

///////// read PLY /////////
typedef CGAL::cpp11::tuple<Point, Vector, Color> PNC;
typedef CGAL::cpp11::tuple<Point, Vector> PN;
typedef CGAL::cpp11::tuple<Point, Sensor> PS;

///////// ray tracing /////////
typedef EPICK::Ray_3                                               Ray;
typedef EPICK::Triangle_3                                          Triangle;
typedef EPICK::Intersect_3                                         Intersect;
typedef EPICK::Segment_3                                           Segment;
typedef EPICK::Tetrahedron_3                                       Tetrahedron;
typedef CGAL::Cartesian_converter<EPICK,EPECK>                     IK_to_EK;
typedef CGAL::Cartesian_converter<EPECK,EPICK>                     EK_to_IK;

///////// Polyhedron mesh /////////
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/tags.h>

typedef CGAL::Polyhedron_3<EPICK, CGAL::Polyhedron_items_with_id_3>      Polyhedron;
//typedef CGAL::Polyhedron_3<EPICK, CGAL::Polyhedron_items_3>      Polyhedron;

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>
namespace SMS = CGAL::Surface_mesh_simplification;



// for PCA / neighborhood search
// for neighborhood search
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Incremental_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/centroid.h>
#include <CGAL/estimate_scale.h>

// for matrix operations, use Eigen lib
#include <Eigen/Core>
#include <Eigen/Dense>
#include<Eigen/Geometry>


//#include <pcl/point_types.h>
//#include <pcl/features/normal_3d.h>


typedef CGAL::Search_traits_3<EPICK> TreeTraits;
typedef CGAL::Search_traits_adapter<point_info,
  CGAL::Nth_of_tuple_property_map<0, point_info>,
  TreeTraits>                                               Traits;

typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

typedef CGAL::Orthogonal_incremental_neighbor_search<Traits> Incremental_neighbor_search;
typedef Incremental_neighbor_search::Tree Incremental_Tree;

////AABB tree
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<EPICK, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> AABB_Tree;
typedef boost::optional< AABB_Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;
typedef boost::optional< AABB_Tree::Intersection_and_primitive_id<Point>::Type > Point_intersection;
typedef AABB_Tree::Primitive_id Primitive_id;


//// Tetrahedron intersection
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_copy_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Tetrahedron_3.h>

typedef EPICK::Plane_3 Plane;

typedef EPECK::Point_3 Point_Exact;
typedef CGAL::Polyhedron_3<EPECK, CGAL::Polyhedron_items_with_id_3> Polyhedron_Exact;
//typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_3> Polyhedron_Exact;
typedef CGAL::Nef_polyhedron_3<EPECK> Nef_polyhedron;


////#include <CGAL/Nef_polyhedron_3.h>
////#include <CGAL/Extended_homogeneous.h>
////#include <CGAL/Exact_integer.h>
////#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>


//////typedef CGAL::Nef_polyhedron_3<CGAL::Extended_homogeneous<CGAL::Exact_integer>> Nef;
////typedef CGAL::Nef_polyhedron_3<EPICK> Nef;
//////typedef Nef::Plane_3  Nef_Plane;
////typedef EPICK::Plane_3 Plane;
#endif // CGAL_TYPEDEFS_H



