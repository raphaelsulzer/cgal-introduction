//
//  exportOFF.cpp
//  cgal1
//
//  Created by Raphael Sulzer on 30/01/2019.
//  Copyright © 2019 Raphael Sulzer. All rights reserved.
//

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K>      Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;
typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef Kernel::Facet_3                              Facet_3;
typedef CGAL::Polyhedron_3<Kernel>                   Polyhedron;
typedef Triangulation::Facet_iterator                   Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

int main() {
    
    // Polyhedron
//    Point_3 p( 0.0, 0.0, 0.0);
//    Point_3 q( 1.0, 0.0, 0.0);
//    Point_3 r( 0.0, 1.0, 0.0);
//    Point_3 s( 0.0, 0.0, 1.0);
//
//    Polyhedron P;
//    P.make_tetrahedron( p, q, r, s);
    
    
    // Tetrahedron
    std::list<Point> L;
    L.push_front(Point(0,0,0));
    L.push_front(Point(1,0,0));
    L.push_front(Point(0,1,0));
    L.push_front(Point(0,0,1));
    Triangulation T(L.begin(), L.end());
    Triangulation::size_type nv = T.number_of_vertices();
    Triangulation::size_type nf = T.number_of_facets();

    
    
    
    
    
    // Write polyhedron in Object File Format (OFF).
    CGAL::set_ascii_mode( std::cout);
    
    
    std::ofstream oFileT("/Users/Raphael/Library/Mobile Documents/com~apple~CloudDocs/Studium/PhD/Paris/data/triangulation2.off",std::ios::out);

    
    oFileT << "OFF" << std::endl << nv << ' '
    << nf << " 0" << std::endl;
    std::copy( T.points_begin(), T.points_end(),
              std::ostream_iterator<Point_3>( oFileT, "\n"));
    
    std::copy( T.facets_begin(), T.facets_end(),
              std::ostream_iterator<Facet_3>( oFileT, "\n"));
    
    
    for (  Facet_iterator i = T.Finite_facets_begin(); i != T.Finite_facets_end(); ++i) {
        Halfedge_facet_circulator j = i->Finite_facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        oFileT << CGAL::circulator_size(j) << ' ';
        do {
            oFileT << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        oFileT << std::endl;
    }
    
    

    
    
    return 0;
}
