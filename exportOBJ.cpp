//
//  exportOBJ.cpp
//  cgal1
//
//  Created by Raphael Sulzer on 27/01/2019.
//  Copyright © 2019 Raphael Sulzer. All rights reserved.
//

#include <stdio.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K>      Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;
typedef Polyhedron::Facet_iterator    Facet_iterator;


// write in obj file format (OBJ).
void write_obj(char *pFilename, std::list<Point> pts,
               int incr  = 1) // 1-based by default
{
    std::ofstream stream(pFilename);
    
    
    // output vertices
    for(std::list<Point>::iterator Point_iterator = pts.begin();
        Point_iterator != pts.end();
        Point_iterator++)
        stream << 'v' << ' ' <<
        Point_iterator->x() << ' ' <<
        Point_iterator->y() << ' ' <<
        Point_iterator->z() << std::endl;
    
    // precompute vertex indices
    this->set_index_vertices();
    
    // this assumes you have defined an index per vertex
    // and the function:
    
    void set_index_vertices()
    {
        unsigned int index = 0;
        for(Vertex_iterator v = vertices_begin();
            v!= vertices_end();
            v++)
            v->index(index++);
    }
    // or something with a STL map to associate one index per vertex
    
    
    // output facets
    for(Facet_iterator pFacet = facets_begin();
        pFacet != facets_end();
        pFacet++)
    {
        stream << 'f';
        Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
        do
            stream << ' ' << pHalfedge->vertex()->index()+incr;
        while(++pHalfedge != pFacet->facet_begin());
        stream << std::endl;
    }
}


int main()
{
    // construction from a list of points :
    std::list<Point> L;
    L.push_front(Point(0,0,0));
    L.push_front(Point(1,0,0));
    L.push_front(Point(0,1,0));
    Triangulation T(L.begin(), L.end());
    Triangulation::size_type n = T.number_of_vertices();
    // insertion from a vector :
    std::vector<Point> V(3);
    V[0] = Point(0,0,1);
    V[1] = Point(1,1,1);
    V[2] = Point(2,2,2);
    n = n + T.insert(V.begin(), V.end());
    assert( n == 6 );       // 6 points have been inserted
    assert( T.is_valid() ); // checking validity of T
    Locate_type lt;
    int li, lj;
    Point p(0,0,0);
    Cell_handle c = T.locate(p, lt, li, lj);
    // p is the vertex of c of index li :
    assert( lt == Triangulation::VERTEX );
    assert( c->vertex(li)->point() == p );
    Vertex_handle v = c->vertex( (li+1)&3 );
    // v is another vertex of c
    Cell_handle nc = c->neighbor(li);
    // nc = neighbor of c opposite to the vertex associated with p
    // nc must have vertex v :
    int nli;
    assert( nc->has_vertex( v, nli ) );
    // nli is the index of v in nc
    std::ofstream oFileT("/Users/Raphael/Library/Mobile Documents/com~apple~CloudDocs/Studium/PhD/Paris/data/triangulation1.off",std::ios::out);
    // writing file output;
    oFileT << T;
    Triangulation T1;
    std::ifstream iFileT("/Users/Raphael/Library/Mobile Documents/com~apple~CloudDocs/Studium/PhD/Paris/data/triangulation1.off",std::ios::in);
    // reading file output;
    iFileT >> T1;
    assert( T1.is_valid() );
    assert( T1.number_of_vertices() == T.number_of_vertices() );
    assert( T1.number_of_cells() == T.number_of_cells() );
    return 0;
    
}
