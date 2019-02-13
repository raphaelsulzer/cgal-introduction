//
//  triangulation3.cpp
//  cgal1
//
//  Created by Raphael Sulzer on 31/01/2019.
//  Copyright © 2019 Raphael Sulzer. All rights reserved.
//

#include <stdio.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
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
int main()
{
    // construction from a list of points :
    std::list<Point> L;
    L.push_front(Point(0,0,0));
    L.push_front(Point(1,0,0));
    L.push_front(Point(0,1,0));
    L.push_front(Point(0,0,1));
    Triangulation T(L.begin(), L.end());
    Triangulation::size_type nv = T.number_of_vertices();
    Triangulation::size_type nf = T.number_of_facets();
    

    std::cout << nv << std::endl;
    std::cout << nf << std::endl;
    
    

}
