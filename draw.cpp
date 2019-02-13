//
//  draw.cpp
//  cgal1
//
//  Created by Raphael Sulzer on 30/01/2019.
//  Copyright © 2019 Raphael Sulzer. All rights reserved.
//

#include <stdio.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/draw_triangulation_3.h>

#define CGAL_USE_BASIC_VIEWER

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K>                   DT3;
typedef CGAL::Creator_uniform_3<double,K::Point_3>          Creator;
int main()
{
    std::vector<K::Point_3> points;
    CGAL::Random_points_in_sphere_3<K::Point_3,Creator> g(1.0);
    CGAL::cpp11::copy_n(g, 50, std::back_inserter(points));
    DT3 dt3(points.begin(), points.end());
    CGAL::draw(dt3);
    return EXIT_SUCCESS;
}
