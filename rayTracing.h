#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <cgal_typedefs.h>

// TODO: check if I can replace the whole thing with the nearest_vertex(const Point& p, Cell_handle start) function,
// which can be found in the Delaunay_triangulation_3.h file in /usr/lib/CGAL

namespace rayTracing{


////// cross product of two vector array.
//std::vector<double> crossProduct(std::vector<double> a, std::vector<double> b){
//    std::vector<double> cp(3);
//    cp[0] = a[1] * b[2] - a[2] * b[1];
//    cp[1] = a[0] * b[2] - a[2] * b[0];
//    cp[2] = a[0] * b[1] - a[1] * b[0];

//    return cp;
//}
//double dotProduct(std::vector<double> a, std::vector<double> b){
//    return  a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
//}

// TODO: maybe update ray triangle intersection with this:
// https://stackoverflow.com/questions/44275153/is-m%C3%B6ller-trumbore-ray-intersection-the-fastest
bool rayTriangleIntersection(Point& rayOrigin,
                           Vector& rayVector,
                           Triangle& inTriangle,
                           Point& outIntersectionPoint);

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
std::pair<float, float> cellScore(float dist2, double eig3, bool inside);


// TODO: why can the Delaunay be const here? I'm changing the cell scores that are saved inside the Delaunay!
int traverseCells(Delaunay& Dt,
                  Cell_handle& current_cell, std::unordered_set<Cell_handle>& processed,
                  Ray ray, double sigma, int oppositeVertex,
                  bool inside);

void firstCell(Delaunay& Dt, Delaunay::Finite_vertices_iterator& vit,
               bool inside,
               int& intersection_count);

void rayTracingFun(Delaunay& Dt);

// end of namespace rayTracing
}










#endif // RAYTRACING_H
