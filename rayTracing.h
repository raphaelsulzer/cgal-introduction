#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <cgal_typedefs.h>

namespace rayTracing{


bool rayTriangleIntersection(Point& rayOrigin,
                           Vector& rayVector,
                           Triangle& inTriangle,
                           Point& outIntersectionPoint);

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
std::pair<double, double> cellScore(double dist2, double eig3, bool inside);


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
