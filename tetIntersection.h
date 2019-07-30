#ifndef TETINTERSECTION_H
#define TETINTERSECTION_H

#include <cgal_typedefs.h>
#include <fileIO.h>




double pointPlaneDistance(Plane plane, Point point);


int tetIntersectionFun(Tetrahedron& tet,
                            std::vector<Plane>& all_planes,
                            double& vol,
                            int plane_count,
                            std::string tet_name,
                            bool exp);

void tetIntersectionTest();







#endif // TETINTERSECTION_H
