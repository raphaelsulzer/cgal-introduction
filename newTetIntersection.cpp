#include <cgal_typedefs.h>
#include <fileIO.h>

void newTetIntersectionTest(){

    int vlen = rand() % 8 + 4;
    std::vector<Point> points(vlen);
    for(int i = 0; i < vlen; i++){

        Point p(rand(), rand(), rand());
        points[i] = p;

    }

    Polyhedron Poly;
    CGAL::convex_hull_3(points.begin(), points.end(),Poly);
}
