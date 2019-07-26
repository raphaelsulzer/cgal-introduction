#include <cgal_typedefs.h>
#include <fileIO.h>




double pointPlaneDistance(Plane plane, Point point){

    // implemented from: http://mathworld.wolfram.com/Point-PlaneDistance.html

    double dist = plane.a()*point.x() + plane.b()*point.y() + plane.c()*point.z() + plane.d() /
            std::sqrt(plane.a()*plane.a() + plane.b()*plane.b() + plane.c()*plane.c());

    return dist;

}



void newTetIntersectionFun(Tetrahedron tet, std::vector<Plane> planes){


    double epsilon = 0.00001;

    for(int p = 0; p < planes.size(); p++){

        std::vector<int> neg;
        std::vector<int> pos;
        for(int v = 0; v < 4; v++){
            // test if point is on positive or negative side
            // only if there are points on both sides of this plane I need to do sth with this plane
            double dist = pointPlaneDistance(planes[p], tet.vertex[v]);
            if(dist < -epsilon)
                neg.push_back(v);
            else if(dist > epsilon)
                pos.push_back(v);
        }
        if(!neg.size() == 0 && !pos.size() == 0){
            // check if a certain triangle has points in neg and pos and intersect it with this plane
            // make a new tet from the result and call this function with it, by removing the plane from the vector
            // e.g.
            // return newTetIntersectionFun(newTet, planes-1)
            for(int t = 0; t < 4; t++){
                tet.triangle[t]
                // check if a certain triangle has points in neg and pos and intersect it with this plane
                        // is hopefully doable with the tet based vertex index


            }
        }
        else
            continue; //to next plane
    }
}







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
