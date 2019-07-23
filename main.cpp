#include <cgal_typedefs.h>
#include <fileIO.h>

#include "surfaceRecon.cpp"


int main(){

    surfaceReconstruction();

//    using namespace Eigen;
//    Vector3d v0, v1;
//    v0 << 51267.609375, 91991.4375, 139.50198364257812;
//    v1 << 51267.59765625, 91991.4375, 139.50042724609375;
//    Vector3d edge1;
//    edge1 = v1 - v0;

//    Point a1(1,1,0);
//    Point a2(1,3,0);
//    Point a3(3,2,0);
//    Point a4(2,2,2);

//    Point b1(3,1,5);
//    Point b2(3,3,5);
//    Point b3(1,2,5);
//    Point b4(2,2,1);

//    Point c1(1,0,1);
//    Point c2(1,2,1);
//    Point c3(3,1,1);
//    Point c4(2,1,3);

//    Polyhedron a;
//    a.make_tetrahedron(a1,a2,a3,a4);

//    Polyhedron b;
//    b.make_tetrahedron(b1,b2,b3,b4);

//    Polyhedron c;
//    c.make_tetrahedron(c1,c2,c3,c4);

//    tetIntersectionFun(a,b);

    return 0;

}
