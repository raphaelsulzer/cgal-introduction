#include <cgal_typedefs.h>
#include <fileIO.h>




double pointPlaneDistance(Plane plane, Point point){

    // implemented from: http://mathworld.wolfram.com/Point-PlaneDistance.html

    double dist = (plane.a()*point.x() + plane.b()*point.y() + plane.c()*point.z() + plane.d()) /
            std::sqrt(plane.a()*plane.a() + plane.b()*plane.b() + plane.c()*plane.c());

    return dist;

}



int tetIntersectionFun(Tetrahedron& tet,
                           std::vector<Plane>& all_planes,
                           double& vol,
                           int plane_count=0,
                           int global_neg=0, int global_pos=0){


    double epsilon = 0.00001;

//    if(plane_count == 4){
//        vol = tet.volume();
//        return 1;
//    }


    Plane plane = all_planes[plane_count];

    std::vector<int> neg;
    std::vector<int> pos;
    std::vector<int> close;
    for(int v = 0; v < 4; v++){
        // test if point is on positive or negative side
        // only if there are points on both sides of this plane I need to do sth with this plane
        Point point = tet.vertex(v);
        double dist = pointPlaneDistance(plane, point);
        if(dist < -epsilon){
            neg.push_back(v);
            global_neg++;
        }
        else if(dist > epsilon){
            pos.push_back(v);
            global_pos++;
        }
        else {
            close.push_back(v);
        }
    }
    // continue according to intersection with the current plane
    int ns = neg.size();
    int ps = pos.size();
    int cs = close.size();
    if(ns != 0 && ps != 0){
        // check if a certain triangle has points in neg and pos and intersect it with this plane
        // make a new tet from the result and call this function with it, by removing the plane from the vector
        // e.g.
        // return newTetIntersectionFun(newTet, planes-1)
        if(ns != ps){
            for(int n = 0; n < ns; n++){
                for(int p = 0; p < ps; p++){
                    Point p1 = tet.vertex(neg[n]);
                    Point p2 = tet.vertex(pos[p]);
                    Line l(p1,p2);
                    CGAL::cpp11::result_of<Intersect(Line, Plane)>::type
                                  result = intersection(l, plane);
                    if(result){
                        if(const Point* intersection_point = boost::get<Point>(&*result)){
                            // make a new tet with all the negative points
                            // and replace the positive one with the new intersection point
                            std::vector<Point> newTetPoints;
                            // add positive points
                            for(int nn = 0; nn < ns; nn++)
                                newTetPoints.push_back(tet.vertex(neg[nn]));
                            // add close points
                            for(int c = 0; c < close.size(); c++)
                                newTetPoints.push_back(tet.vertex(close[c]));
                            newTetPoints.push_back(*intersection_point);
                            Tetrahedron newTet(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            tetIntersectionFun(tet, all_planes, vol, ++plane_count);
                        }
                    }
                    break;
                }
                break;
            }
        }
    }
    vol = tet.volume();
    return 1;
//    else if(ps == 0){           // if positive side is empty, means no intersection with this plane
//        if(plane_count == 3){           // if it was the last plane, (current) tet is contained in all planes
//            vol+=tet.volume();
//            return 1;}
//        else                            // else continue with next plane
//            tetIntersectionFun(tet, all_planes, vol, ++plane_count);
//    }
//    else if(ns == 0)            // if negative side is empty, means no intersection at all
//        return 0;

}





void tetIntersectionTest(){



    Point sp0(1,1,1);
    Point sp1(1,3,1);
    Point sp2(3,2,1);
    Point sp3(2,2,4);
    Point spc = CGAL::centroid(sp0,sp1,sp2,sp3);
    Tetrahedron stet(sp0,sp1,sp2,sp3);
    Polyhedron sp;
    sp.make_tetrahedron(sp0,sp1,sp2,sp3);
    exportOFF(sp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/sp");
    std::vector<Plane> planes(4);
    planes[0] = Plane(sp0,sp2,sp1);
    planes[1] = Plane(sp0,sp1,sp3);
    planes[2] = Plane(sp1,sp2,sp3);
    planes[3] = Plane(sp0,sp3,sp2);
    for(int i = 0; i<4; i++){
        if(!planes[i].has_on_negative_side(spc)){
            planes[i]=planes[i].opposite();
        }
    }

    Point dp0(3,1,5);
    Point dp1(3,3,5);
    Point dp2(1,2,5);
    Point dp3(2,2,2);
    Tetrahedron dtet(dp0,dp1,dp2,dp3);
    Polyhedron dp;
    dp.make_tetrahedron(dp0,dp1,dp2,dp3);
    exportOFF(dp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/dp");

    double cvol = 0.0;
    tetIntersectionFun(dtet, planes, cvol);

    Polyhedron_Exact sp_exact;
    CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> sensor_modifier(sp);
    sp_exact.delegate(sensor_modifier);
    Nef_polyhedron snef(sp_exact);

    Polyhedron_Exact dp_exact;
    CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> sensor_modifier(dp);
    dp_exact.delegate(sensor_modifier);
    Nef_polyhedron dnef(dp_exact);

    Nef_polyhedron fullnef = dnef*snef;

    Polyhedron_Exact intersection_tet_exact;
    fullnef.convert_to_polyhedron(intersection_tet_exact);
    double vol_full1 = CGAL::Polygon_mesh_processing::volume(intersection_tet_exact);

    int a=5;


//    int vlen = rand() % 8 + 4;
//    std::vector<Point> points(vlen);
//    for(int i = 0; i < vlen; i++){

//        Point p(rand(), rand(), rand());
//        points[i] = p;

//    }
//    Polyhedron Poly;
//    CGAL::convex_hull_3(points.begin(), points.end(),Poly);
}
