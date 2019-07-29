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
                            std::string tet_name = ""){


    std::cout << plane_count << std::endl;

    if(plane_count == 4){

        double newVol = abs(tet.volume());
        std::cout << "vol counted tet: " << tet_name << std::endl;
        vol+=newVol;
        return 1;
    }


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
        }
        else if(dist > epsilon){
            pos.push_back(v);
        }
        else {
            close.push_back(v);
        }
    }
    // continue according to intersection with the current plane
    int ns = neg.size();
    int ps = pos.size();
    int cs = close.size();
    Tetrahedron newTet;
    if(ns != 0 && ps != 0){
        // check if a certain triangle has points in neg and pos and intersect it with this plane
        // make a new tet from the result and call this function with it, by removing the plane from the vector
        // e.g.
        // return newTetIntersectionFun(newTet, planes-1)
        if(ps > ns){
            std::set<Point> newTetPoints;
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
                            // add positive points
                            for(int nn = 0; nn < ns; nn++)
                                newTetPoints.insert(tet.vertex(neg[nn]));
                            // add close points
                            for(int c = 0; c < close.size(); c++)
                                newTetPoints.insert(tet.vertex(close[c]));
                            newTetPoints.insert(*intersection_point);
                            if(newTetPoints.size() == 4){
                                auto newTetPointIt = newTetPoints.begin();
                                Point t1 = *newTetPointIt++;
                                Point t2 = *newTetPointIt++;
                                Point t3 = *newTetPointIt++;
                                Point t4 = *newTetPointIt;
                                newTet = Tetrahedron(t1,t2,t3,t4);
//                                Tetrahedron newTet(t1,t2,t3,t4);
                                Polyhedron sp;
                                sp.make_tetrahedron(t1,t2,t3,t4);
                                std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/tet_";
                                std::string tet_name = std::to_string(plane_count)+"_"+std::to_string(n)+std::to_string(p);
                                exportOFF(sp, path+tet_name);
                                tet = newTet;
                                tetIntersectionFun(tet, all_planes, vol, ++plane_count, tet_name);
                            }

                        }
                    }
                }
            }
        }
        else if(ns > ps)
//        if(ns > ps)
        {
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
//                            Tetrahedron newTet(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            newTet = Tetrahedron(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            Polyhedron sp;
                            sp.make_tetrahedron(newTet.vertex(0), newTet.vertex(1), newTet.vertex(2), newTet.vertex(3));
                            std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/tet_";
                            std::string tet_name = std::to_string(plane_count)+"_"+std::to_string(n)+std::to_string(p);
                            exportOFF(sp, path+tet_name);
                            tet = newTet;
                            tetIntersectionFun(tet, all_planes, vol, ++plane_count, tet_name);
                        }
                    }
                }
            }
        }
        else{
//        if(ns==ps){
            for(int n = 0; n < ns; n++){
                for(int p = 0; p < ps; p++){
                    if(n<p){
                        continue;
                    }
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
//                            Tetrahedron newTet(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            newTet = Tetrahedron(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            Polyhedron sp;
                            sp.make_tetrahedron(newTet.vertex(0), newTet.vertex(1), newTet.vertex(2), newTet.vertex(3));
                            std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/tet_";
                            std::string tet_name = std::to_string(plane_count)+"_"+std::to_string(n)+std::to_string(p);
                            exportOFF(sp, path+tet_name);
                            tet = newTet;
                            tetIntersectionFun(tet, all_planes, vol, ++plane_count, tet_name);
                        }
                    }
                }
            }
        }
//        double newVol = abs(newTet.volume());
//        vol+=newVol;
//        tetIntersectionFun(newTet, all_planes, vol, ++plane_count);
    }


    if(ps == 0 && ns != 0){           // if positive side is empty, means no intersection with this plane
        if(plane_count != 3){           // if it was the last plane, (current) tet is contained in all planes
            tetIntersectionFun(tet, all_planes, vol, ++plane_count);
    }}
//    else if(ns == 0 && ps != 0)            // if negative side is empty, means no intersection at all
//        return 0;

//    if(ps == 0 && ns == 0){


    return 1;

}




void tetIntersectionTest(){



        Point sp0(1,1,1);
        Point sp1(0,4,1);
    //    Point sp1(1,3,1);
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
    //    Point dp0(0,2,5);
    //    Point dp1(0,4,5);
    //    Point dp2(2,3,5);
    //    Point dp3(2,2,0);
        Tetrahedron dtet(dp0,dp1,dp2,dp3);
        Polyhedron dp;
        dp.make_tetrahedron(dp0,dp1,dp2,dp3);
        exportOFF(dp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/dp");

        double cvol = 0.0;
        int pc = 0;
        tetIntersectionFun(dtet, planes, cvol, pc);

        CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> sensor_modifier(sp);
        Polyhedron_Exact sp_exact;
        sp_exact.delegate(sensor_modifier);
        Nef_polyhedron snef(sp_exact);

        CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> dp_modifier(dp);
        Polyhedron_Exact dp_exact;
        dp_exact.delegate(dp_modifier);
        Nef_polyhedron dnef(dp_exact);

        Nef_polyhedron fullnef = dnef*snef;

        Polyhedron_Exact intersection_tet_exact;
        fullnef.convert_to_polyhedron(intersection_tet_exact);

        CGAL::Polyhedron_copy_3<Polyhedron_Exact, Polyhedron::HalfedgeDS> modifier_rev(intersection_tet_exact);
        Polyhedron full_poly;
        full_poly.delegate(modifier_rev);
        exportOFF(full_poly, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/intersection");


        double vol_full1 = CGAL::Polygon_mesh_processing::volume(full_poly);

        std::cout << "my vol: " << cvol << "    nef vol: " << vol_full1 << std::endl;

        int a=5;


    //    int vlen = rand() % 8 + 4;
    //    std::vector<Point> points(vlen);
    //    for(int i = 0; i < vlen; i++){

    //        Point p(rand(), rand(), rand());
    //        points[i] = p;

    //    }
    //    Polyhedron Poly;²
    //    CGAL::convex_hull_3(points.begin(), points.end(),Poly);
    }
