#include <cgal_typedefs.h>
#include <fileIO.h>
#include <tetIntersection.h>



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


    bool exp = false;
//    bool exp = true;

    // path for intermediate tet export
    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/tet_";
//    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/tet";

//    std::cout << plane_count << std::endl;

    // when tet was recursivly called it goes to plane 4
    // which doesn't exist; hence algorithm is over and the current tets
    // volume should be counted
    if(plane_count == 4){
        double newVol = abs(tet.volume());
//        std::cout << "vol counted tet: " << tet_name << std::endl;
        vol+=newVol;
        return 1;
    }

    // get the current plane from the vector of planes
    Plane plane = all_planes[plane_count];
    // check if there are points of the Dt tet on the positive AND negative side of the
    // current plane
    std::vector<int> neg;
    std::vector<int> pos;
    std::vector<int> close;
    double epsilon = 0.00001;
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
        else{
            close.push_back(v);
        }
    }

    // continue according to intersection with the current plane
    int ns = neg.size();
    int ps = pos.size();
    int cs = close.size();
    // if there was an intersection left and right of the plane
    if(ns != 0 && ps != 0){
        // initialize the new tet
        Tetrahedron newTet;
        if(ps > ns)
        {
            std::vector<Tetrahedron> splitTets;
            // first make the 3 outside split tets, which are not really needed,
            // but its intersection points are needed, for afterwards making the single inside tet
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
                            // first add the intersection points of the previous tetrahedron
                            if(splitTets.size() > 0){
                                for(int st = 0; st < splitTets.size(); st++){
                                    newTetPoints.push_back(splitTets[st].vertex(3));
                                }
                            }
                            // add the positive point of the line that was split
                            newTetPoints.push_back(p2);
                            // add other positive points (will get less as the iteration moves on)
                            for(int pp = p+1; pp < ps; pp++)
                                newTetPoints.push_back(tet.vertex(pos[pp]));
                            // add close points
                            for(int c = 0; c < close.size(); c++)
                                newTetPoints.push_back(tet.vertex(close[c]));
                            newTetPoints.push_back(*intersection_point);
//                            Tetrahedron newTet(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            newTet = Tetrahedron(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            ////export
                            Polyhedron sp;
                            sp.make_tetrahedron(newTet.vertex(0), newTet.vertex(1), newTet.vertex(2), newTet.vertex(3));
                            std::string tet_name = std::to_string(plane_count)+"_"+std::to_string(n)+std::to_string(p);
//                            if(exp)
//                                exportOFF(sp, path+tet_name);
                            splitTets.push_back(newTet);
//                            tetIntersectionFun(newTet, all_planes, vol, plane_count+1, tet_name);
                        }
                    }
                }
            }
            // now form the remaining (inside) polyhedron and only call the tetIntersectionFun with that one
            std::vector<Point> newTetPoints;
            // add all the intersection points
            for(int st = 0; st < splitTets.size(); st++){
                newTetPoints.push_back(splitTets[st].vertex(3));
            }
            // need to add close points here as well,
            // because a shared edge of a split tet and Dt tet
            // leeds to only having 2 split tets, which in this positive sided
            // tetrahedron case here would not leed to a valid tetrahedron
            // because there are only two points coming from the loop above,
            // and not 3
            for(int c = 0; c < close.size(); c++)
                newTetPoints.push_back(tet.vertex(close[c]));
            // add the only negative point of this plane-cut
            newTetPoints.push_back(tet.vertex(neg[0]));
            newTet = Tetrahedron(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
            Polyhedron sp;
            sp.make_tetrahedron(newTet.vertex(0), newTet.vertex(1), newTet.vertex(2), newTet.vertex(3));
            std::string tet_name = std::to_string(plane_count)+"_pos";
            if(exp)
                exportOFF(sp, path+tet_name);
            // this is now a SINGLE tet that is completely inside the negative side of one plane
            // note: the difference here to the next else (ns > ps) is that there it is not possible
            // to form a single tet from the remaining cut lying to the inside of the one plane
            tetIntersectionFun(newTet, all_planes, vol, plane_count+1, tet_name);
        }
        else if(ns > ps)
        {
            std::vector<Tetrahedron> splitTets;
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
                            // and always replace the positive one with the new intersection point
                            std::vector<Point> newTetPoints;
                            // first add the intersection points of the previous tetrahedron
                            // in the first iteration there will be no points yet
                            if(splitTets.size() > 0){
                                for(int st = 0; st < splitTets.size(); st++){
                                    newTetPoints.push_back(splitTets[st].vertex(3));
                                }
                            }
                            // add the negative point of the line that was split as the first point of the split tet
                            newTetPoints.push_back(p1);
                            // add other negative points (will get less as the iteration moves on)
                            for(int nn = n+1; nn < ns; nn++)
                                newTetPoints.push_back(tet.vertex(neg[nn]));
                            // add close points
                            for(int c = 0; c < close.size(); c++)
                                newTetPoints.push_back(tet.vertex(close[c]));
                            newTetPoints.push_back(*intersection_point);
//                            Tetrahedron newTet(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            newTet = Tetrahedron(newTetPoints[0], newTetPoints[1], newTetPoints[2], newTetPoints[3]);
                            // ////export
                            Polyhedron sp;
                            sp.make_tetrahedron(newTet.vertex(0), newTet.vertex(1), newTet.vertex(2), newTet.vertex(3));
                            std::string tet_name = std::to_string(plane_count)+"_"+std::to_string(n)+std::to_string(p);
                            if(exp)
                                exportOFF(sp, path+tet_name);
                            splitTets.push_back(newTet);
                            tetIntersectionFun(newTet, all_planes, vol, plane_count+1, tet_name);
                        }
                    }
                }
            }
        }
        else{ // meaning ns == ps
            // make a vector where all 4 intersection points are saved
            std::vector<Point> intersectionPoints;
            for(int n = 0; n < ns; n++){
                // push back the two negative points
                for(int p = 0; p < ps; p++){
                    Point p1 = tet.vertex(neg[n]);
                    Point p2 = tet.vertex(pos[p]);
                    Line l(p1,p2);
                    CGAL::cpp11::result_of<Intersect(Line, Plane)>::type
                                  result = intersection(l, plane);
                    if(result){
                        if(const Point* intersection_point = boost::get<Point>(&*result)){
                            // close points
//                            for(int c = 0; c < close.size(); c++)
//                                newTetPoints.push_back(tet.vertex(close[c]));
                            intersectionPoints.push_back(*intersection_point);
                        }
                    }
                }
            }
            // now that all intersection points are calculated, for the 3 inside tets
            std::vector<Tetrahedron> splitTets;
            std::vector<std::string> tetNames;
            // first make the big one with two intersection points and two negative points
            newTet = Tetrahedron(intersectionPoints[0], intersectionPoints[1], tet.vertex(neg[0]), tet.vertex(neg[1]));
            // //export
            Polyhedron sp;
            sp.make_tetrahedron(newTet.vertex(0), newTet.vertex(1), newTet.vertex(2), newTet.vertex(3));
            tetNames.push_back(std::to_string(plane_count)+"_big");
            if(exp)
                exportOFF(sp, path+tetNames[0]);
            splitTets.push_back(newTet);
            // now the first small one
            newTet = Tetrahedron(intersectionPoints[0], intersectionPoints[1], intersectionPoints[2], tet.vertex(neg[1]));
            Polyhedron sp1;
            sp1.make_tetrahedron(newTet.vertex(0), newTet.vertex(1), newTet.vertex(2), newTet.vertex(3));
            tetNames.push_back(std::to_string(plane_count)+"_small1");
            if(exp)
                exportOFF(sp1, path+tetNames[1]);
            splitTets.push_back(newTet);
            // now the second small one
            newTet = Tetrahedron(intersectionPoints[1], intersectionPoints[2], intersectionPoints[3], tet.vertex(neg[1]));
            Polyhedron sp2;
            sp2.make_tetrahedron(newTet.vertex(0), newTet.vertex(1), newTet.vertex(2), newTet.vertex(3));
            tetNames.push_back(std::to_string(plane_count)+"_small2");
            if(exp)
                exportOFF(sp2, path+tetNames[2]);
            splitTets.push_back(newTet);
            // now continue with looping over planes for all three tets
            for(int t = 0; t < splitTets.size(); t++){
                tetIntersectionFun(splitTets[t], all_planes, vol, plane_count+1, tetNames[t]);
            }
        } // end of 2-2 case
    } // end of n != 0 && p != 0
    else if(ps == 0 && ns != 0){           // if positive side is empty, means no intersection with this plane
        tetIntersectionFun(tet, all_planes, vol, plane_count+1, tet_name);
    }
//    else if(ns == 0 && ps != 0)
//        return 0;
    // if ANY!! plane has negative side is empty (with the original tet!), means no intersection at all
    // of course it can be with cut tets; but then it will continue in another recursion loop
    // also works if one (ore more) points are close, and the rest empty
    else{
        return 0;
    }
    return 1;


}




void tetIntersectionTest(){

//    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/";
    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/tet_";


    system("exec rm -r /home/raphael/Dropbox/Studium/PhD/data/sampleData/tetTest/*");


//    Point sp01(1.3,1.3,1.3);
//    Point sp02(0.3,3.3,1.3);
//    Point sp03(2.4,2,1.3);


    Point sp0(1,1,1);
    Point sp1(0,4,1); // toy example 1
//    Point sp1(1,3,1);   // toy example 2
    Point sp2(3,2,1);
    Point sp3(2,2,4);

//    Polyhedron six;
//    six.make_tetrahedron(sp0,sp1,sp2,sp01);
//    six.make_tetrahedron(sp01,sp02,sp03,sp1);
//    six.make_tetrahedron(sp02,sp03,sp0,sp2);
//    exportOFF(six, path+"six");

    Point spc = CGAL::centroid(sp0,sp1,sp2,sp3);
    Tetrahedron stet(sp0,sp1,sp2,sp3);
    Polyhedron sp;
    sp.make_tetrahedron(sp0,sp1,sp2,sp3);
    exportOFF(sp, path+"sp");
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

//    Point dp0(3,1,5);
//    Point dp1(3,3,5);
//    Point dp2(1,2,5);
//    Point dp3(2,2,2);
    // example 3
    Point dp0(0,2,5);
    Point dp1(0,4,5);
    Point dp2(2,3,5);
    Point dp3(2,2,0);
    Tetrahedron dtet(dp0,dp1,dp2,dp3);
    Polyhedron dp;
    dp.make_tetrahedron(dp0,dp1,dp2,dp3);
    exportOFF(dp, path+"dp");

    double cvol = 0.0;
    int pc = 0;
    int ib = tetIntersectionFun(dtet, planes, cvol, pc);

    // NEF intersection
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
    exportOFF(full_poly, path+"intersection");

    double vol_full1 = CGAL::Polygon_mesh_processing::volume(full_poly);

    // compare
    std::cout << "intersection: " << ib << " my vol: " << cvol << "    nef vol: " << vol_full1 << std::endl;

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
