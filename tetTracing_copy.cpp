#include <cgal_typedefs.h>
#include <fileIO.h>
#include <rayTracing.cpp>

namespace tetTracingCopy{

double tetIntersectionFun(Polyhedron& a, Polyhedron& b){

    if(!a.is_closed() && !b.is_closed() && !a.is_valid() && !b.is_valid())
        return 0.0;


//    typedef CGAL::Cartesian_converter<EPICK,EPECK>                         IK_to_EK;
//    typedef CGAL::Cartesian_converter<EPECK,EPICK>                         EK_to_IK;
    IK_to_EK to_exact;
    EK_to_IK to_inexact;


    Polyhedron::Facet_iterator sfi;
    std::vector<EPECK::Triangle_3> a_triangles;
    for(sfi = a.facets_begin(); sfi != a.facets_end(); sfi++){
        Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
        std::vector<EPECK::Point_3> triangle_points;
//        std::cout << "new face" << std::endl;
        do{
            triangle_points.push_back(to_exact(circ->vertex()->point()));
//            std::cout << circ->vertex()->point() << std::endl;
        }
        while (++circ != sfi->facet_begin());
        EPECK::Triangle_3 tri(triangle_points[0], triangle_points[1], triangle_points[2]);
        a_triangles.push_back(tri);
    }

    std::vector<EPECK::Triangle_3> b_triangles;
    for(sfi = b.facets_begin(); sfi != b.facets_end(); sfi++){
        Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
        std::vector<EPECK::Point_3> triangle_points;
//        std::cout << "new face" << std::endl;
        do{
            triangle_points.push_back(to_exact(circ->vertex()->point()));
//            std::cout << circ->vertex()->point() << std::endl;
        }
        while (++circ != sfi->facet_begin());
        EPECK::Triangle_3 tri(triangle_points[0], triangle_points[1], triangle_points[2]);
        b_triangles.push_back(tri);
    }

    // double loop over the triangles of each polyhedron
    std::set<Point> intersection_points;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){


//            if(!a_triangles[i].is_degenerate() && !b_triangles[j].is_degenerate())
//                continue;

            std::cout << a_triangles[i][0][0] << std::endl;
            std::cout << a_triangles[i][1][1] << std::endl;
            std::cout << b_triangles[i][0][0] << std::endl;
            std::cout << b_triangles[i][1][1] << std::endl;

            if(!do_intersect(a_triangles[i], b_triangles[j]))
                continue;


            EPECK::Triangle_3 a_tri = a_triangles[i];
            EPECK::Triangle_3 b_tri = b_triangles[j];


            CGAL::cpp11::result_of<EPECK::Intersect_3(EPECK::Triangle_3, EPECK::Triangle_3)>::type result;
            result = intersection(a_tri, b_tri);

            if (result){
                if (const EPECK::Segment_3* seg = boost::get<EPECK::Segment_3>(&*result)){
                    Point p1 = to_inexact(seg->point(0));
                    Point p2 = to_inexact(seg->point(1));
                    intersection_points.insert(p1);
                    intersection_points.insert(p2);
                }
            }
        }
    }

    double vol;
    Polyhedron intersection_poly;
    if(intersection_points.size() > 3){
        CGAL::convex_hull_3(intersection_points.begin(), intersection_points.end(), intersection_poly);
        if(intersection_poly.is_closed() && intersection_poly.is_valid())
            vol = CGAL::Polygon_mesh_processing::volume(intersection_poly);
    }
    else{
        vol = 0.0;
    }
    return vol;



}


/////////////////////////////////////////////////////////////////
/////////////////////// Tetrahedron tracing /////////////////////
/////////////////////////////////////////////////////////////////
///
///
///
///
/// cannot simply look for the sensor tetrahedron from the same 3 points, because the connection is broken when combining different sensors. Meaning I will have
/// Delaunay surface triangles that are formed by points from different sensors.
// TODO:
// 1. do the full ray tracing to the outside, but not to the inside
// simply replace the one_cell if-statement with one_cell && inside
// 2. get the correct sensor orientation from COLMAP or from seperate depth map from each image from MicMac
// 3. intersect a sensor topology tetrahedron (formed by 3 pixels next to each other, or LiDAR points next to each other and their (almost common -> barycenter) ray source
// use this for outside vote of the Delaunay tetrahedra, and keep the ray for inside votes for now
int traverseCells(const Delaunay& Dt,
                   Cell_handle& current_cell, std::unordered_set<Cell_handle>& processed,
                   std::vector<Plane>& planes, Polyhedron& sp){

    // if current cell is not infinite then go on
    // processed state is already checked before calling this function, so no need to check again
    if(!Dt.is_infinite(current_cell))
    {

//        for(int ci = 0; ci < 4; ci++){
//            Triangle tri = Dt.triangle(current_cell, 0);
//            if(do_intersect(tri,st)){
//                break;
//            }
//            return 0;
//        }
        // Dt tet
        Point tp0 = current_cell->vertex(0)->point();
        Point tp1 = current_cell->vertex(1)->point();
        Point tp2 = current_cell->vertex(2)->point();
        Point tp3 = current_cell->vertex(3)->point();
        Polyhedron tp;
        tp.make_tetrahedron(tp0,tp1,tp2,tp3);
//        Point tp_centroid = CGAL::centroid(tp0, tp1, tp2, tp3);
//        planes[4] = Plane(tp0,tp2,tp1);
//        planes[5] = Plane(tp0,tp1,tp3);
//        planes[6] = Plane(tp1,tp2,tp3);
//        planes[7] = Plane(tp0,tp3,tp2);
//        for(int i = 4; i<8; i++){
//            if(!planes[i].has_on_negative_side(tp_centroid)){
//                planes[i]=planes[i].opposite();
//            }
//        }
        // intersection
//        Polyhedron P_full;
//        double vol_full;
//        try{
//            CGAL::halfspace_intersection_with_constructions_3(std::begin(planes), std::end(planes),
//                                                              P_full);
//            vol_full = CGAL::Polygon_mesh_processing::volume(P_full);
//        }
//        catch(...){
////            std::cout<<"sp is valid: " << sp.is_valid() << "    sp is closed: " << sp.is_closed() << std::endl;
////            std::cout<<"tp is valid: " << tp.is_valid() << "    tp is closed: " << tp.is_closed() << std::endl;
//            std::cout << "is closed: "  << tp.is_closed() << "  is valid: " << tp.is_valid() << std::endl;
////            exportOFF(P_full,"/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetras/fail_case");
////            exportOFF(sp,"/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetras/sp");
////            exportOFF(tp,"/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetras/tp");
//            vol_full = 0.0;
//        }
        double vol_full = tetIntersectionFun(sp, tp);
        current_cell->info().outside_score+=vol_full;
        processed.insert(current_cell);

        // get neighbouring cells
        for(int ci = 0; ci < 4; ci++){
            Facet fac = std::make_pair(current_cell, ci);
            Facet mirror_fac = Dt.mirror_facet(fac);
            Cell_handle newCell = mirror_fac.first;
            if(processed.find(newCell) == processed.end()){
                traverseCells(Dt, newCell, processed, planes, sp);
            }
        }
    }
    else{
        // give infinite cell high value
        current_cell->info().outside_score+=10;
    }
    return 0;
}


void firstCell(const Delaunay& Dt, std::vector<Point>& points, std::vector<vertex_info>& infos, std::vector<std::vector<int>>& sensor_polys){

    auto start = std::chrono::high_resolution_clock::now();

    // save Delaunay vertex handle in sensor_infos vector for each sensor poly
    std::vector<std::vector<Vertex_handle>> sensor_infos(sensor_polys.size());
    Delaunay::Finite_vertices_iterator vit;
    int nv = Dt.number_of_vertices();
    for(vit = Dt.finite_vertices_begin(); vit != Dt.finite_vertices_end(); vit++){
        std::vector<int> sensor_tets = vit->info().sensor_tet;
        for(int i = 0; i < sensor_tets.size(); i++){
            int current_sensor_tet = sensor_tets[i];
            sensor_infos[current_sensor_tet].push_back(vit);
            int a = 5;
        }
    }

    // iterate over all sensor triangles/tetrahedrons
    for(int k = 0; k < sensor_polys.size(); k++){
        std::cout << k << std::endl;
        std::unordered_set<Cell_handle> processed;

        if(sensor_infos[k].size()<3)
            continue;

        // iterate over all 3 points of the sensor triangle
        for(int j = 0; j < 3; j++){

            // get the corresponding Dt vertex handle for the current point
            // so we are now considering the Dt vertex sensor_infos[k][j]
            Vertex_handle current_vertex = sensor_infos[k][j];
            Ray ray(current_vertex->point(), current_vertex->info().sensor_vec);
            // vector of incident cells to the current vertex (from vertex iterator vit)
            std::vector<Cell_handle> inc_cells;
            Dt.incident_cells(current_vertex, std::back_inserter(inc_cells));
            for(std::size_t i=0; i < inc_cells.size(); i++){
                Cell_handle current_cell = inc_cells[i];
                // if current cell is not infinite and is not already processed (for this sensor poly) then go on
                if(!Dt.is_infinite(current_cell) && processed.find(current_cell) == processed.end())
                {
                    int cellBasedVertexIndex = current_cell->index(current_vertex);
                    Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
//                    Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
                    // intersection from here: https://doc.cgal.org/latest/Kernel_23/group__intersection__linear__grp.html
                    // get the intersection of the ray and the triangle
//                    CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
//                      result = intersection(tri, ray);
                    // check if there is an intersection between the current ray and current triangle
//                    if(result){
//                        if(const Point* p = boost::get<Point>(&*result)){
                    Point intersectionPoint;
                    Point source = vit->point();
                    Vector rayV = ray.to_vector();
//                    bool result = rayTracing::rayTriangleIntersection(source, rayV, tri, intersectionPoint);

                    bool dointersect = CGAL::do_intersect(ray,tri);
//                    std::cout << result << " " << dointersect << std::endl;

                    if(dointersect){
//                        std::cout << "hit" << std::endl;

                        std::vector<Plane> planes(8);

//                            // calc volume with halfspace intersections
//                            // sensor tet
                        Point sp0 = sensor_infos[k][0]->point();
                        Point sp1 = sensor_infos[k][1]->point();
                        Point sp2 = sensor_infos[k][2]->point();
                        Point sp3 = sensor_infos[k][j]->info().sensor_pos;
                        Point spc = CGAL::centroid(sp0,sp1,sp2,sp3);
//                            Polyhedron sp;
//                            sp.make_tetrahedron(sp0,sp1,sp2,sp3);

//                            // Dt tet
                        Point tp0 = current_cell->vertex(0)->point();
                        Point tp1 = current_cell->vertex(1)->point();
                        Point tp2 = current_cell->vertex(2)->point();
                        Point tp3 = current_cell->vertex(3)->point();
                        Point tpc = CGAL::centroid(tp0,tp1,tp2,tp3);
//                            Polyhedron tp;
//                            tp.make_tetrahedron(tp0,tp1,tp2,tp3);

                        // calc volume with halfspace intersections

                        // The plane is oriented such that p, q and r are oriented in a positive sense (that is counterclockwise) when seen from the positive side of h.
                        // from: https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html
                        // and tetrahedron orientation can be found here: https://doc.cgal.org/latest/Triangulation_3/index.html
                        // and the cell centroid has to be on the NEGATIVE side of the plane
                        // DT planes
                        planes[0] = Plane(sp0,sp2,sp1);
                        planes[1] = Plane(sp0,sp1,sp3);
                        planes[2] = Plane(sp1,sp2,sp3);
                        planes[3] = Plane(sp0,sp3,sp2);
                        for(int i = 0; i<4; i++){
                            if(!planes[i].has_on_negative_side(spc)){
                                planes[i]=planes[i].opposite();
                            }
                        }

                        // sensor planes
                        planes[4] = Plane(tp0,tp2,tp1);
                        planes[5] = Plane(tp0,tp1,tp3);
                        planes[6] = Plane(tp1,tp2,tp3);
                        planes[7] = Plane(tp0,tp3,tp2);
                        for(int i = 4; i<8; i++){
                            if(!planes[i].has_on_negative_side(tpc)){
                                planes[i]=planes[i].opposite();
                            }
                        }


                        Polyhedron P_full;
                        CGAL::halfspace_intersection_3(std::begin(planes), std::end(planes), P_full);
                        bool close = P_full.is_closed();
                        Polyhedron convex_hull;
                        CGAL::convex_hull_3(P_full.points_begin(), P_full.points_end(), convex_hull);



//                            double vol_full = 0.0;
//                            if(P_full.is_closed())
                        double vol_full = CGAL::Polygon_mesh_processing::volume(convex_hull);



                        // intersection
//                            double vol_full = tetIntersectionFun(sp, tp);
                        current_cell->info().outside_score+=vol_full;
                        processed.insert(current_cell);

//                            // get neighbouring cells
//                            for(int ci = 0; ci < 4; ci++){
//                                Facet fac = std::make_pair(current_cell, ci);
//                                Facet mirror_fac = Dt.mirror_facet(fac);
//                                Cell_handle newCell = mirror_fac.first;
//                                if(processed.find(newCell) == processed.end()){
//                                    traverseCells(Dt, newCell, processed, planes, sp);
//                                }
//                            }

                        }

                    }

                }


        }

    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Tet intersection done in " << full_duration.count() << "s" << std::endl;

}



}

