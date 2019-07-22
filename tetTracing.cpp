#include <cgal_typedefs.h>
#include <fileIO.h>


namespace tetTracing{

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

            if(!do_intersect(a_triangles[i], b_triangles[j]))
                continue;

//            EPECK::Triangle_3 a_tri;
//            EPECK::Triangle_3 b_tri;
//            try{
//                a_tri = to_exact(a_triangles[i]);
//                b_tri = to_exact(b_triangles[j]);
//            }
//            catch(...){
//                continue;
//            }

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
//    std::cout << vol << std::endl;
    return vol;

//    exportOFF(a, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/a");
//    exportOFF(b, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/b");
//    exportOFF(c, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/c");

}

//double tetIntersectionFun(Polyhedron& a, Polyhedron& b){


//    Polyhedron::Facet_iterator sfi;
//    std::vector<Triangle> a_triangles;
//    for(sfi = a.facets_begin(); sfi != a.facets_end(); sfi++){
//        std::cout << "new face" << std::endl;
//        Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
//        std::vector<Point> triangle_points;
//        do{
//            std::cout << circ->vertex()->point() << std::endl;
//            triangle_points.push_back(circ->vertex()->point());
//        }
//        while (++circ != sfi->facet_begin());
//        Triangle tri(triangle_points[0], triangle_points[1], triangle_points[2]);
//        a_triangles.push_back(tri);
//    }

//    std::vector<Triangle> b_triangles;
//    for(sfi = b.facets_begin(); sfi != b.facets_end(); sfi++){
//        std::cout << "new face" << std::endl;
//        Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
//        std::vector<Point> triangle_points;
//        do{
//            std::cout << circ->vertex()->point() << std::endl;
//            triangle_points.push_back(circ->vertex()->point());
//        }
//        while (++circ != sfi->facet_begin());
//        Triangle tri(triangle_points[0], triangle_points[1], triangle_points[2]);
//        b_triangles.push_back(tri);
//    }

//    // double loop over the triangles of each polyhedron
//    std::set<Point> intersection_points;
//    for(int i=0; i<3; i++){
//        for(int j=0; j<3; j++){

//            CGAL::cpp11::result_of<Intersect(Triangle, Triangle)>::type result;

//            if(!a_triangles[i].is_degenerate() && !b_triangles[j].is_degenerate())
//                result = intersection(a_triangles[i], b_triangles[j]);

//            if (result){
//                if (const Segment* seg = boost::get<Segment>(&*result)){
//                    intersection_points.insert(seg->point(0));
//                    intersection_points.insert(seg->point(1));
//                }
//            }
//        }
//    }

//    double vol;
//    Polyhedron intersection_poly;
//    if(intersection_points.size() > 3){
//        CGAL::convex_hull_3(intersection_points.begin(), intersection_points.end(), intersection_poly);
//        if(intersection_poly.is_closed() && intersection_poly.is_valid())
//            vol = CGAL::Polygon_mesh_processing::volume(intersection_poly);
//    }
//    else{
//        vol = 0.0;
//    }
////    std::cout << vol << std::endl;
//    return vol;

////    exportOFF(a, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/a");
////    exportOFF(b, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/b");
////    exportOFF(c, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/c");

//}

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

        // iterate over all 3 points of the sensor triangle
        if(sensor_infos[k].size()<3)
            continue;
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
                    if(do_intersect(tri,ray)){

                            std::vector<Plane> planes(8);

                            // calc volume with halfspace intersections
                            // sensor tet
                            Point sp0 = sensor_infos[k][0]->point();
                            Point sp1 = sensor_infos[k][1]->point();
                            Point sp2 = sensor_infos[k][2]->point();
                            Point sp3 = sensor_infos[k][0]->info().sensor_pos;
                            Polyhedron sp;
                            sp.make_tetrahedron(sp0,sp1,sp2,sp3);

                            // Dt tet
                            Point tp0 = current_cell->vertex(0)->point();
                            Point tp1 = current_cell->vertex(1)->point();
                            Point tp2 = current_cell->vertex(2)->point();
                            Point tp3 = current_cell->vertex(3)->point();
                            Polyhedron tp;
                            tp.make_tetrahedron(tp0,tp1,tp2,tp3);

                            // intersection
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

//                            // calc volume with halfspace intersections
//                            // sensor tet
//                            Point sp0 = sensor_infos[k][0]->point();
//                            Point sp1 = sensor_infos[k][1]->point();
//                            Point sp2 = sensor_infos[k][2]->point();
//                            Point sp3 = sensor_infos[k][0]->info().sensor_pos;
//                            Polyhedron sp;
//                            sp.make_tetrahedron(sp0,sp1,sp2,sp3);
//                            Tetrahedron st(sp0,sp1,sp2,sp3);
//                            Point sp_centroid = CGAL::centroid(sp0,sp1,sp2,sp3);
//                            std::vector<Plane> planes(8);
//                            planes[0] = Plane(sp0,sp2,sp1);
//                            planes[1] = Plane(sp0,sp1,sp3);
//                            planes[2] = Plane(sp1,sp2,sp3);
//                            planes[3] = Plane(sp0,sp3,sp2);
//                            for(int i = 0; i<4; i++){
//                                if(!planes[i].has_on_negative_side(sp_centroid)){
//                                    planes[i]=planes[i].opposite();
//                                }
//                            }

//                            // Dt tet
//                            Point tp0 = current_cell->vertex(0)->point();
//                            Point tp1 = current_cell->vertex(1)->point();
//                            Point tp2 = current_cell->vertex(2)->point();
//                            Point tp3 = current_cell->vertex(3)->point();
//                            Polyhedron tp;
//                            tp.make_tetrahedron(tp0,tp1,tp2,tp3);
//                            Point tp_centroid = CGAL::centroid(tp0, tp1, tp2, tp3);
//                            planes[4] = Plane(tp0,tp2,tp1);
//                            planes[5] = Plane(tp0,tp1,tp3);
//                            planes[6] = Plane(tp1,tp2,tp3);
//                            planes[7] = Plane(tp0,tp3,tp2);
//                            for(int i = 4; i<8; i++){
//                                if(!planes[i].has_on_negative_side(tp_centroid)){
//                                    planes[i]=planes[i].opposite();
//                                }
//                            }
//                            // intersection
//                            Polyhedron P_full;
//                            double vol_full;
//                            try{
//                                CGAL::halfspace_intersection_with_constructions_3(std::begin(planes), std::end(planes),
//                                                                                  P_full);
//                                vol_full = CGAL::Polygon_mesh_processing::volume(P_full);
//                            }
//                            catch(...){
////                                std::cout<<"sp is valid: " << sp.is_valid() << "    sp is closed: " << sp.is_closed() << std::endl;
////                                std::cout<<"tp is valid: " << tp.is_valid() << "    tp is closed: " << tp.is_closed() << std::endl;
////                                std::cout << "is closed: "  << P_full.is_closed() << "  is valid: " << P_full.is_valid() << std::endl;
////                                exportOFF(P_full,"/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetras/fail_case");
////                                exportOFF(sp,"/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetras/sp");
////                                exportOFF(tp,"/home/raphael/Dropbox/Studium/PhD/data/sampleData/tetras/tp");
//                                vol_full = 0.0;
//                            }
//                            current_cell->info().outside_score+=vol_full;
//                            processed.insert(current_cell);

//                            // get neighbouring cells
//                            for(int ci = 0; ci < 4; ci++){
//                                Facet fac = std::make_pair(current_cell, ci);
//                                Facet mirror_fac = Dt.mirror_facet(fac);
//                                Cell_handle newCell = mirror_fac.first;
//                                if(processed.find(newCell) == processed.end()){
//                                    traverseCells(Dt, newCell, processed, planes, st, sp);
//                                }

//                            }

//                            // calc volume with nef
//                            // make nef from current sensor polygon
//                            Polyhedron sensor_poly;
//                            sensor_poly.make_tetrahedron(sensor_infos[k][0]->point(),       // which should also be equal to points[sensor_polys[k][0]]
//                                                     sensor_infos[k][1]->point(),
//                                                     sensor_infos[k][2]->point(),
//                                                    // for now just take the sensor position of the first point, but can also take a barycenter later
//                                                     sensor_infos[k][0]->info().sensor_pos);// which should also be equal to infos[sensor_polys[k][0].sensor_pos]
//                            Polyhedron_Exact sensor_poly_exact;
//                            CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> sensor_modifier(sensor_poly);
//                            sensor_poly_exact.delegate(sensor_modifier);
//                            bool closed = sensor_poly_exact.is_closed();
//                            bool valid = sensor_poly_exact.is_valid();
//    //                        exportOFF(sensor_poly_exact, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/poly_exact"+std::to_string(i));
//                            Nef_polyhedron sensor_nef(sensor_poly_exact);
//                            // make nef from current cell
//                            Polyhedron dt_poly;
//                            dt_poly.make_tetrahedron(current_cell->vertex(0)->point(),
//                                                     current_cell->vertex(1)->point(),
//                                                     current_cell->vertex(2)->point(),
//                                                     current_cell->vertex(3)->point());
//                            Polyhedron_Exact dt_poly_exact;
//                            CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> dt_modifier(dt_poly);
//                            dt_poly_exact.delegate(dt_modifier);
//                            Nef_polyhedron dt_nef(dt_poly_exact);
//                            // intersection
//                            Nef_polyhedron intersection_nef = dt_nef*sensor_nef;
////                            Nef_polyhedron::Volume_const_handle volit;
////                            volit->shells_begin();
////                            Polyhedron_Exact intersection_tet_exact;
////                            volit->shells_begin().convert_to_polyhedron(intersection_tet_exact);




//                            //                        typedef typename Nef_polyhedron::SFace_const_handle SFace_const_handle;
//                            //                        typedef typename Nef_polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;
//                            //                        Shell_entry_const_iterator seci;
//                            //                        intersection_nef.convert_inner_shell_to_polyhedron(SFace_const_handle(seci), intersection_tet_exact);

//                            Polyhedron_Exact intersection_tet_exact;
//                            if(intersection_nef.is_simple()){
//                                intersection_nef.convert_to_polyhedron(intersection_tet_exact);
//                                Polyhedron intersection_tet;
//                                CGAL::Polyhedron_copy_3<Polyhedron_Exact, Polyhedron::HalfedgeDS> tet_modifier(intersection_tet_exact);
//                                intersection_tet.delegate(tet_modifier);
//                                auto vol = CGAL::Polygon_mesh_processing::volume(intersection_tet);
////                                std::cout << double(vol) << std::endl;
//                                //   take volume of intersection nef and save it in the cell
//                                //   mark as traversed (for the current sensor_tet - thus needs to be unset for next sensor_tet iteration step)
//                                //   go to neighbours of current cell and check intersection there
//                                //   go to next cell
//                                //   go to next sensor polygon/tet
//                            }

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




//void iterateOverTetras(const Delaunay& Dt, std::vector<Point>& points, std::vector<vertex_info> infos, std::vector<std::vector<int>>& polys){

//    Polyhedron sensor_mesh;
//    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polys);
//    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polys, sensor_mesh);
////    std::vector<Polyhedron::Facet> degenerate_faces;
////    CGAL::Polygon_mesh_processing::remove_degenerate_faces(faces(sensor_mesh), sensor_mesh, std::back_inserter(degenerate_faces));
//    CGAL::Polygon_mesh_processing::remove_degenerate_faces(sensor_mesh);




//    // constructs AABB tree
//    AABB_Tree tree(faces(sensor_mesh).first, faces(sensor_mesh).second, sensor_mesh);
//    // maybe make the tree with a boost zip iterator to have an id
//    // for each facet that i can use to get the polygon in the polygon vector
//    // from the primitive_id

//    Delaunay::Finite_vertices_iterator vit;
//    for(vit = Dt.finite_vertices_begin(); vit != Dt.finite_vertices_end(); vit++){

////        int index = vit->info().idx;

//        Point current_point = vit->point();
//        std::vector<Primitive_id> primitives;
//        tree.all_intersected_primitives(current_point, std::back_inserter(primitives));

//        Polyhedron::Halfedge_around_facet_circulator fac;

//        // here needs to go an iterator at some point that iterates over all the sensor primitives of the primitive vector
//        auto prim_id = primitives[0]->facet_begin();
//        Point p1 = prim_id->vertex()->point();
//        int id1 = std::distance(points.begin(),std::find(points.begin(), points.end(), p1));
//        Point sensor_pos = infos[id1].sensor_pos;


////        unsigned long id1 = prim_id->vertex()->id();
//        prim_id++;
//        Point p2 = prim_id->vertex()->point();
//        prim_id++;
//        Point p3 = prim_id->vertex()->point();
//        // now form a nef tetra with p1,p2,p3 and sensor position
//        // and intersect with the current cell tetrahedra

////        prim_id++;
////        Point p4 = prim_id->vertex()->point();
//        // point p4 equals p1 again, which is correct.
//        // problem, segfault at some point and id = -1
//        // next step. get the correct ID, which is hopefully the same ID as the vertex_id
//        // of the vertex_info vector from which I can get the sensor point to form a sensor cell
//        // intersect the sensor cell with the current 3DT cell and all it's sourrounding 3DT
//        // cells and put the score on them
//        // that sould be ALL?!!

//        int a = 5;

////        for(fac = primitives[0]->facet_begin(); fac != primitives[0]->facet_end(); fac++){
////            Point p = fac->vertex()->point();
////            std::cout << p << " ";
////        }
////        std::endl;

//    }

//}







//void iterateOverTetras(const Delaunay& Dt, std::vector<Point>& points, std::vector<vertex_info> infos, std::vector<std::vector<int>>& polys){

//    // make Nef polyhedron from the Delaunay tetra
//    Delaunay::Finite_cells_iterator cit;
//    int j = 0;
//    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){

//        j++;
//        Point p0 = cit->vertex(0)->point();
//        Point p1 = cit->vertex(1)->point();
//        Point p2 = cit->vertex(2)->point();
//        Point p3 = cit->vertex(3)->point();

//        // The plane is oriented such that p, q and r are oriented in a positive sense (that is counterclockwise) when seen from the positive side of h.
//        // from: https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html
//        // and tetrahedron orientation can be found here: https://doc.cgal.org/latest/Triangulation_3/index.html
//        // and the cell centroid has to be on the NEGATIVE side of the plane
//        // DT planes
//        Plane planes[8];
//        planes[0] = Plane(p0,p2,p1);
//        planes[1] = Plane(p0,p1,p3);
//        planes[2] = Plane(p1,p2,p3);
//        planes[3] = Plane(p0,p3,p2);

//        // sensor planes
//        for(int i = 0; i < polys.size(); i++){

//            int i0 = polys[i][0]; int i1 = polys[i][1]; int i2 = polys[i][2];
//            Point p0 = points[i0];
//            Point p1 = points[i1];
//            Point p2 = points[i2];
//            Point p3 = infos[i0].sensor_pos;

//            planes[4] = Plane(p0,p2,p1);
//            planes[5] = Plane(p0,p1,p3);
//            planes[6] = Plane(p1,p2,p3);
//            planes[7] = Plane(p0,p3,p2);

//            Polyhedron P_full;
//            try{
//                CGAL::halfspace_intersection_with_constructions_3(std::begin(planes), std::end(planes), P_full);
//                double vol_full = CGAL::Polygon_mesh_processing::volume(P_full);
//                std::cout << "is closed: "  << P_full.is_closed() << "    full volume: " << vol_full << std::endl;
//            }
//            catch(...){}
//        }
//    }


////    for(int i = 0; i < sensor_planes.size(); i++){

////        Polyhedron P;
////        CGAL::halfspace_intersection_with_constructions_3(sensor_planes.begin(), sensor_planes.end(), P);



////    }
////    for(int i = 0; i < dt_planes.size(); i++){

////        Polyhedron P;
////        CGAL::halfspace_intersection_with_constructions_3(dt_planes[i].begin(), dt_planes[i].end(), P);



////    }





////    Polyhedron_Exact target;
//////    Nef_polyhedron target;
////    CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> modifier(P);
////    target.delegate(modifier);

////    Nef_polyhedron newNef(target);
//}

// COLMAP Delaunay meshing stems from
// P. Labatut, J‐P. Pons, and R. Keriven. "Robust and efficient surface
// reconstruction from range data". Computer graphics forum, 2009.
//
// and can be found in meshing.h in colmap/src/mvs

}

