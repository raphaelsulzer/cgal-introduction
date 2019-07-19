#include <cgal_typedefs.h>
#include "fileIO.cpp"

// TODO: check if I can replace the whole thing with the nearest_vertex(const Point& p, Cell_handle start) function,
// which can be found in the Delaunay_triangulation_3.h file in /usr/lib/CGAL


////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
std::pair<float, float> cellScore(float dist2, double eig3, bool inside){

    float score_inside;
    float score_outside;
    // noise
    float sigma_d = eig3;
    // scene thickness // good for fontaine dataset is 0.1
    float sigma_o = 1.0;
    // scale of the outside area?? // good for fontaine dataset is 1.0
    float sigma_e = 1.0;
    // not to be confused with the following, which means if I am walking inside/outside
    if(inside){
        //        sigma = 0.05;
        score_inside = (1 - 0.5*exp(-pow((sqrt(dist2)/sigma_d),2)))*exp(-pow((sqrt(dist2)/sigma_o),2));
        score_outside = 0.5*exp(-pow((sqrt(dist2)/sigma_d),2));
    }
    else {
        score_outside = (1 - 0.5*exp(-pow((sqrt(dist2)/sigma_d),2)))*exp(-pow((sqrt(dist2)/sigma_e),2));
        score_inside = 0.5*exp(-pow((sqrt(dist2)/sigma_d),2));
    }
//    if(score_inside != score_inside || score_outside != score_outside)
//        std::cout << score_outside << "  " << score_inside << std::endl;

    // first element is the outside score, second the inside score
    std::pair<float,float> sigmas(score_outside, score_inside);
    return sigmas;
}

//float cellScore(float dist2, double eig3, bool inside){
      // implementation according to Efficient volumetric fusion paper

//    float sigma;
//    if(inside){
//        //        sigma = 0.05;
//        sigma = eig3;
//        // truncate inside ray
//        // in efficient volumetric fusion paper they also introduce outside limit of 3*sigma, simply for having "a shorter walk in the 3DT".
////        if(dist2 > 8*sigma)
////            dist2=0;
//    }
//    else {
////        sigma = 0.25;
//        // good with normals
////        sigma = eig3*5;
//        sigma = eig3;
//    }

//    float S = 1 - exp(-dist2/(2*sigma*sigma));
//    return S;
//}


int traverseCells(const Delaunay& Dt, double sigma, Ray ray, Cell_handle current_cell, int oppositeVertex, Point source, bool inside)
{

    // input::
    // &Delaunay            Reference to a Delaunay triangulation
    // ray                  Current ray
    // source               The "Delaunay" point of the current ray
    // current_cell         The cell that has just been entered (in the global context)
    // oppositeVertex       The opposite vertex of the facet where the current_cell was entered
    // inside               Bool that says if the current cell is before or after the point

    if(!Dt.is_infinite(current_cell)){
        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int idx = (oppositeVertex+i)%4;

            Triangle tri = Dt.triangle(current_cell, idx);
            double tri_area = tri.squared_area();
            if(tri_area < 1e-100){
                std::cout << "triangle with " << tri_area << " skipped." << std::endl;
//                return 0;
            }
            Facet fac = std::make_pair(current_cell, idx);

            // btw, here I don't have the problem of ray intersecting multiple cells, because I'm only checking in the current cell
            CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
              result = intersection(tri, ray);

            // check if there is an intersection between the current ray and current triangle
            if (result){
                // check if ray triangle intersection is a point (probably in most cases)
                // or a line segment (if ray lies inside the triangle)
                // if result is a point
                std::pair<float,float> score;
                if (const Point* p = boost::get<Point>(&*result))
                {//std::cout << "point of ray-triangle-intersection :  " << *p << std::endl;
                    // get the distance of this point to the current source:
//                    double dist = sqrt(CGAL::squared_distance(*p, source));
                    float dist2 = CGAL::squared_distance(*p, source);
                    //std::cout << dist2 << std::endl;
                    // dist2 = squared distance from intersection to the point; sigma = noise of the point; inside = bool
                    score = cellScore(dist2, sigma, inside);

                }
                else{
                    const Segment* s = boost::get<Segment>(&*result);
                    std::cout << "segment 3 intersection behind the first cell:  " << *s << std::endl;

                    // get the three edges of the current triangle
                    // check how they intersect with the current ray
                    // since the ray passes straight through the triangle in this case

    //                Point end_point = s->target();
    //                int vertexIndexOfEdgeCrossedRay = all_vertices.find(end_point)->second;

    //                std::cout << vertexIndexOfEdgeCrossedRay << std::cout;

                    // for now just return in this case, until it is solved
                    return 0;
                }
                // now locate the current cell in the global context of the triangulation,
                // so I can mark that it is crossed by a ray
                current_cell->info().outside_score += score.first;
                current_cell->info().inside_score += score.second;

                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                Facet mirror_fac = Dt.mirror_facet(fac);
                // now from this new cell that I am in (get it from mirror_fac), iterate over all the triangles that are not the mirror triangle
                // and check if there is an intersection
                // this should be entered again at if(!Dt.is_infinite(current_cell)), since like this I can check if the cell is not already the infinite cell
                // so start from there to put this into a function
                Cell_handle newCell = mirror_fac.first;
                int newIdx = mirror_fac.second;

                // for every vertex of the cell, check if there is an intersection, because it means that one cell has multiple triangles that are intersected by
                // the ray (since a vertex is shared by three facets).
                // That will generate an infinite loop, so it will always re-enter the same cell at some point
                // that's why I am - for now - just returning from that cell if it happens
                // what could be done is make this intersection point the new source of the ray
                // simply say ray source = this point, and ray target = new point in the opposite direction of the previous point
                // and than start from RayTracingFun again, because there it will just look for the opposite facet of the intersected vertex
                for(int i=0; i<4; i++){
                    Point pt = newCell->vertex(i)->point();
                    CGAL::cpp11::result_of<Intersect(Point, Ray)>::type
                      point_intersection = intersection(pt, ray);

                    if(point_intersection){
//                        const Point* p = boost::get<Point>(&*point_intersection);
//                        std::cout << "intersection with a vertex of the cell: " << *p << std::endl;
                        return 0;
                    }
                }
                traverseCells(Dt, sigma, ray, newCell, newIdx, source, inside);
            }
        }
    }
    // put outside score of infinite cell very high
    else{
        current_cell->info().outside_score+=1;
        current_cell->info().inside_score+=0;

    }
    return 0;
}

void firstCell(const Delaunay& Dt, Delaunay::Finite_vertices_iterator& vit, bool inside, bool one_cell){


    // get sigma of the current vertex
    double sigma = vit->info().sigma;


    // ray constructed from point origin to (end of) normal
    // introduces a ray r with source p and with a direction given by v.
    Ray ray(vit->point(), vit->info().sensor_vec);

    // make the inside ray
    // in fact MicMac saves the camera normals pointing away from the camera,
    // so I take the opposite ray for outside traversal and the normal ray for inside
    if(inside){
        ray = ray.opposite();
    }


    // vector of incident cells to the current vertex (from vertex iterator vit)
    std::vector<Cell_handle> inc_cells;
    Dt.incident_cells(vit, std::back_inserter(inc_cells));
    // for every cell of incident cells, check if facet(cell, vertex) intersects with vertex normal
    // so this is checking in all directions of a vertex, but we will only have an intersection in (maximum) one direction
    // why only in one direction? because I'm only checking the OPPOSITE facade. It of course also intersects with the bordering facets
    // of all the neighbouring cells
    for(std::size_t i=0; i < inc_cells.size(); i++){

        Cell_handle current_cell = inc_cells[i];


        if(!Dt.is_infinite(current_cell))
        {
            int cellBasedVertexIndex = current_cell->index(vit);
            Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
            Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);

            // intersection from here: https://doc.cgal.org/latest/Kernel_23/group__intersection__linear__grp.html
            // get the intersection of the ray and the triangle
            CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
              result = intersection(tri, ray);


            // check if there is an intersection between the current ray and current triangle
            if(result){
                // check if ray triangle intersection is a point (probably in most cases) or a line segment (if ray lies inside the triangle).
                // so far this is not used (and not needed) since I am not handling the unlikely case where the ray goes through a triangle
                // if result is a point
                std::pair<float,float> score;
                Point source = vit->point();
                if(const Point* p = boost::get<Point>(&*result)){

                //std::cout << "point of ray-triangle-intersection :  " << *p << std::endl;

                    // get the distance of this point to the current source:
//                    double dist = sqrt(CGAL::squared_distance(*p, source));
                    float dist2 = CGAL::squared_distance(*p, source);
                    // dist2 = squared distance from intersection to the point; sigma = noise of the point; inside = bool
                    score = cellScore(dist2, sigma, inside);
                }
                // else result is a line
                else{
                    const Segment* s = boost::get<Segment>(&*result);
                    std::cout << "segment 3 intersection in first cell:  " << *s << std::endl;
                    // for now just return in this case, until it is solved
                    //continue;
                    // TOOD: simply calculate the distance to this edge, and then I can also get a score from cellScore()
                    score = std::make_pair(0.0, 0.0);
                }
                // now locate the current cell in the global context of the triangulation,
                // so I can mark that it is crossed by a ray
                current_cell->info().outside_score += score.first;
                current_cell->info().inside_score += score.second;

                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                // now from this new cell that I am in (get it from mirror_fac), iterate over all the triangles that are not the mirror triangle
                // and check if there is an intersection
                // this should be entered again at if(!Dt.is_infinite(current_cell)), since like this I can check if the cell is not already the infinite cell
                // so start from there to put this into a function
                if(!inside){
                    Facet mirror_fac = Dt.mirror_facet(fac);
                    Cell_handle newCell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    // go to next cell
                    traverseCells(Dt, sigma, ray, newCell, newIdx, source, inside);
                }
            // if there was a match, break the loop around this vertex, so it can go to the next one
            // this is only done for speed reasons, it shouldn't have any influence on the label
                // because a ray can only hit more than one facet of a cell if it hits another point of the cell
                // in which case it goes THROUGH a facet of the cell, in which case the intersection is not a point
                // but an edge.
                // it does however make a difference if this is turn on or not. why??
            break;
            }
        }
        // put outside score of infinite cell very high
        else{
//            std::cout << "infinite cell score set" << std::endl;
            current_cell->info().outside_score+=1;
            current_cell->info().inside_score+=0;
        }
    }


}

void rayTracingFun(const Delaunay& Dt, bool one_cell){

    std::cout << "Start tracing rays to every point..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // from Efficient volumetric fusion paper, it follows that I need to sum over all rays to get the score for a tetrahedron
    // iterate over every vertices = iterate over every ray
    // TODO: go in here and map the following to a function
    // this function can than be called with the vft2 iterator
    // and I pass this iterator around and also pass it from within the traversal function
    // when I want to restart from a point that is hit
    Delaunay::Finite_vertices_iterator vit;
    for(vit = Dt.finite_vertices_begin() ; vit != Dt.finite_vertices_end() ; vit++){

        // collect outside votes
        firstCell(Dt, vit, 0, one_cell);    // one_cell currently not used in the correct way
        // collect inside votes
        firstCell(Dt, vit, 1, one_cell);    // one_cell currently not used in the correct way
    }
    // now that all rays have been traced, apply the last function to all the cells:
//    float gamma = 2.0;
//    Cell_map::iterator it;
//    for(it = all_cells.begin(); it!=all_cells.end(); it++)
//    {
//        std::get<1>(it->second) = 1 - exp(-std::get<1>(it->second)/gamma);
//        std::get<2>(it->second) = 1 - exp(-std::get<2>(it->second)/gamma);
//    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Ray tracing done in " << full_duration.count() << "s" << std::endl;
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


void iterateOverTetras(const Delaunay& Dt, std::vector<Point>& points, std::vector<vertex_info>& infos, std::vector<std::vector<int>>& sensor_polys){

//    Polyhedron sensor_mesh;
//    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, sensor_polys);
//    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, sensor_polys, sensor_mesh);
//    exportOFF(sensor_mesh, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/oriented_sensor_mesh");
////    std::vector<Polyhedron::Facet> degenerate_faces;
////    CGAL::Polygon_mesh_processing::remove_degenerate_faces(faces(sensor_mesh), sensor_mesh, std::back_inserter(degenerate_faces));
//    CGAL::Polygon_mesh_processing::remove_degenerate_faces(sensor_mesh);
//    exportOFF(sensor_mesh, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/cleaned_sensor_mesh");


//    std::vector<Point> new_points;
//    Polyhedron::Vertex_iterator svi;
//    int id = 0;
//    for(svi = sensor_mesh.vertices_begin(); svi != sensor_mesh.vertices_end(); svi++){
//        new_points.push_back(svi->point());
//        svi->id() = id++;
//    }
//    Polyhedron::Facet_iterator sfi;
//    id = 0;
//    for(sfi = sensor_mesh.facets_begin(); sfi != sensor_mesh.facets_end(); sfi++){
//        sfi->id() = id++;
//    }
//    Polyhedron::Edge_iterator sei;
//    id = 0;
//    for(sei = sensor_mesh.edges_begin(); sei != sensor_mesh.edges_end(); sei++){
//        sei->id() = id++;
//    }
//    std::vector<std::vector<int>> new_polys;
//    for(sfi = sensor_mesh.facets_begin(); sfi != sensor_mesh.facets_end(); sfi++){
//        Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
//        std::vector<int> ids;
//        do{ids.push_back(circ->vertex()->id());}
//        while ( ++circ != sfi->facet_begin());
//        new_polys.push_back(ids);
//    }


//    float stop_ratio=0.7;
//    SMS::LindstromTurk_params LT_params;
////    if(i_arg < argc) LT_params.BoundaryWeight = atof(argv[i_arg++]);
////    if(i_arg < argc) LT_params.VolumeWeight = atof(argv[i_arg++]);
////    if(i_arg < argc) LT_params.ShapeWeight = atof(argv[i_arg++]);
//    LT_params.BoundaryWeight = 0.3;
//    LT_params.VolumeWeight = 0.1;
//    LT_params.ShapeWeight = 0.1;

//    // This is a stop predicate (defines when the algorithm terminates).
//    SMS::Count_ratio_stop_predicate<Polyhedron> stop(stop_ratio);

//    // This the actual call to the simplification algorithm.
//    // The surface mesh and stop conditions are mandatory arguments.
//    // The index maps are needed because the vertices and edges
//    // of this surface mesh lack an "id()" field.
//    int r = SMS::edge_collapse
//              (sensor_mesh
//              ,stop
//              ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,sensor_mesh)) // for CGAL <4.6, remove ::parameters
//               .halfedge_index_map  (get(CGAL::halfedge_external_index,sensor_mesh))
//               .get_cost (SMS::Edge_length_cost<Polyhedron>())
//              );

//    std::cout << "\nFinished...\n" << r << " edges removed.\n"
//              << (sensor_mesh.size_of_halfedges()/2) << " final edges.\n";
//    exportOFF(sensor_mesh, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/decimated_sensor_mesh");


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
                if(!Dt.is_infinite(current_cell))
                {
                    int cellBasedVertexIndex = current_cell->index(current_vertex);
                    Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
//                    Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
                    // intersection from here: https://doc.cgal.org/latest/Kernel_23/group__intersection__linear__grp.html
                    // get the intersection of the ray and the triangle
                    CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
                      result = intersection(tri, ray);
                    // check if there is an intersection between the current ray and current triangle
                    if(result){

                        Point sp0 = sensor_infos[k][0]->point();
                        Point sp1 = sensor_infos[k][1]->point();
                        Point sp2 = sensor_infos[k][2]->point();
                        Point sp3 = sensor_infos[k][0]->info().sensor_pos;

                        Plane planes[8];
                        planes[0] = Plane(sp0,sp2,sp1);
                        planes[1] = Plane(sp0,sp1,sp3);
                        planes[2] = Plane(sp1,sp2,sp3);
                        planes[3] = Plane(sp0,sp3,sp2);

                        Point tp0 = current_cell->vertex(0)->point();
                        Point tp1 = current_cell->vertex(1)->point();
                        Point tp2 = current_cell->vertex(2)->point();
                        Point tp3 = current_cell->vertex(3)->point();

                        planes[4] = Plane(tp0,tp2,tp1);
                        planes[5] = Plane(tp0,tp1,tp3);
                        planes[6] = Plane(tp1,tp2,tp3);
                        planes[7] = Plane(tp0,tp3,tp2);

                        Polyhedron P_full;
                        try{
                            CGAL::halfspace_intersection_with_constructions_3(std::begin(planes), std::end(planes), P_full);
                            double vol_full = CGAL::Polygon_mesh_processing::volume(P_full);
                            std::cout << "is closed: "  << P_full.is_closed() << "    full volume: " << vol_full << std::endl;
                        }
                        catch(...){}



//                        // calc volume
//                        // make nef from current sensor polygon
//                        Polyhedron sensor_poly;
//                        sensor_poly.make_tetrahedron(sensor_infos[k][0]->point(),       // which should also be equal to points[sensor_polys[k][0]]
//                                                 sensor_infos[k][1]->point(),
//                                                 sensor_infos[k][2]->point(),
//                                                // for now just take the sensor position of the first point, but can also take a barycenter later
//                                                 sensor_infos[k][0]->info().sensor_pos);// which should also be equal to infos[sensor_polys[k][0].sensor_pos]
//                        Polyhedron_Exact sensor_poly_exact;
//                        CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> sensor_modifier(sensor_poly);
//                        sensor_poly_exact.delegate(sensor_modifier);
//                        bool closed = sensor_poly_exact.is_closed();
//                        bool valid = sensor_poly_exact.is_valid();
////                        exportOFF(sensor_poly_exact, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/poly_exact"+std::to_string(i));
//                        Nef_polyhedron sensor_nef(sensor_poly_exact);
//                        // make nef from current cell
//                        Polyhedron dt_poly;
//                        dt_poly.make_tetrahedron(current_cell->vertex(0)->point(),
//                                                 current_cell->vertex(1)->point(),
//                                                 current_cell->vertex(2)->point(),
//                                                 current_cell->vertex(3)->point());
//                        Polyhedron_Exact dt_poly_exact;
//                        CGAL::Polyhedron_copy_3<Polyhedron, Polyhedron_Exact::HalfedgeDS> dt_modifier(dt_poly);
//                        dt_poly_exact.delegate(dt_modifier);
//                        Nef_polyhedron dt_nef(dt_poly_exact);
//                        // intersection
//                        Nef_polyhedron intersection_nef = dt_nef*sensor_nef;
//                        Polyhedron_Exact intersection_tet_exact;
//                        intersection_nef.convert_to_polyhedron(intersection_tet_exact);
//                        Polyhedron intersection_tet;
//                        CGAL::Polyhedron_copy_3<Polyhedron_Exact, Polyhedron::HalfedgeDS> tet_modifier(intersection_tet_exact);
//                        intersection_tet.delegate(tet_modifier);
//                        auto vol = CGAL::Polygon_mesh_processing::volume(intersection_tet);
//                        std::cout << double(vol) << std::endl;
                        //   take volume of intersection nef and save it in the cell
                        //   mark as traversed (for the current sensor_tet - thus needs to be unset for next sensor_tet iteration step)
                        //   go to neighbours of current cell and check intersection there
                        //   go to next cell
                        //   go to next sensor polygon/tet
                        }

                    }

                }


        }

    }



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










