//#include <cgal_typedefs.h>
//#include <fileIO.h>
//#include <rayTracing.cpp>
//#include <tetIntersection.cpp>

//namespace tetTracing{

///////////////////////////////////////////////////////////////////
///////////////////////// Tetrahedron tracing /////////////////////
///////////////////////////////////////////////////////////////////
/////
/////
/////
/////
///// cannot simply look for the sensor tetrahedron from the same 3 points, because the connection is broken when combining different sensors. Meaning I will have
///// Delaunay surface triangles that are formed by points from different sensors.
/////
/////
//void idxCells(Delaunay& Dt){
//    Delaunay::Finite_cells_iterator fci;
//    int i=0;
//    for(fci=Dt.finite_cells_begin();fci!=Dt.finite_cells_end();fci++){
//        fci->info().idx = i++;
//    }
//}
////void idxCells(Delaunay& Dt){
////    Delaunay::All_cells_iterator fci;
////    int i=0;
////    for(fci=Dt.all_cells_begin();fci!=Dt.all_cells_end();fci++){
////        if(!Dt.is_infinite(fci))
////            fci->info().idx = i++;
////        else
////            fci->info().outside_score = 0.0;
////    }
////}



//// TODO:
//// 1. do the full ray tracing to the outside, but not to the inside
//// simply replace the one_cell if-statement with one_cell && inside
//// 2. get the correct sensor orientation from COLMAP or from seperate depth map from each image from MicMac
//// 3. intersect a sensor topology tetrahedron (formed by 3 pixels next to each other, or LiDAR points next to each other and their (almost common -> barycenter) ray source
//// use this for outside vote of the Delaunay tetrahedra, and keep the ray for inside votes for now
//int traverseCells(Delaunay& Dt,
//                  Cell_handle& first_cell,
//                  Cell_handle& current_cell, std::unordered_set<Cell_handle>& processed,
//                  std::vector<Plane>& planes,
//                  Ray& ray){

//    // if current cell is not infinite then go on
//    // processed state is already checked before calling this function, so no need to check again
//    if(!Dt.is_infinite(current_cell))
//    {

////        int fidx = first_cell->info().idx;
////        int cidx = current_cell->info().idx;
////        std::cout << "first cell: " << fidx << "    second cell: " << cidx << std::endl;

////        for(int t = 0; t < 4; t++){
////            Triangle tri = Dt.triangle(current_cell, t);
////            bool dointersect = CGAL::do_intersect(ray,tri);
////            if(dointersect)

////        }
//        bool dointersect = false;
//        int t=0;
//        while(!dointersect && t < 4)
//        {
//            Triangle tri = Dt.triangle(current_cell, t++);
//            dointersect = CGAL::do_intersect(ray,tri);
//        }
//        if(!dointersect){
////            processed.insert(current_cell);
//            return 0;
//        }

//        Tetrahedron current_tet = Dt.tetrahedron(current_cell);

//        double vol = 0.0;
//        tetIntersectionFun(current_tet, planes, vol);
//        if(!isnan(vol)){
//            double score = vol/abs(current_tet.volume());
//            current_cell->info().outside_score+=score;
//        }
//        else{
//            std::cout << "NaN hit. Intersection? " << std::endl;
////            exportOFF(sp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/sp");
//            exportOFF(current_tet, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dp");
//            processed.insert(current_cell);
//            return 0;
//        }
//        processed.insert(current_cell);

//        // get neighbouring cells
//        for(int ci = 0; ci < 4; ci++){
//            Facet fac = std::make_pair(current_cell, ci);
//            Facet mirror_fac = Dt.mirror_facet(fac);
//            Cell_handle newCell = mirror_fac.first;
//            if(processed.find(newCell) == processed.end()){
//                traverseCells(Dt, first_cell, newCell, processed, planes, ray);
//            }
//        }
//    }
//    else{
//        // give infinite cell high value
//        current_cell->info().outside_score+=1;
//        processed.insert(current_cell);
//    }
//    return 1;
//}


//void firstCell(Delaunay& Dt, std::vector<std::vector<int>>& sensor_polys){
//    std::cout << "Start marking crossed tets..." << std::endl;
//    auto start = std::chrono::high_resolution_clock::now();

//    // save Delaunay vertex handle in sensor_infos vector for each sensor poly
//    std::vector<std::vector<Vertex_handle>> sensor_infos(sensor_polys.size());
//    Delaunay::Finite_vertices_iterator vit;
//    for(vit = Dt.finite_vertices_begin(); vit != Dt.finite_vertices_end(); vit++){
//        std::vector<int> sensor_tets = vit->info().sensor_tet;
//        for(int i = 0; i < sensor_tets.size(); i++){
//            int current_sensor_tet = sensor_tets[i];
//            sensor_infos[current_sensor_tet].push_back(vit);
//            int a = 5;
//        }
//    }
//    idxCells(Dt);


//    // iterate over all sensor triangles/tetrahedrons
//    int hit = 0;
//    for(int k = 0; k < sensor_polys.size(); k++){
//        std::unordered_set<Cell_handle> processed;
////        std::cout << "sensor poly " << k << std::endl;

//        if(sensor_infos[k].size()<3)
//            continue;

//        // iterate over all 3 points of the sensor triangle
//        for(int j = 0; j < 3; j++){

//            // get the corresponding Dt vertex handle for the current point
//            // so we are now considering the Dt vertex sensor_infos[k][j]
//            Vertex_handle current_vertex = sensor_infos[k][j];
//            Ray ray(current_vertex->point(), current_vertex->info().sensor_vec);
//            // vector of incident cells to the current vertex
//            std::vector<Cell_handle> inc_cells;
//            Dt.incident_cells(current_vertex, std::back_inserter(inc_cells));
//            for(std::size_t c=0; c < inc_cells.size(); c++){
//                Cell_handle current_cell = inc_cells[c];
////                int fcidx = current_cell->info().idx;
////                std::cout << "first cell: " << fcidx << std::endl;
//                // if current cell is not infinite and is not already processed (for this sensor poly) then go on
//                if(!Dt.is_infinite(current_cell)){
//                    if(processed.find(current_cell) == processed.end())
//                    {
//                        int cellBasedVertexIndex = current_cell->index(current_vertex);
//                        Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);

//                        bool dointersect = CGAL::do_intersect(ray,tri);
//                        if(dointersect){
//    //                        if(processed.find(current_cell) == processed.end()){
//                                // make a polyhedron of the sensor tet
//                                // not actually necessary anymore, just for outputting it, while debugging
//                                Point sp0 = sensor_infos[k][0]->point();
//                                Point sp1 = sensor_infos[k][1]->point();
//                                Point sp2 = sensor_infos[k][2]->point();
//                                Point sp3 = sensor_infos[k][j]->info().sensor_pos;
//                                Point spc = CGAL::centroid(sp0,sp1,sp2,sp3);
//                                Polyhedron sp;
//                                sp.make_tetrahedron(sp0,sp1,sp2,sp3);

//                                // The plane is oriented such that p, q and r are oriented in a positive sense (that is counterclockwise) when seen from the positive side of h.
//                                // from: https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html
//                                // and tetrahedron orientation can be found here: https://doc.cgal.org/latest/Triangulation_3/index.html
//                                // and the cell centroid has to be on the NEGATIVE side of the plane
//                                // sensor planes
//                                std::vector<Plane> planes(4);
//                                planes[0] = Plane(sp0,sp2,sp1);
//                                planes[1] = Plane(sp0,sp1,sp3);
//                                planes[2] = Plane(sp1,sp2,sp3);
//                                planes[3] = Plane(sp0,sp3,sp2);
//                                for(int i = 0; i<4; i++){
//                                    if(!planes[i].has_on_negative_side(spc)){
//                                        planes[i]=planes[i].opposite();
//                                    }

//                                }

//                                Tetrahedron current_tet = Dt.tetrahedron(current_cell);

//                                // intersection
//                                double vol = 0;
//                                tetIntersectionFun(current_tet, planes, vol);

//                                // TODO: investigate why there are sometimes NaNs!!
//                                if(!isnan(vol)){
//                                    double score = vol/abs(current_tet.volume());
//                                    current_cell->info().outside_score+=score;
//    //                                std::cout << "inside score: " << current_cell->info().inside_score <<
//    //                                             "  outside score " << current_cell->info().outside_score << std::endl;
//    //                                exportOFF(sp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/sp");
//    //                                exportOFF(current_tet, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dp");
//                                }
//                                else{
//                                    std::cout << "NaN hit. Intersection? " << std::endl;
//    //                                exportOFF(sp, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/sp");
//    //                                exportOFF(current_tet, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/TLS/failureCases/dp");
//                                    processed.insert(current_cell);
//                                    continue;
//                                }
//                                processed.insert(current_cell);

//                                // everything works quite well, but one problem is also that sensor mesh has wholes!!
//                                // thats why I have the errors in the door example
//                                // in fact, bad cells have inside and outside scores
//                                // so change to probability somehow?
//                                // furthermore also try with theory that all can be carved out with outside score
//                                // and triangle behind it is inside

//                                // update:
//                                // working well with carving out outside triangles. However, is it possible that
//                                // I am still missing outside triangles, due to some weird circumstances.
//                                // could try something like TODO:
//                                // if(close && vol==0)
//                                //      still check neighbouring cells

//                                //// traverse neighbouring cells
//                                for(int ci = 0; ci < 4; ci++){
//                                    Facet fac = std::make_pair(current_cell, ci);
//                                    Facet mirror_fac = Dt.mirror_facet(fac);
//                                    Cell_handle newCell = mirror_fac.first;
//                                    if(processed.find(newCell) == processed.end()){
//                                        traverseCells(Dt, current_cell, newCell, processed, planes, ray);
//                                    }
//                                }
//    //                        }//end of IF-already-processed
//                        }
//                    }
//                } // if cell is not infinite end
//                else{
//                    // give infinite cell high value
//                    current_cell->info().outside_score+=1;
//                    processed.insert(current_cell);
//                }
//            }
//        }//end of iteration over all three triangles of a sensor tetrahedron
//    }//end of iteration over all sensor tetrahedrons
//    auto stop = std::chrono::high_resolution_clock::now();
//    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    std::cout << "Tet intersection done in " << full_duration.count() << "s" << std::endl;
//}// end of firstCell function

//}// end of namespace


