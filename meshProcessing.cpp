#include <cgal_typedefs.h>
#include <fileIO.h>

////////////////////////////////////////////////////////////
////////////////////// Triangulation ///////////////////////
////////////////////////////////////////////////////////////
Delaunay makeDelaunayWithInfo(std::vector<Point>& points, std::vector<vertex_info>& info)
{
    auto start = std::chrono::high_resolution_clock::now();

    // make the triangulation
    Delaunay Dt( boost::make_zip_iterator(boost::make_tuple(points.begin(), info.begin() )),
    boost::make_zip_iterator(boost::make_tuple(points.end(), info.end() ) )  );

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Triangulation done in " << duration.count() << "s" << std::endl;
    return Dt;
}


std::vector<Point> sampleTriangulation(std::vector<Point>& all_points)
{

    std::tuple<int, double> x_min;
    std::tuple<int, double> y_min;
    std::tuple<int, double> z_min;

    std::tuple<int, double> x_max;
    std::tuple<int, double> y_max;
    std::tuple<int, double> z_max;
    for(int i = 0; i < all_points.size(); i++)
    {
        double x_coord = all_points[i].x();
        double y_coord = all_points[i].y();
        double z_coord = all_points[i].z();



    }

}


////////////////////////////////////////////////////////////
/////////////////// preprocessing functions ////////////////
////////////////////////////////////////////////////////////

void preprocessSensorMesh(std::vector<Point>& points,
                          std::vector<vertex_info>& infos,
                          std::vector<std::vector<int>>& polys){


        auto start = std::chrono::high_resolution_clock::now();


        assert(points.size() == infos.size());
        std::map<Point, vertex_info> pim;
        for(int i = 0; i < points.size(); i++){
            pim[points[i]] = infos[i];
        }

        Polyhedron sensor_mesh;
        CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polys);
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polys, sensor_mesh);
//        exportOFF(sensor_mesh, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/oriented_sensor_mesh");
    //    std::vector<Polyhedron::Facet> degenerate_faces;
    //    CGAL::Polygon_mesh_processing::remove_degenerate_faces(faces(sensor_mesh), sensor_mesh, std::back_inserter(degenerate_faces));
        CGAL::Polygon_mesh_processing::remove_degenerate_faces(sensor_mesh);
//        exportOFF(sensor_mesh, "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/cleaned_sensor_mesh");


        std::vector<Point> new_points;
        std::vector<vertex_info> new_infos;
        Polyhedron::Vertex_iterator svi;
        int id = 0;
        for(svi = sensor_mesh.vertices_begin(); svi != sensor_mesh.vertices_end(); svi++){
            new_points.push_back(svi->point());
            vertex_info info = pim.find(svi->point())->second;
            new_infos.push_back(info);
            svi->id() = id++;
        }
        Polyhedron::Facet_iterator sfi;
        id = 0;
        for(sfi = sensor_mesh.facets_begin(); sfi != sensor_mesh.facets_end(); sfi++){
            sfi->id() = id++;
        }
//        Polyhedron::Edge_iterator sei;
//        id = 0;
//        for(sei = sensor_mesh.edges_begin(); sei != sensor_mesh.edges_end(); sei++){
//            sei->id() = id++;
//        }
        std::vector<std::vector<int>> new_polys;
        for(sfi = sensor_mesh.facets_begin(); sfi != sensor_mesh.facets_end(); sfi++){
            Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
            std::vector<int> ids;
            do{ids.push_back(circ->vertex()->id());}
            while ( ++circ != sfi->facet_begin());
            new_polys.push_back(ids);
        }

        points = new_points;
        infos = new_infos;
        polys = new_polys;


        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "Mesh preprocessing done in " << duration.count() << "s" << std::endl;

}


//void preprocessSensorMesh(){

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

//}











void decimateSurfaceMesh(std::string input, std::string output)
{
    // optional
    float stop_ratio=0.5;


    SMS::LindstromTurk_params LT_params;
//    if(i_arg < argc) LT_params.BoundaryWeight = atof(argv[i_arg++]);
//    if(i_arg < argc) LT_params.VolumeWeight = atof(argv[i_arg++]);
//    if(i_arg < argc) LT_params.ShapeWeight = atof(argv[i_arg++]);

    clock_t start = clock();
    Polyhedron surface_mesh;

    std::ifstream is(input.c_str());
    std::ofstream os(output.c_str());
    std::string line;
    getline(is, line);
    if(line != "OFF") std::cout << "ERROR: Not an OFF file, starts with " << line << std::endl;
    os << "OFF" << std::endl;
    getline(is, line);
    while(line[0] == '#') {os << line << std::endl; getline(is, line);}
    is.seekg(0);
    is >> surface_mesh;
    // This is a stop predicate (defines when the algorithm terminates).
    SMS::Count_ratio_stop_predicate<Polyhedron> stop(stop_ratio);

    // This the actual call to the simplification algorithm.
    // The surface mesh and stop conditions are mandatory arguments.
    // The index maps are needed because the vertices and edges
    // of this surface mesh lack an "id()" field.
    int r = SMS::edge_collapse
              (surface_mesh
              ,stop
              ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface_mesh)) // for CGAL <4.6, remove ::parameters
               .halfedge_index_map  (get(CGAL::halfedge_external_index  ,surface_mesh))
               .get_cost (SMS::LindstromTurk_cost<Polyhedron>(LT_params))
               .get_placement(SMS::LindstromTurk_placement<Polyhedron>(LT_params))
              );

    std::cout << "\nFinished...\n" << r << " edges removed.\n"
              << (surface_mesh.size_of_halfedges()/2) << " final edges.\n" ;

    os << '#'; // we copied the header to keep the comments
    os << surface_mesh ;
    std::cout << "Done in " << (double)(clock() - start) / (double) CLOCKS_PER_SEC << "s" << std::endl;


}





////////////////////////////////////////////////////////////
/////////////////// postprocessing functions ////////////////
////////////////////////////////////////////////////////////
void createSurfaceMesh(const Delaunay& Dt, std::vector<Point>& points, std::vector<std::vector<int>>& polygons,
                       Polyhedron& out_mesh,
                       int orient, int nb_components_to_keep)
{

    // POINTS:

    // TODO:
    // in fact here I am still using all the original points from the point set, and not only the points that are left after the
    // GC optimization. this should be changed
    int index = 0;
    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        vft->info().idx = index;
        index++;
        points.push_back(vft->point());
    }


    // FACETS:

    // TODO:
    // could also add the orientation from the exportSurfacePLY function
    // maybe it helps to speed up the surface creation below.

    // initialise cell and vertex handle
    int vidx;
    Cell_handle c;
    Vertex_handle v;
    Delaunay::Finite_facets_iterator fft;
    for(fft = Dt.finite_facets_begin() ; fft != Dt.finite_facets_end() ; fft++){

        // get vertex and cell index that describes theurface_mesh facet
        // facet fft is represented by std::pair(cell c, int vidx). vidx is the vertex opposite to the cell.
        // even though some of the facets may be described by infinite cells, the facet is still has a neighbouring cell that is finite.
        // see: https://doc.cgal.org/latest/Triangulation_3/index.html
        c = fft->first;         // cell
        vidx = fft->second;     // vertex index

        //////// check which faces to prune:

//        // with GC labelling:
        int clabel = c->info().final_label;
        Cell_handle m = Dt.mirror_facet(*fft).first;
        int mlabel = m->info().final_label;

        // check if two neighbouring cells have the same label, and if so (and the prunce_faces export function is active) continue to next face
        if(clabel == mlabel){
            continue;
        }

        std::vector<int> polygon_indecis;
        // if label of neighbouring cells is not the same...
        // start printed facet line with a 3
        // if opposite vertex vidx is 2, we start at j = vidx + 1 = 3, 3%4 = 3
        // next iteration: j = 4, 4%4 = 0, next iteration: j = 5, 5%4 = 1;
        // so we exactely skip 2 - the opposite vertex.
        for(int j = vidx + 1 ; j <= vidx + 3 ; j++){
            // print the indicies of each cell to the file
            // vertices is a map of all vertices of the triangulation to an index
            v = c->vertex(j%4);
            // add the stuff to a list of Polygons
//            polygon_indecis.push_back(Vertices.find(v)->second);
            polygon_indecis.push_back(v->info().idx);
        }
        polygons.push_back(polygon_indecis);

    }

    int oriented;
    if(orient > 0){oriented = CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);}
    else{oriented = -1;}

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points,polygons,out_mesh);

    std::cout << "Surface mesh created..." << std::endl;
    std::cout << "Surface mesh oriented: " << oriented << std::endl;
    std::cout << "... 0 meaning additional vertices were added for orientation" << std::endl;

//    Kernel kernel;
//    double max_dist = CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(points,polygons,0.35,kernel);


    // export a Polyhedron surface mesh (as .OFF)
    int erased_components;
    if(nb_components_to_keep > 0){erased_components = out_mesh.keep_largest_connected_components(nb_components_to_keep);}
    else{erased_components = -1;}

    std::cout << "Number of erased components from surface mesh: " << erased_components << std::endl;


}
