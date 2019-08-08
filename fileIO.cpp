#include <cgal_typedefs.h>
#include <fileIO.h>

#include "plyDefinition.cpp"
#include "colmapPLY.cpp"


//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
//////////////////////////////////////////////////////////
// generate a Delaunay triangulation from a PLY file
Delaunay triangulationFromFile(std::string ifn, std::vector<PNC> ply_lines)
{
    auto start = std::chrono::high_resolution_clock::now();


    //    std::vector<PC> points; // store points
    std::ifstream in(ifn);

    // read the ply
    typedef CGAL::Nth_of_tuple_property_map<0, PNC> Point_map;
    typedef CGAL::Nth_of_tuple_property_map<1, PNC> Normal_map;
    typedef CGAL::Nth_of_tuple_property_map<2, PNC> Color_map;
    CGAL::read_ply_points_with_properties
      (in,
       std::back_inserter (ply_lines),
       CGAL::make_ply_point_reader (Point_map()),
       CGAL::make_ply_normal_reader (Normal_map()),
       std::make_tuple (Color_map(),
                               CGAL::Construct_array(),
                               CGAL::PLY_property<unsigned char>("red"),
                               CGAL::PLY_property<unsigned char>("green"),
                               CGAL::PLY_property<unsigned char>("blue"))
       );

    std::vector<Point> points;
//    std::vector<IdxSigNormCol> infos;
    std::vector<vertex_info> infos;
    for (int i = 0; i < ply_lines.size (); ++ i)
    {
        // make vector of points
        points.push_back(get<0>(ply_lines[i]));
        // make vector of infos as: tuple(idx, sigma, normal, color)
        vertex_info inf;
        inf.idx = i;
        inf.color = get<2>(ply_lines[i]);
        inf.normal = get<1>(ply_lines[i]);
        inf.sigma = 0.0;
        infos.push_back(inf);
    }

    // make the triangulation
    Delaunay Dt( boost::make_zip_iterator(boost::make_tuple(points.begin(), infos.begin() )),
    boost::make_zip_iterator(boost::make_tuple(points.end(), infos.end() ) )  );

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Triangulation done in " << duration.count() << "s" << std::endl;
    return Dt;
}

Delaunay triangulationFromFile(std::string ifn, std::vector<PN> ply_lines)
{
    auto start = std::chrono::high_resolution_clock::now();

    //    std::vector<PC> points; // store points
    std::ifstream in(ifn);

    // read the ply
    typedef CGAL::Nth_of_tuple_property_map<0, PN> Point_map;
    typedef CGAL::Nth_of_tuple_property_map<1, PN> Normal_map;
    CGAL::read_ply_points_with_properties
      (in,
       std::back_inserter (ply_lines),
       CGAL::make_ply_point_reader (Point_map()),
       CGAL::make_ply_normal_reader (Normal_map())
       );

    std::vector<Point> points;
//    std::vector<IdxSigNormCol> infos;
    std::vector<vertex_info> infos;
    for (int i = 0; i < ply_lines.size (); ++ i)
    {
        // make vector of points
        points.push_back(get<0>(ply_lines[i]));
        // make vector of infos as: tuple(idx, sigma, normal, color)
        vertex_info inf;
        inf.idx = i;
        unsigned char red = 255; unsigned char green = 255; unsigned char blue = 255;
        Color col = CGAL::make_array(red, green, blue);
        inf.color = col;
        inf.normal = get<1>(ply_lines[i]);
        inf.sigma = 0.0;
        infos.push_back(inf);
    }

    // make the triangulation
    Delaunay Dt( boost::make_zip_iterator(boost::make_tuple(points.begin(), infos.begin() )),
    boost::make_zip_iterator(boost::make_tuple(points.end(), infos.end() ) )  );

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Triangulation done in " << duration.count() << "s" << std::endl;
    return Dt;
}


void readColmapPLY(){


    auto start = std::chrono::high_resolution_clock::now();

    std::string path = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/reconstruction/dense_matching/fused";

    std::string ifn=path+".ply";

    ////// get the sensor infos
    // this reads the visibility info file, which has indicies for each point that correspond to the "index" of the camera
    std::string sensor_file = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/reconstruction/dense_matching/fused.ply.vis";
    std::vector<std::vector<int>> vis_info;
    loadVisibilityFile(sensor_file.c_str(), vis_info);
    // this reads the actual position of the camera corresponding to its index
    std::map<int, Point> sensor_map;
    std::string image_file = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/reconstruction/sfm_result/images.bin";
    loadImageFile(image_file.c_str(), sensor_map);

    ///// read Binary PLY, and connect to the sensor infos
    Mesh_ply aMesh;
    Import_PLY(ifn.c_str(), &aMesh);
    std::vector<Point> points;
    std::vector<vertex_info> infos;
    for(int i = 0; i < aMesh.mVertices.size(); i++){
        // save point
        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
        points.push_back(pt);
        // sensor
        vertex_info vec_inf;
        // get the first of the cameras of this point
        // TODO: in the future I could make sensor pos a vector of points, and not only save one camera per point, but all of them
        Point sensor_location = sensor_map.find(vis_info[i][0])->second;
        // sensor posiion/point
        vec_inf.sensor_pos = sensor_location;
        // sensor vector
        Vector vec = sensor_location - pt;
        vec_inf.sensor_vec = vec;
        // color
        unsigned char r,g,b;
        if(aMesh.mvColors.size() > 0){
            r = aMesh.mvColors[i].x * 256;
            g = aMesh.mvColors[i].y * 256;
            b = aMesh.mvColors[i].z * 256;
        }
        std::array<unsigned char, 3> col = {r,g,b};
        vec_inf.color = col;
        // index
        vec_inf.idx = i;
        // save vertex info
        infos.push_back(vec_inf);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "File " << ifn << " read in " << duration.count() << "s" << std::endl;

    // export the PLY with points and sensor position
    exportPoints(path, points, infos);

}


void readAP(std::string ifn,
                   std::vector<Point>& points, std::vector<vertex_info>& infos)
{
    auto start = std::chrono::high_resolution_clock::now();

    // read Binary PLY with sensor
    Mesh_ply aMesh;
    Import_PLY(ifn.c_str(), &aMesh);

    for(int i = 0; i < aMesh.mVertices.size(); i++){
        // save points
        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
        points.push_back(pt);
        // sensor
        vertex_info vec_info;
        Point sensor_position(aMesh.mvCapture[i].x, aMesh.mvCapture[i].y, aMesh.mvCapture[i].z);
        // apply Transformation from cloud compare registration to TLS to the sensor location
        Eigen::Vector4d slh(sensor_position.x(), sensor_position.y(), sensor_position.z(), 1.0);
        Eigen::Matrix4d RT;
        // this rotation+translation matrix is coming from Cloud Compare where I registered the
        // AP to the TLS point cloud.
        // now this also needs to be applied to the sensor position.
        // last column has the global translation shift added, that Cloud Compare applies when opening the file,
        // but of course not adds to the translation matrix
        RT << -114.147835, -6.538993, 0.458483, -13.759159+51200,
              -6.528348, 114.125534, 2.332502, 117.205711+91900,
              -0.591047, 2.302480, -114.311195, 849.723450,
              0.000000, 0.000000, 0.000000, 1.000000;
        Eigen::Vector4d RTslh = RT*slh;
        vec_info.sensor_pos = Point(RTslh[0] / RTslh[3],
                                    RTslh[1] / RTslh[3],
                                    RTslh[2] / RTslh[3]);
        Vector vec(vec_info.sensor_pos.x() - pt.x(), vec_info.sensor_pos.y() - pt.y(), vec_info.sensor_pos.z() - pt.z());
        vec_info.sensor_vec = vec;
        // color
        unsigned char r,g,b;
        if(aMesh.mvColors.size() > 0){
            r = aMesh.mvColors[i].x * 256;
            g = aMesh.mvColors[i].y * 256;
            b = aMesh.mvColors[i].z * 256;
        }
        std::array<unsigned char, 3> col = {r,g,b};
        vec_info.color = col;
        // index, e.g. used in the concatenate data function
        vec_info.idx = i;
        // save vertex info
        infos.push_back(vec_info);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "File " << ifn << " read in " << duration.count() << "s" << std::endl;
}





void readTLS(std::string ifn,
                   std::vector<Point>& points, std::vector<vertex_info>& infos,
                   std::vector<std::vector<int>>& sensor_triangle){

    auto start = std::chrono::high_resolution_clock::now();

    // read Binary PLY with sensor
    Mesh_ply aMesh;
    Import_PLY(ifn.c_str(), &aMesh);

    for(int i = 0; i < aMesh.mVertices.size(); i++){
        // save points
        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
        points.push_back(pt);
        // sensor
        Vector vec(aMesh.mvCapture[i].x - pt.x(), aMesh.mvCapture[i].y - pt.y(), aMesh.mvCapture[i].z - pt.z());
        vertex_info vec_inf;
        vec_inf.sensor_vec = vec;
        vec_inf.sensor_pos = Point(aMesh.mvCapture[i].x, aMesh.mvCapture[i].y, aMesh.mvCapture[i].z);
        // color
        unsigned char r,g,b;
        if(aMesh.mvColors.size() > 0){
            r = aMesh.mvColors[i].x * 256;
            g = aMesh.mvColors[i].y * 256;
            b = aMesh.mvColors[i].z * 256;
        }
        std::array<unsigned char, 3> col = {r,g,b};
        vec_inf.color = col;
        // save vertex info
        infos.push_back(vec_inf);
    }

    // save incident sensor triangles in each 3DT vertex
    for(int i = 0; i < aMesh.mIndices.size()/3; i++){
        std::vector<int> poly(3);
        int id0 = aMesh.mIndices[(i*3)+0];
        int id1 = aMesh.mIndices[(i*3)+1];
        int id2 = aMesh.mIndices[(i*3)+2];
        poly[0] = id0;
        poly[1] = id1;
        poly[2] = id2;
        sensor_triangle.push_back(poly);

        // save incident sensor triangles in each 3DT vertex
        infos[id0].sensor_tet.push_back(i);
        infos[id1].sensor_tet.push_back(i);
        infos[id2].sensor_tet.push_back(i);
        // basically just integrating the loop below in here.

    }
//    // save incident sensor triangles in each 3DT vertex
//    for(int i = 0; i < sensor_triangle.size(); i++){
//        for(int j = 0; j < 3; j++){
//            int current_vertex = sensor_triangle[i][j];
//            infos[current_vertex].sensor_tet.push_back(i);
//        }
//    }


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "File " << ifn << " read in " << duration.count() << "s" << std::endl;
}


void readASCIIPLY(std::string ifn, std::vector<Point>& points, std::vector<vertex_info>& infos)
{

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<PNC> ply_lines;

    //    std::vector<PC> points; // store points
    std::ifstream in(ifn);

    // read the ply
    typedef CGAL::Nth_of_tuple_property_map<0, PNC> Point_map;
    typedef CGAL::Nth_of_tuple_property_map<1, PNC> Normal_map;
    typedef CGAL::Nth_of_tuple_property_map<2, PNC> Color_map;
    CGAL::read_ply_points_with_properties
      (in,
       std::back_inserter (ply_lines),
       CGAL::make_ply_point_reader (Point_map()),
       CGAL::make_ply_normal_reader (Normal_map()),
       std::make_tuple (Color_map(),
                               CGAL::Construct_array(),
                               CGAL::PLY_property<unsigned char>("red"),
                               CGAL::PLY_property<unsigned char>("green"),
                               CGAL::PLY_property<unsigned char>("blue"))
       );

    for (int i = 0; i < ply_lines.size (); ++ i)
    {
        // make vector of points
        points.push_back(get<0>(ply_lines[i]));
        // make vector of infos as: tuple(idx, sigma, normal, color)
        vertex_info inf;
        inf.idx = i;
        inf.color = get<2>(ply_lines[i]);
        inf.sensor_vec = get<1>(ply_lines[i]);
        inf.sigma = 0.0;
        infos.push_back(inf);
    }


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "File " << ifn << " read in " << duration.count() << "s" << std::endl;
}


void concatenateData(std::vector<Point>& a_points, std::vector<vertex_info>& a_info,
                     std::vector<Point>& t_points, std::vector<vertex_info>& t_info,
                     bool copyInfo = 1)
{

    // two versions of this function, one with copy info, and one where I give a specific color per dataset
    if(copyInfo){
        Incremental_Tree tree(boost::make_zip_iterator(boost::make_tuple( a_points.begin(),a_info.begin() )),
                  boost::make_zip_iterator(boost::make_tuple( a_points.end(), a_info.end() ) ));

        for(int i = 0; i < t_points.size(); i++)
        {
            Incremental_neighbor_search search(tree, t_points[i]);
            Incremental_neighbor_search::iterator it = search.begin();

            // index of nearest neighbor in AP
            int idx = boost::get<1>(it->first).idx;

    //        t_info[i].normal = a_info[idx].normal;
            t_info[i].color = a_info[idx].color;

            a_points.push_back(t_points[i]);
            a_info.push_back(t_info[i]);
        }

    }
    else{
        std::array<unsigned char, 3> red = {240,128,128};
        std::array<unsigned char, 3> blue = {30,144,255};

        for(int i = 0; i < a_info.size(); i++)
        {
            a_info[i].color = red;
        }
        for(int i = 0; i < t_info.size(); i++)
        {
            t_info[i].color = blue;
            a_points.push_back(t_points[i]);
            a_info.push_back(t_info[i]);
        }
    }
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////// OUTPUT //////////////////////////////
/////////////////////////////////////////////////////////////////////
void exportEdges(std::fstream& fo, const Delaunay& Dt, const Cell_map& all_cells, const Vertex_map& all_vertices)
{

    std::cout << "Export edges start" << std::endl;

//    std::fstream fo;
//    fo.open(ofn, std::fstream::out);

    int edge_count = 0;

    Delaunay::Finite_edges_iterator edge;
    for(edge = Dt.finite_edges_begin(); edge != Dt.finite_edges_end(); edge++){

        Cell_handle current_cell = edge->first;

        int current_label = std::get<3>(all_cells.find(current_cell)->second);
        int current_index = std::get<0>(all_cells.find(current_cell)->second);

        for(int i = 0; i < 4; i++){

            Cell_handle neighbour_cell = current_cell->neighbor(i);

            int neighbour_label = std::get<3>(all_cells.find(current_cell)->second);
            int neighbour_index = std::get<0>(all_cells.find(current_cell)->second);

            if(neighbour_index > current_index)
                continue;

            Vertex_handle fv = current_cell->vertex(edge->second);
            Vertex_handle fs = current_cell->vertex(edge->third);

            int first_vertex = all_vertices.find(fv)->second;
            int second_vertex = all_vertices.find(fs)->second;

            if(current_label == 1 && neighbour_label == 1){
                fo << first_vertex << " " << second_vertex << " 0 191 255" << std::endl;
                edge_count++;
            }
            else if(current_label == 0 && neighbour_label == 0){
                fo << first_vertex << " " << second_vertex << " 255 0 0" << std::endl;
                edge_count++;
            }
        }
    }
    std::cout << "edge count is: " << edge_count << std::endl;
}

void exportPoints(std::string path, std::vector<Point>& points, std::vector<vertex_info>& infos){

    auto start = std::chrono::high_resolution_clock::now();

    path+="_fixedSensor.ply";

    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = points.size();

    std::fstream fo;
    fo.open(path, std::fstream::out);

    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << nv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property float scalar_x0" << std::endl;
    fo << "property float scalar_y0" << std::endl;
    fo << "property float scalar_z0" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(13);


    Delaunay::Finite_vertices_iterator vft;
    for(int i = 0; i < points.size(); i++){
        // print data to file
        fo << points[i] << " "                           // coordinates
           << infos[i].sensor_pos << " "
           << int(infos[i].color[0]) << " " << int(infos[i].color[1]) << " " << int(infos[i].color[2]) << std::endl;                     // normal
    }

    fo.close();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "File exported to " << path << " in " << duration.count() << "s" << std::endl;

}

void exportOFF(Polyhedron& out_mesh, std::string path)
{
    path = path + ".off";
    std::ofstream out(path);
    out << std::setprecision(15);
    out << out_mesh;
    out.close();
//    std::cout << "Exported to " << path << std::endl;
}
void exportOFF(Polyhedron_Exact& out_mesh, std::string path)
{
    path = path + ".off";
    std::ofstream out(path);
    out << std::setprecision(15);
    out << out_mesh;
    out.close();
    std::cout << "Exported to " << path << std::endl;
}
void exportOFF(Tetrahedron& in_tet, std::string path)
{
    Polyhedron out_poly;
    out_poly.make_tetrahedron(in_tet.vertex(0), in_tet.vertex(1), in_tet.vertex(2), in_tet.vertex(3));
    path = path + ".off";
    std::ofstream out(path);
    out << std::setprecision(15);
    out << out_poly;
    out.close();
//    std::cout << "Exported to " << path << std::endl;
}
void importOff(std::string path, Polyhedron& import_poly){

    std:ifstream in(path);
    in >> import_poly;

}
void importOff(std::string path, Tetrahedron& import_tet){

    Polyhedron import_poly;
    std:ifstream in(path);
    in >> import_poly;
    Polyhedron::Vertex_iterator vit;
    vit = import_poly.vertices_begin();
    import_tet = Tetrahedron(vit->point(), vit++->point(), vit++->point(), vit++->point());
}
void importOff(std::string path, std::vector<Point>& points){

    Polyhedron import_poly;
    std:ifstream in(path);
    in >> import_poly;
    Polyhedron::Vertex_iterator vit;
    for(vit = import_poly.vertices_begin(); vit != import_poly.vertices_end(); vit++){
        points.push_back(vit->point());
    }
}
void importOff(std::string path, std::vector<Plane>& planes){

    Polyhedron import_poly;
    std:ifstream in(path);
    in >> import_poly;
    Polyhedron::Vertex_iterator vit;
    std::vector<Point> points;
    for(vit = import_poly.vertices_begin(); vit != import_poly.vertices_end(); vit++){
        points.push_back(vit->point());
    }
    planes.push_back(Plane(points[0], points[2], points[1]));
    planes.push_back(Plane(points[0], points[1], points[3]));
    planes.push_back(Plane(points[1], points[2], points[3]));
    planes.push_back(Plane(points[0], points[3], points[2]));

    Point spc = CGAL::centroid(points[0], points[1], points[2], points[3]);
    for(int i = 0; i<4; i++){
        if(!planes[i].has_on_negative_side(spc)){
            planes[i]=planes[i].opposite();
            std::cout << i << " has wrong orientation!" << std::endl;
        }
    }
}


void fixSensorCenter(std::string path, std::vector<Point>& points, std::vector<vertex_info>& infos){

    auto start = std::chrono::high_resolution_clock::now();


    path=path+"_fixedNormals.ply";

    std::fstream fo;
    fo.open(path, std::fstream::out);

    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << points.size() << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property float nx" << std::endl;
    fo << "property float ny" << std::endl;
    fo << "property float nz" << std::endl;
//    fo << "property uchar red" << std::endl;
//    fo << "property uchar green" << std::endl;
//    fo << "property uchar blue" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(8);

    for(int i=0; i < points.size(); i++){

        Vector sensor(infos[i].sensor_vec[0] - points[i].x(), infos[i].sensor_vec[1] - points[i].y(), infos[i].sensor_vec[2] - points[i].z());
        sensor=sensor/std::sqrt(sensor.squared_length());
        // print data to file
        // coordinates
        fo  << points[i] << " "
       // normal
            << sensor << std::endl;
//            << infos[i].normal << std::endl;
        // color
//            << int(infos[i].color[0]) <<  " " << int(infos[i].color[1]) <<  " " << int(infos[i].color[2]) <<  std::endl;

    }
    fo.close();


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "File exported to " << path << " in " << duration.count() << "s" << std::endl;


}


void exportCellCenter(std::string path, const Delaunay& Dt){

    Delaunay::size_type nc = Dt.number_of_finite_cells();

    path+="_cellScore.ply";
    std::fstream fo;
    fo.open(path, std::fstream::out);
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << nc << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    fo << "property float inside_score" << std::endl;
    fo << "property float outside_score" << std::endl;
//    fo << "property float nz" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(8);


    Delaunay::Finite_cells_iterator cit;
    std::vector<double> inside_scores;
    std::vector<double> outside_scores;
    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        inside_scores.push_back(cit->info().inside_score);
        outside_scores.push_back(cit->info().outside_score);
    }

//    double inside_scale = 255/(*std::max_element(inside_scores.begin(), inside_scores.end()));
//    double outside_scale = 255/(*std::max_element(outside_scores.begin(), outside_scores.end()));
//    double inside_average = std::accumulate(inside_scores.begin(), inside_scores.end(), 0.0)/inside_scores.size();
//    double outside_average = std::accumulate(outside_scores.begin(), outside_scores.end(), 0.0)/outside_scores.size();
//    double inside_scale = 128/inside_average;
//    double outside_scale = 128/outside_average;
//    double inside_scale = *std::max_element(Dt.finite_cells_begin()->info().inside_score, Dt.finite_cells_end()->info().inside_score);
//    double outside_scale = *std::max_element(Dt.finite_cells_begin()->info().outside_score, Dt.finite_cells_end()->info().outside_score);

    double inside_min = *std::min_element(inside_scores.begin(), inside_scores.end());
    double outside_min = *std::min_element(outside_scores.begin(), outside_scores.end());
    double inside_max = *std::max_element(inside_scores.begin(), inside_scores.end());
    double outside_max = *std::max_element(outside_scores.begin(), outside_scores.end());

    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        Point p1 = cit->vertex(0)->point();
        Point p2 = cit->vertex(1)->point();
        Point p3 = cit->vertex(2)->point();
        Point p4 = cit->vertex(3)->point();

        Point centroid = CGAL::centroid(p1,p2,p3,p4);

        double inside_score = cit->info().inside_score;
        double outside_score = cit->info().outside_score;
        int green = 0;
        int red  = int(255*(inside_score-inside_min)/(inside_max-inside_min));
        int blue = int(255*(outside_score-outside_min)/(outside_max-outside_min));
        if(inside_score == outside_score){
            blue = 0;
            red = 0;
            green = 128;
        }

        fo << centroid << " " << red << " " << green << " " << blue << " " << inside_score << " " << outside_score << std::endl;
//        fo << centroid << " " << 0 << " " << 0 << " " << blue << " " << outside_score << std::endl;
    }


    fo.close();
    std::cout << "Exported to " << path << std::endl;

}

void exportSurfacePLY(const Delaunay& Dt,
                std::string path,
                bool optimized, bool prune_or_color)
{
    auto start = std::chrono::high_resolution_clock::now();

    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();
    Delaunay::size_type nc = Dt.number_of_finite_cells();


    // calculate how many faces to print
    Delaunay::Finite_facets_iterator fft;
    int deletedFaceCount = 0;    
    if(!optimized)
    {
        for(fft = Dt.finite_facets_begin(); fft != Dt.finite_facets_end(); fft++){
            Cell_handle c = fft->first;
            int clabel;
            float c1 = c->info().outside_score;
            float c2 = c->info().inside_score;
            if(c1 > c2)
                clabel = 1;
            else {
                clabel = 0;
            }
            Cell_handle m = Dt.mirror_facet(*fft).first;
            int mlabel;
            float m1 = m->info().outside_score;
            float m2 = m->info().inside_score;
            if(m1 > m2)
                mlabel = 1;
            else {
                mlabel = 0;
            }
            if(clabel == mlabel){deletedFaceCount++;}
        }
    }
    else{
        for(fft = Dt.finite_facets_begin(); fft != Dt.finite_facets_end(); fft++){
            Cell_handle c = fft->first;
            int clabel = c->info().final_label;
            Cell_handle m = Dt.mirror_facet(*fft).first;
            int mlabel = m->info().final_label;
            if(clabel == mlabel){deletedFaceCount++;}
        }
    }

    int sub = nf - deletedFaceCount;

    // create PLY output file for outputting the triangulation, with point coordinates, color, normals and triangle facets
    if(optimized)
        path+="_optimized";
    else
        path+="_initial";
    if(prune_or_color)
        path+="_pruned";
    else
        path+="_colored";
            
    path+=".ply";
    
    std::fstream fo;
    fo.open(path, std::fstream::out);

    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << nv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
//    if(normals){
        fo << "property float nx" << std::endl;
        fo << "property float ny" << std::endl;
        fo << "property float nz" << std::endl;
//    }
//    else
//        fo << "property int camera_index" << std::endl;
    if(prune_or_color)
        fo << "element face " << sub << std::endl;
    else
        fo << "element face " << nf << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
//    if(!prune_or_color){
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
//    }
    fo << "end_header" << std::endl;
    fo << std::setprecision(8);

    // give every vertex from the triangulation an index starting at 0
    // and already print the point coordinates, color and normal of the vertex to the PLY file
    int index = 0;
    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        // reset the vertex index here, because I need to know the order of exactly this loop here
        // for the indexing of the facets in the PLY file
        vft->info().idx = index;
        // print data to file
        // coordinates
        fo  << vft->point() << " "
        // color
            << int(vft->info().color[0]) <<  " " << int(vft->info().color[1]) <<  " " << int(vft->info().color[2]) <<  " "
        // normal
            << vft->info().sensor_vec << std::endl;
        index++;
    }

    // Save the facets to the PLY file
    int vidx;
    // initialise cell and vertex handle
    Cell_handle c;
    Vertex_handle v;
//    unsigned int deletedFaceCount = 0;
    for(fft = Dt.finite_facets_begin() ; fft != Dt.finite_facets_end() ; fft++){

        // get vertex and cell index that describes the facet
        // facet fft is represented by std::pair(cell c, int vidx). vidx is the vertex opposite to the cell.
        // even though some of the facets may be described by infinite cells, the facet is still has a neighbouring cell that is finite.
        // see: https://doc.cgal.org/latest/Triangulation_3/index.html
        c = fft->first;         // cell
        vidx = fft->second;     // vertex index

//        if(Dt.is_infinite(c))
//            std::cout << "infinite cell from finite facet" << std::endl;

        //////// check which faces to prune:



        int clabel;
        int mlabel;
        if(!optimized){
            float c1 = c->info().outside_score;
            float c2 = c->info().inside_score;
            if(c1 > c2)
                clabel = 1;
            else {
                clabel = 0;
            }
            Cell_handle m = Dt.mirror_facet(*fft).first;
            float m1 = m->info().outside_score;
            float m2 = m->info().inside_score;
            if(m1 > m2)
                mlabel = 1;
            else {
                mlabel = 0;
            }
//            // infinite stuff
//            if(Dt.is_infinite(m)){
//                mlabel == clabel;
//            }
//            else if(Dt.is_infinite(c)){
//                clabel == mlabel;
//            }

//            if(clabel==1 && Dt.is_infinite(m))
//                mlabel = 1;
//            else if(clabel==0 && Dt.is_infinite(m))
//                mlabel = 0;

        }
        else{
    //        // with GC labelling:
            clabel = c->info().final_label;
            Cell_handle m = Dt.mirror_facet(*fft).first;
            mlabel = m->info().final_label;
        }


        // check if two neighbouring cells have the same label, and if so (and the prunce_faces export function is active) continue to next face
        if(clabel == mlabel && prune_or_color){
//            std::cout << clabel << " == " << mlabel << std::endl;
            continue;
        }
//        else{
//            std::cout << clabel << " != " << mlabel << std::endl;
//        }

        // if label of neighbouring cells is not the same...
        // if opposite vertex vidx is 2, we start at j = vidx + 1 = 3, 3%4 = 3
        // next iteration: j = 4, 4%4 = 0, next iteration: j = 5, 5%4 = 1;
        // so we exactely skip 2 - the opposite vertex.
        // start printed facet line with a 3
        fo << 3 << ' ';
        // fix the orientation
        std::vector<Vertex_handle> tri;
        for(int j = vidx + 1; j <= vidx + 3; j++){
            tri.push_back(c->vertex(j%4));
//            Point p = c->vertex(j%4)->point();
        }
        // add up all the sensor positions
        Point avSensor = CGAL::centroid(tri[0]->info().sensor_pos,
                                        tri[1]->info().sensor_pos,
                                        tri[2]->info().sensor_pos);
        // check if sensor position is on the positive side of the triangle
        // otherwise change order
        Plane p(tri[0]->point(), tri[1]->point(), tri[2]->point());
        if(p.has_on_negative_side(avSensor)){
//        if(p.has_on_positive_side(avSensor)){
            Vertex_handle temp = tri[1];
            tri[1] = tri[2];
            tri[2] = temp;
        }
        for(int j = 0; j < 3; j++){
            // print the indicies of each cell to the file
            fo << tri[j]->info().idx << ' ';
        }

        // color the facets (if !prune_faces, meaning coloring is active)
        if(clabel == 1 && mlabel == 1){
            fo << " 0 191 255";
        }
        else if(clabel == 0 && mlabel == 0){
            fo << "255 0 0";
        }
        else{
            fo << "0 255 0";
        }

        fo << std::endl;
    }

    fo.close();

    std::cout << "before face count: " << nf << std::endl;
    std::cout << "remaining faces: " << sub << std::endl;
    std::cout << "Exported to " << path << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "...in " << duration.count() << "s" << std::endl;
}
