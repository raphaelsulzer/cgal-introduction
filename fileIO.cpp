#include <cgal_typedefs.h>
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



void readBinaryPLY(std::string ifn, std::vector<Point>& points, std::vector<vertex_info>& infos, bool colmap)
{
    auto start = std::chrono::high_resolution_clock::now();

    // read Binary PLY with sensor
    Mesh_ply aMesh;
    Import_PLY(ifn.c_str(), &aMesh);

    if(colmap){

        std::string sensor_file = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/fused.ply.vis";
        std::vector<std::vector<int>> vis_info;
        loadVisibilityFile(sensor_file.c_str(), vis_info);

        std::map<int, Point> sensor_map;
        std::string image_file = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/sfm_result/images.bin";
        loadImageFile(image_file.c_str(), sensor_map);

        for(int i = 0; i < aMesh.mVertices.size(); i++){
            // save point
            Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
            points.push_back(pt);
            // sensor
            // get the first of the cameras of this point
            Point sensor_location = sensor_map.find(vis_info[i][0])->second;
    //        Vector vec(sensor_location.x() - pt.x(), aMesh.mvCapture[i].y - pt.y(), aMesh.mvCapture[i].z - pt.z());
            Vector vec = sensor_location - pt;
            vertex_info vec_inf;
            vec_inf.sensor = vec;
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
    }
    else{
        for(int i = 0; i < aMesh.mVertices.size(); i++){
            // save points
            Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
            points.push_back(pt);
            // sensor
//            Vector vec(aMesh.mNormals[i].x - pt.x(), aMesh.mNormals[i].y - pt.y(), aMesh.mNormals[i].z - pt.z());
            Vector vec(aMesh.mvCapture[i].x - pt.x(), aMesh.mvCapture[i].y - pt.y(), aMesh.mvCapture[i].z - pt.z());
            vertex_info vec_inf;
            vec_inf.sensor = vec;
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
    }





//    std::vector<std::vector<int>> sensor_triangle;
//    for(int i = 0; i < aMesh.mIndices.size()/3; i++){
//        std::vector<int> poly(3);
//        poly[0] = aMesh.mIndices[(i*3)%3+0];
//        poly[1] = aMesh.mIndices[(i*3)%3+1];
//        poly[2] = aMesh.mIndices[(i*3)%3+2];
//        sensor_triangle.push_back(poly);
//    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "File " << ifn << " read in " << duration.count() << "s" << std::endl;


}


//void readPLYWithSensor(std::string ifn, std::vector<Point>& points, std::vector<vertex_info>& infos)
//{

//    auto start = std::chrono::high_resolution_clock::now();

//    std::vector<PS> ply_lines;

//    //    std::vector<PC> points; // store points
//    std::ifstream in(ifn);

//    // read the ply
//    typedef CGAL::Nth_of_tuple_property_map<0, PS> Point_map;
//    typedef CGAL::Nth_of_tuple_property_map<1, PS> Sensor_map;
//    CGAL::read_ply_points_with_properties
//      (in,
//       std::back_inserter (ply_lines),
//       CGAL::make_ply_point_reader (Point_map()),
//       std::make_tuple (Sensor_map(),
//                               CGAL::Construct_array(),
//                               CGAL::PLY_property<float>("x0"),
//                               CGAL::PLY_property<float>("y0"),
//                               CGAL::PLY_property<float>("z0"))
//       );


//    for (int i = 0; i < ply_lines.size (); ++ i)
//    {
//        // make vector of points
//        points.push_back(get<0>(ply_lines[i]));
//        // make vector of infos as: tuple(idx, sigma, normal, color)
//        vertex_info inf;
//        inf.idx = i;
//        inf.sensor = get<1>(ply_lines[i]);
//        inf.sigma = 0.0;
//        infos.push_back(inf);
//    }


//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    std::cout << "File " << ifn << " read in " << duration.count() << "s" << std::endl;
//}


void readPLY(std::string ifn, std::vector<Point>& points, std::vector<vertex_info>& infos)
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
        inf.sensor = get<1>(ply_lines[i]);
        inf.sigma = 0.0;
        infos.push_back(inf);
    }


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "File " << ifn << " read in " << duration.count() << "s" << std::endl;
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
    fo << "property float nx" << std::endl;
    fo << "property float ny" << std::endl;
    fo << "property float nz" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(8);


    Delaunay::Finite_vertices_iterator vft;
    for(int i = 0; i < points.size(); i++){
        // print data to file
        fo << points[i] << " "                           // coordinates
           << infos[i].sensor << " "
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
    out << out_mesh;
    out.close();
    std::cout << "Exported to " << path << std::endl;
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

        Vector sensor(infos[i].sensor[0] - points[i].x(), infos[i].sensor[1] - points[i].y(), infos[i].sensor[2] - points[i].z());
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


//void exportCellCenter(std::string path, const Delaunay& Dt){

//    Delaunay::size_type nc = Dt.number_of_finite_cells();

//    path+="_cellScore.ply";
//    std::fstream fo;
//    fo.open(path, std::fstream::out);
//    fo << "ply" << std::endl;
//    fo << "format ascii 1.0" << std::endl;
//    fo << "element vertex " << nc << std::endl;
//    fo << "property float x" << std::endl;
//    fo << "property float y" << std::endl;
//    fo << "property float z" << std::endl;
//    fo << "property float inside" << std::endl;
//    fo << "property float outside" << std::endl;
//    fo << "end_header" << std::endl;
//    fo << std::setprecision(8);


//    Delaunay::Finite_cells_iterator cft;
////    std::vector<double> inside_scores;
////    std::vector<double> outside_scores;
////    for(cft = Dt.finite_cells_begin(); cft!=Dt.finite_cells_end(); cft++)
////    {
////        inside_scores.push_back(cft->info().inside_score);
////        outside_scores.push_back(cft->info().outside_score);
////    }

////    double inside_max = *std::max_element(inside_scores.begin(), inside_scores.end());
////    double outside_max = *std::max_element(outside_scores.begin(), outside_scores.end());

////    double inside_scale = 255.0/inside_max;
////    double outside_scale = 255.0/outside_max;

//    // iterate over the all_cells map
//    for(cft = Dt.finite_cells_begin(); cft!=Dt.finite_cells_end(); cft++)
//    {

//        Point p1 = cft->vertex(0)->point();
//        Point p2 = cft->vertex(1)->point();
//        Point p3 = cft->vertex(2)->point();
//        Point p4 = cft->vertex(3)->point();

//        Point centroid = CGAL::centroid(p1,p2,p3,p4);

//        double inside = cft->info().inside_score;
//        double outside = cft->info().outside_score;

////        std::cout << centroid << " " << int(red*inside_scale) << " " << 0 <<  " " << int(blue*outside_scale) << std::endl;

////        fo  << centroid << " " << int(red*inside_scale) << " " << 0 <<  " " << int(blue*outside_scale) << std::endl;
//        fo  << centroid << " " << inside << " " << outside << std::endl;
//    }
//    fo.close();

//}


void exportCellCenter(std::fstream& fo, const Delaunay& Dt){


    Delaunay::Finite_cells_iterator cit;

    std::vector<double> inside_scores;
    std::vector<double> outside_scores;
    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        inside_scores.push_back(cit->info().inside_score);
        outside_scores.push_back(cit->info().outside_score);
    }

    double inside_scale = 255/(*std::max_element(inside_scores.begin(), inside_scores.end()));
    double outside_scale = 255/(*std::max_element(outside_scores.begin(), outside_scores.end()));
//    double inside_scale = *std::max_element(Dt.finite_cells_begin()->info().inside_score, Dt.finite_cells_end()->info().inside_score);
//    double outside_scale = *std::max_element(Dt.finite_cells_begin()->info().outside_score, Dt.finite_cells_end()->info().outside_score);

    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        Point p1 = cit->vertex(0)->point();
        Point p2 = cit->vertex(1)->point();
        Point p3 = cit->vertex(2)->point();
        Point p4 = cit->vertex(3)->point();

        Point centroid = CGAL::centroid(p1,p2,p3,p4);

        double inside_score = cit->info().inside_score;
        double outside_score = cit->info().outside_score;

        fo << centroid << " " << int(inside_score*inside_scale) << " 0 " << int(outside_score*outside_scale) << " 0 0 0" << std::endl;
    }

}






void exportPLY(const Delaunay& Dt,
                std::string path,
                bool normals, bool optimized, bool prune_or_color)
{
    auto start = std::chrono::high_resolution_clock::now();

    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();
    Delaunay::size_type nc = Dt.number_of_finite_cells();


    // calculate how many faces to print
    Delaunay::Finite_facets_iterator fft;
    int deletedFaceCount = 0;
    for(fft = Dt.finite_facets_begin(); fft != Dt.finite_facets_end(); fft++){
        Cell_handle c = fft->first;
        int clabel = c->info().final_label;
        Cell_handle m = Dt.mirror_facet(*fft).first;
        int mlabel = m->info().final_label;
        if(clabel == mlabel){deletedFaceCount++;}
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
    fo << "element vertex " << nv+nc << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    if(normals){
        fo << "property float nx" << std::endl;
        fo << "property float ny" << std::endl;
        fo << "property float nz" << std::endl;
    }
    else
        fo << "property int camera_index" << std::endl;
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
            << vft->info().sensor << std::endl;
        index++;
    }


    // export cell center
    exportCellCenter(fo, Dt);




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

        //////// check which faces to prune:

//        // with GC labelling:
        int clabel = c->info().final_label;
        Cell_handle m = Dt.mirror_facet(*fft).first;
        int mlabel = m->info().final_label;


        if(!optimized){
            float c1 = c->info().outside_score;
            float c2 = c->info().inside_score;
            int clabel;
            if(c1 > c2)
                clabel = 1;
            else {
                clabel = 0;
            }
            Cell_handle m = Dt.mirror_facet(*fft).first;
            float m1 = m->info().outside_score;
            float m2 = m->info().inside_score;
            int mlabel;
            if(m1 > m2)
                mlabel = 1;
            else {
                mlabel = 0;
            }
        }


        // check if two neighbouring cells have the same label, and if so (and the prunce_faces export function is active) continue to next face
        if(clabel == mlabel && prune_or_color){
            continue;
        }

        // if label of neighbouring cells is not the same...
        // start printed facet line with a 3
        fo << 3 << ' ';
        // if opposite vertex vidx is 2, we start at j = vidx + 1 = 3, 3%4 = 3
        // next iteration: j = 4, 4%4 = 0, next iteration: j = 5, 5%4 = 1;
        // so we exactely skip 2 - the opposite vertex.
        for(int j = vidx + 1 ; j <= vidx + 3 ; j++){
            // print the indicies of each cell to the file
            // vertices is a map of all vertices of the triangulation to an index
            v = c->vertex(j%4);
            fo << v->info().idx << ' ';
        }


//        // color the facets (if !prune_faces, meaning coloring is active)
//        if(clabel == 1 && mlabel == 1 && !prune_or_color){
//            fo << " 0 191 255";
//        }
//        else if(clabel == 0 && mlabel == 0 && !prune_or_color){
//            fo << "255 0 0";
//        }
//        else if(!prune_or_color){
//            fo << "0 255 0";
//        }

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
