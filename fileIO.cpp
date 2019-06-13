#include <cgal_typedefs.h>


//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
//////////////////////////////////////////////////////////
// generate a Delaunay triangulation from a PLY file
Delaunay triangulationFromFile(std::string ifn, std::vector<PNC> ply_lines)
{

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
    std::vector<IdxSigNormCol> infos;
    for (int i = 0; i < ply_lines.size (); ++ i)
    {
        // make vector of points
        points.push_back(get<0>(ply_lines[i]));
        // make vector of infos as: tuple(idx, sigma, normal, color)
        infos.push_back(std::make_tuple(i, 0.0, get<1>(ply_lines[i]), get<2>(ply_lines[i])));
    }

    // make the triangulation
    Delaunay Dt( boost::make_zip_iterator(boost::make_tuple(points.begin(), infos.begin() )),
    boost::make_zip_iterator(boost::make_tuple(points.end(), infos.end() ) )  );
    std::cout << "Triangulation done.." << std::endl;
    return Dt;
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

void exportSimple(const Delaunay& Dt, std::map<Vertex_handle, std::pair<Point, double>>& all_vertices, std::string path){

    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();

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
    fo << "property float curvature" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(3);


    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        Point normal = all_vertices.find(vft)->second.first;
        double curvature = all_vertices.find(vft)->second.second;
        // print data to file
        fo << vft->point() << " "                           // coordinates
           << normal << " "
           << curvature << std::endl;                     // normal
    }

    fo.close();

}

void exportOFF(Polyhedron& out_mesh, std::string path)
{

    path = path + ".off";
    std::ofstream out(path);
    out << out_mesh;
    out.close();
    std::cout << "Exported to " << path << std::endl;


}

void exportPLY(const Delaunay& Dt,
                std::string path,
                bool normals, bool optimized, bool prune_or_color)
{
    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();

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
    fo << "element vertex " << nv << std::endl;
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
    if(!prune_or_color){
        fo << "property uchar red" << std::endl;
        fo << "property uchar green" << std::endl;
        fo << "property uchar blue" << std::endl;
    }
    fo << "end_header" << std::endl;
    fo << std::setprecision(3);

    // give every vertex from the triangulation an index starting at 0
    // and already print the point coordinates, color and normal of the vertex to the PLY file
    int index = 0;
    Vertex_map Vertices;
    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        Vertices[vft] = index;
        // print data to file
        fo  << vft->point() << " "                           // coordinates
            << int(get<3>(vft->info())[0]) << " " << int(get<3>(vft->info())[1]) << " " << int(get<3>(vft->info())[2]) << " "             // color
            << get<2>(vft->info()) << std::endl;                    // normal

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
            fo << Vertices.find(v)->second << ' ';
        }


        // color the facets (if !prune_faces, meaning coloring is active)
        if(clabel == 1 && mlabel == 1 && !prune_or_color){
            fo << " 0 191 255";
        }
        else if(clabel == 0 && mlabel == 0 && !prune_or_color){
            fo << "255 0 0";
        }
        else if(!prune_or_color){
            fo << "0 255 0";
        }

        fo << std::endl;

    }

    fo.close();


    std::cout << "before face count: " << nf << std::endl;
    std::cout << "remaining faces: " << sub << std::endl;
    std::cout << "Exported to " << path << std::endl;
}
