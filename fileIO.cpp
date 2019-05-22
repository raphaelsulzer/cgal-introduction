#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/intersections.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Euclidean_distance.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;

// vertext base for point + info (=vector, color, intensity)
typedef Kernel::Vector_3                                            Vector;
typedef CGAL::cpp11::array<unsigned char, 3>                        Color;
//typedef std::tuple<Vector, Color>                                   VColb;
// vertext base for point + info (=vector)
typedef CGAL::Triangulation_vertex_base_with_info_3<Vector, Kernel> VNb;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, Kernel> VCb;


typedef CGAL::Delaunay_triangulation_cell_base_3<Kernel>            Cb;         // cell base
typedef CGAL::Triangulation_data_structure_3<VNb, Cb>               Tds;        // triangulation data structure
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>                 Delaunay;   // delaunay triangulation based on triangulation data structure
typedef Delaunay::Point                                             Point;
typedef Delaunay::Edge                                              Edge;
typedef Delaunay::Facet                                             Facet;
typedef Delaunay::Cell_handle                                       Cell_handle;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
// map cell of the Dt, to an index, the outside score, the inside score, and the final label
typedef std::map<Cell_handle, std::tuple<int, float, float, int>>   Cell_map;
typedef std::map<Vertex_handle, int>                                Vertex_map;

// for reading PLY file
typedef CGAL::cpp11::tuple<Point, int> PC;
typedef CGAL::Nth_of_tuple_property_map<0, PC> PointC_map;
typedef CGAL::Nth_of_tuple_property_map<1, PC> Camera_map;

typedef CGAL::cpp11::tuple<Point, Vector> PN;
typedef CGAL::Nth_of_tuple_property_map<0, PN> PointN_map;
typedef CGAL::Nth_of_tuple_property_map<1, PN> Normal_map;




//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
//////////////////////////////////////////////////////////
std::vector<Point> readPlyWithO(std::string fname)
{
    // Reads a .ply point set file
    std::vector<Point> points; // store points
    std::ifstream in(fname);
    CGAL::read_ply_points(in, std::back_inserter (points));

    std::cout << "PLY file read..." << std::endl;
    return points;
}

std::vector<PN> readPlyWithN(std::string fname)
{
    std::vector<PN> points; // store points
    std::ifstream in(fname);

    CGAL::read_ply_points_with_properties
      (in,
       std::back_inserter (points),
       CGAL::make_ply_point_reader (PointN_map()),
       CGAL::make_ply_normal_reader (Normal_map()));

    std::cout << "PLY file with normals read..." << std::endl;
    return points;
}

std::vector<PC> readPlyWithC(std::string fname)
{
    std::vector<PC> points; // store points
    std::ifstream in(fname);

    CGAL::read_ply_points_with_properties
      (in,
       std::back_inserter (points),
       CGAL::make_ply_point_reader (PointC_map()),
       std::make_pair (Camera_map(), CGAL::PLY_property<int>("camera_index")));

    std::cout << "PLY file with camera index read..." << std::endl;
    return points;
}

//// generate a Delaunay triangulation from a PLY file
//Delaunay triangulationFromFile(const char* ifn)
//{
//   // get data as vector of tuples(point, normal, color, intensity)
//   auto ply = readPlyWithCnFun(ifn);

//   std::vector<Point> points;
//   std::vector<Vector> infos;
//   for (std::size_t i = 0; i < ply.size (); ++ i)
//   {
//       // make vector of points
//       points.push_back(get<0>(ply[i]));
////       // make vector of infos as: tuple(normal, color, intensity)
////       infos.push_back(std::make_tuple(get<1>(ply[i]), get<2>(ply[i]), get<3>(ply[i])));
//       // make vector of infos as: tuple(normal, color, intensity)
//       infos.push_back(get<1>(ply[i]));
//   }

//   // make the triangulation
//   Delaunay Dt( boost::make_zip_iterator(boost::make_tuple( points.begin(),infos.begin() )),
//             boost::make_zip_iterator(boost::make_tuple( points.end(),infos.end() ) )  );
//   std::cout << "Triangulation done.." << std::endl;
//   return Dt;
//}

// generate a Delaunay triangulation from a PLY file
//Delaunay triangulationFromFile(const char* ifn, const char* option)
Delaunay triangulationFromFile(std::string ifn)
{

//    if(std::strcmp(option,"O"))
//        std::vector<Point> ply = readPlyWithO(ifn);
//    else if(std::strcmp(option,"N"))
//        std::vector<PN> ply = readPlyWithN(ifn);
//    else if(std::strcmp(option,"C"))
//        std::vector<PC> ply = readPlyWithC(ifn);

//    std::vector<Point> ply = readPlyWithO(ifn);

    std::vector<PN> ply = readPlyWithN(ifn);
    std::vector<Vector> infos;

//    std::vector<PC> ply = readPlyWithC(ifn);
//    std::vector<int> infos;

    std::vector<Point> points;
    for (std::size_t i = 0; i < ply.size (); ++ i)
    {
        // make vector of points
        points.push_back(get<0>(ply[i]));
        //       // make vector of infos as: tuple(normal, color, intensity)
        //       infos.push_back(std::make_tuple(get<1>(ply[i]), get<2>(ply[i]), get<3>(ply[i])));
        // make vector of infos as: tuple(normal, color, intensity)
        infos.push_back(get<1>(ply[i]));
    }

    // make the triangulation
    Delaunay Dt( boost::make_zip_iterator(boost::make_tuple( points.begin(),infos.begin() )),
    boost::make_zip_iterator(boost::make_tuple( points.end(),infos.end() ) )  );
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



void exportSoup(const Delaunay& Dt, Cell_map& all_cells, std::string path, bool optimized, bool prune_faces)
{
    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();

    // calculate how many faces to print
    Delaunay::Finite_facets_iterator fft;
    int deletedFaceCount = 0;
    for(fft = Dt.finite_facets_begin(); fft != Dt.finite_facets_end(); fft++){
        Cell_handle c = fft->first;
        int clabel = std::get<3>(all_cells.find(c)->second);
        Cell_handle m = Dt.mirror_facet(*fft).first;
        int mlabel = std::get<3>(all_cells.find(m)->second);
        if(clabel == mlabel){deletedFaceCount++;}
    }
    int sub = nf - deletedFaceCount;

    // create PLY output file for outputting the triangulation, with point coordinates, color, normals and triangle facets
    if(optimized)
        path+="_optimized";
    else
        path+="_initial";
    if(prune_faces)
        path+="_pruned";
    else
        path+="_colored";
            
    path+=".ply";
                
    
    std::fstream fo;
    fo.open(path, std::fstream::out);

    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "comment VCGLIB generated" << std::endl;
    fo << "element vertex " << nv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
//    fo << "property float nx" << std::endl;
//    fo << "property float ny" << std::endl;
//    fo << "property float nz" << std::endl;
    fo << "property int camera_index" << std::endl;
    if(prune_faces)
        fo << "element face " << sub << std::endl;
    else
        fo << "element face " << nf << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
    if(!prune_faces){
        fo << "property uchar red" << std::endl;
        fo << "property uchar green" << std::endl;
        fo << "property uchar blue" << std::endl;
    }
    fo << "end_header" << std::endl;
    fo << std::setprecision(3);

    // give every vertex from the triangulation an index starting at 0
    // and already print the point coordinates, color and normal of the vertex to the PLY file
    Vertex_map Vertices;
    int index = 0;
    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        Vertices[vft] = index;
        // print data to file
        fo << vft->point() << " "                           // coordinates
           << vft->info() << std::endl;                     // normal
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
        int clabel = std::get<3>(all_cells.find(c)->second);
        Cell_handle m = Dt.mirror_facet(*fft).first;
        int mlabel = std::get<3>(all_cells.find(m)->second);


        if(!optimized){
            float c1 = std::get<1>(all_cells.find(c)->second);
            float c2 = std::get<2>(all_cells.find(c)->second);
            int clabel;
            if(c1 > c2)
                clabel = 0;
            else {
                clabel = 1;
            }
            Cell_handle m = Dt.mirror_facet(*fft).first;
            float m1 = std::get<1>(all_cells.find(m)->second);
            float m2 = std::get<2>(all_cells.find(m)->second);
            int mlabel;
            if(m1 > m2)
                clabel = 0;
            else {
                clabel = 1;
            }
        }



        if(clabel == mlabel && prune_faces){
            continue;
        }

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


        if(clabel == 1 && mlabel == 1 && !prune_faces){
            fo << " 0 191 255";
        }
        else if(clabel == 0 && mlabel == 0 && !prune_faces){
            fo << "255 0 0";
        }
        else if(!prune_faces){
            fo << "0 255 0";
        }

        fo << std::endl;
    }

//    fo.clear();
//    fo.seekg(170,std::ios_base::beg);

//    // rewrite the number of faces line in the file, by substracting the number of pruned faces
//    int sub = nf - deletedFaceCount;
//    unsigned int currentLine = 0;
//    while ( currentLine < 10 )
//    {
//         fo.ignore( std::numeric_limits<std::streamsize>::max(), '\n');
//         ++currentLine;
//    }
//    fo << "element face " << sub << std::endl;

//    exportEdges(fo, Dt, all_cells, Vertices);

    fo.close();

    std::cout << "before face count: " << nf << std::endl;
    std::cout << "remaining faces: " << sub << std::endl;
    std::cout << "Delaunay triangulation done and exported to PLY file!" << std::endl;
}
