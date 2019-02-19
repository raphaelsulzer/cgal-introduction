#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <vector>

#include <readPlyWithCn.cpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;
//typedef Kernel::FT FT;
typedef Kernel::Vector_3                                            Vector;
typedef CGAL::cpp11::array<unsigned char, 3>                        Color;

typedef std::tuple<Vector, Color, int>                              In;
typedef CGAL::Triangulation_vertex_base_with_info_3<In, Kernel>     Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<Kernel>            Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>                 Delaunay;
typedef Delaunay::Point                                             Point;
typedef Delaunay::Cell_handle                                       Cell_handle;
typedef Delaunay::Vertex_handle                                     Vertex_handle;
typedef CGAL::cpp11::tuple<Point, Vector, Color, int>               PNCI;


int exportTriWithCnFun(const char* ifn, const char* ofn)
{

    // get points with colors and normals
    auto ply = readPlyWithCnFun(ifn);

    std::vector<Point> points;
    std::vector<In> infos;
    for (std::size_t i = 0; i < ply.size (); ++ i)
    {
        points.push_back(get<0>(ply[i]));
        // make info as: tuple(normal, color, intensity)
        infos.push_back(std::make_tuple(get<1>(ply[i]), get<2>(ply[i]), get<3>(ply[i])));
    }

    Delaunay Dt( boost::make_zip_iterator(boost::make_tuple( points.begin(),infos.begin() )),
              boost::make_zip_iterator(boost::make_tuple( points.end(),infos.end() ) )  );

    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();

    std::cout << "Triangulation done.." << std::endl;

    std::ofstream fo;
    fo.open (ofn);

    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "comment VCGLIB generated" << std::endl;
    fo << "element vertex " << nv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    fo << "property float nx" << std::endl;
    fo << "property float ny" << std::endl;
    fo << "property float nz" << std::endl;
    fo << "element face " << nf << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(3);

    // give every vertex from the triangulation an index starting at 0
    // and already print the point coordinates of the vertex to the PLY file
    std::map<Vertex_handle, int> Vertices;
    int index = 0;
    Delaunay::Finite_vertices_iterator vft;
    //std::cout << "iterate over finite vertices: " << std::endl;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        Vertices[vft] = index;  // assign index to vertex handle
        //std::cout << "index: " <<index << ": " << vft->point() << std::endl;
        fo << vft->point() << " "
           << int(std::get<1>(vft->info())[0]) << " "
           << int(std::get<1>(vft->info())[1]) << " "
           << int(std::get<1>(vft->info())[2]) << " "
           << std::get<0>(vft->info()) << std::endl;
        index++;
    }

    // Save the facets to the PLY file
    //std::cout << "iterate over finite triangles: " << std::endl;
    Delaunay::Finite_facets_iterator fft;
    int vidx;
    // initialise cell and vertex handle
    Cell_handle c;
    Vertex_handle v;
    for(fft = Dt.finite_facets_begin() ; fft != Dt.finite_facets_end() ; fft++){

        // facet fft is represented by std::pair(cell c, int vidx). vidx is the vertex opposite to the cell.
        // even though some of the facets may be described by infinite cells, the facet is still has a neighbouring cell that is finite.
        // see: https://doc.cgal.org/latest/Triangulation_3/index.html


        c = fft->first;         // cell
        vidx = fft->second;     // vertex index


        //std::cout << "is infinite: " << T.is_infinite(v) << std::endl;
        //std::cout << "opposite vertex: " << vidx << std::endl;
        fo << 3 << ' ';
        // if opposite vertex vidx is 2, we start at j = vidx + 1 = 3, 3%4 = 3
        // next iteration: j = 4, 4%4 = 0, next iteration: j = 5, 5%4 = 1;
        // so we exactely skip 2 - the opposite vertex.
        for(int j = vidx + 1 ; j <= vidx + 3 ; j++){


            //std::cout << "modulo: " << j%4;

            // in the first and second iteration I am calling vertex() with the same value for j%4,
            // but get a different vertex idx. Reason is that the cell c is different (see std::cout of cell points).
            // for some reason the first facet is seen from the infinite cell.
            v = c->vertex(j%4);

            //std::cout << ' ' << "idx: " << Vertices.find(v)->second << ' ';
            fo << Vertices.find(v)->second << ' ';

        }
        //std::cout << ' ' << std::endl;
        fo << std::endl;
    }
    fo.close();
    return 0;

}



