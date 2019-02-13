#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/make_surface_mesh.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include "readPLY.cpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K>       Triangulation;
typedef Triangulation::Cell_handle              Cell_handle;
typedef Triangulation::Vertex_handle            Vertex_handle;
typedef Triangulation::Locate_type              Locate_type;
typedef Triangulation::Point                    Point;
typedef Triangulation::Facet                    Facet;

int main()
{
    
    
    auto filename = "/Users/Raphael/Library/Mobile Documents/com~apple~CloudDocs/Studium/PhD/Paris/data/EMS/Est1 - Cloud_subsampled.ply";
    auto L = readPLYfun(filename);
    
    
//    // construction from a list of points :
//    std::list<Point> L;
//
//
//
//
//    L.push_front(Point(10,10,10));
//    L.push_front(Point(11,10,10));
//    L.push_front(Point(10,11,10));
//    L.push_front(Point(10,10,11));
//    L.push_front(Point(10,14,11));
    //    L.push_front(Point(11,10,13));
    //    L.push_front(Point(22,13,11));
    //    L.push_front(Point(10,31,41));
    //    L.push_front(Point(12,10,17));
    
    Triangulation T(L.begin(), L.end());
    Triangulation::size_type nv = T.number_of_vertices();
    Triangulation::size_type nf = T.number_of_finite_facets();
    
    std::cout << "Triangulation done.." << std::endl;
    
    std::ofstream fo;
    fo.open ("/Users/Raphael/Library/Mobile Documents/com~apple~CloudDocs/Studium/PhD/Paris/data/triangulation.ply");
    
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "comment VCGLIB generated" << std::endl;
    fo << "element vertex " << nv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "element face " << nf << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(12);
    
    std::cout << nv << std::endl;
    std::cout << nf << std::endl;
    std::cout << "number of facets: " << nf << std::endl;
    
    
    Triangulation::Finite_facets_iterator fft;
    Triangulation::Finite_vertices_iterator vft;
    int vidx;
    
    
    Cell_handle c;
    Vertex_handle v;
    
    std::map<Vertex_handle, int> Vertices;
    int index = 0;
    
    std::cout << "iterate over finite vertices: " << std::endl;
    for (vft = T.finite_vertices_begin() ; vft != T.finite_vertices_end() ; vft++){
        
        Vertices[vft] = index;
        //std::cout << "index: " <<index << ": " << vft->point() << std::endl;
        fo << vft->point() << std::endl;
        index++;
    }
    
    //std::cout << "iterate over finite triangles: " << std::endl;
    for(fft = T.finite_facets_begin() ; fft != T.finite_facets_end() ; fft++){
        
        // facet fft is represented by pair (cell c, int vidx). vidx is the vertex opposite to the cell.
        // even though some of the facets may be described by infinite cells, the facet is still has a neighbouring cell that is finite. See: https://doc.cgal.org/latest/Triangulation_3/index.html
        
        c = fft->first;
        vidx = fft->second;
        
        //std::cout << "cell points: ";
        for(int k = 0 ; k <= 3 ; k++){
            
            //std::cout << c->vertex(k)->point() << "; ";
        }
        //std::cout << std::endl;
        
        
        
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
