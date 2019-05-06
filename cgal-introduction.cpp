#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <vector>

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/intersections.h>

#include <CGAL/IO/read_ply_points.h>


//#include <readPlyWithCn.cpp>
//#include <exportTri.cpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;

// vertext base for point + info (=vector, color, intensity)
typedef Kernel::Vector_3                                            Vector;
typedef CGAL::cpp11::array<unsigned char, 3>                        Color;
//typedef std::tuple<Vector, Color>                                   VC;
//typedef CGAL::Triangulation_vertex_base_with_info_3<VC, Kernel>     Vb;
// vertext base for point + info (=vector)
typedef CGAL::Triangulation_vertex_base_with_info_3<Vector, Kernel> Vb;


typedef CGAL::Delaunay_triangulation_cell_base_3<Kernel>            Cb;         // cell base
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;        // triangulation data structure
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>                 Delaunay;   // delaunay triangulation based on triangulation data structure
typedef Delaunay::Point                                             Point;
typedef Delaunay::Facet                                             Facet;
typedef Delaunay::Cell_handle                                       Cell_handle;
typedef Delaunay::Vertex_handle                                     Vertex_handle;


typedef Kernel::Ray_3 Ray;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Intersect_3 Intersect;
typedef Kernel::Segment_3 Segment;

typedef std::pair<Point, Vector> PointVectorPair;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// generate a simple point set as an example
std::vector<Point> makeSimplePointSet()
{
    std::vector<Point> L(18);
    L[0]=Point(0,2,2);
    L[1]=Point(0,7,2);
    L[2]=Point(0,7,4);
    L[3]=Point(0,4,4);
    L[4]=Point(0,4,7);
    L[5]=Point(0,2,7);

    L[6]=Point(3,2,2);
    L[7]=Point(3,7,2);
    L[8]=Point(3,7,4);
    L[9]=Point(3,4,4);
    L[10]=Point(3,4,7);
    L[11]=Point(3,2,7);

    L[12]=Point(1.5,5.5,4);
    L[13]=Point(1.5,4,5.5);
    L[14]=Point(0,4.5,3);
    L[15]=Point(0,3,4.5);
    L[16]=Point(3,4.5,3);
    L[17]=Point(3,3,4.5);

    return L;
}

// estimate normals of a point set
std::vector<PointVectorPair> estimateNormalsFun(const std::vector<Point>& points)
{

    // initialise a vector of point-vector-pairs
    std::vector<PointVectorPair> pointVectorPairs(points.size());

    // add the points as the first element of the point vector pair
    for(std::size_t i=0; i < points.size(); ++i)
    {
        pointVectorPairs[i].first = points[i];
    }

    // following two blocks are from this CGAL example: https://doc.cgal.org/latest/Point_set_processing_3/Point_set_processing_3_2normals_example_8cpp-example.html
    // Estimates normals direction.
    // Note: pca_estimate_normals() requires a range of points
    // as well as property maps to access each point's position and normal.
    const int nb_neighbors = 4;
    CGAL::pca_estimate_normals<Concurrency_tag>
      (pointVectorPairs, nb_neighbors,
       CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
       normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

    // Orients normals.
    // Note: mst_orient_normals() requires a range of points
    // as well as property maps to access each point's position and normal.
    std::vector<PointVectorPair>::iterator unoriented_points_begin =    // this returns unoriented points that could be removed in the next step
            CGAL::mst_orient_normals(pointVectorPairs, nb_neighbors,
                               CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                               normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

    std::cout << "Normals calculated!" << std::endl;
    return pointVectorPairs;
};


// // generate a Delaunay triangulation from a PLY file
//Delaunay triangulationFromFile(const char* ifn)
//{
//    // get data as vector of tuples(point, normal, color, intensity)
//    auto ply = readPlyWithCnFun(ifn);

//    std::vector<Point> points;
//    std::vector<IVCI> infos;
//    for (std::size_t i = 0; i < ply.size (); ++ i)
//    {
//        // make vector of points
//        points.push_back(get<0>(ply[i]));
//        // make vector of infos as: tuple(normal, color, intensity)
//        infos.push_back(std::make_tuple(get<1>(ply[i]), get<2>(ply[i]), get<3>(ply[i])));
//    }

//    // make the triangulation
//    Delaunay Dt( boost::make_zip_iterator(boost::make_tuple( points.begin(),infos.begin() )),
//              boost::make_zip_iterator(boost::make_tuple( points.end(),infos.end() ) )  );
//    std::cout << "Triangulation done.." << std::endl;


//    return Dt;

//}



// generate a simple Delaunay triangulation
Delaunay triangulationSimple()
{
    // get data as vector of tuples(point, normal, color, intensity)
    std::vector<Point> points = makeSimplePointSet();

    // estimate normals on the point set
    std::vector<PointVectorPair> pVP = estimateNormalsFun(points);

    // make the triangulation
    Delaunay Dt(pVP.begin(), pVP.end());
    std::cout << "Triangulation done.." << std::endl;

    return Dt;

}



//Delaunay triangulatePNC(std::vector<Point>)
//{
//    // get data as vector of tuples(point, normal, color, intensity)
//    std::vector<Point> points = makeSimplePointSet();

//    // estimate normals on the point set
//    std::vector<PointVectorPair> pVP = estimateNormalsFun(points);

//    // make the triangulation
//    Delaunay Dt(pVP.begin(), pVP.end());
//    std::cout << "Triangulation done.." << std::endl;

//    return Dt;

//}



//int exportTriWithCnFun(const char* ifn, const char* ofn, std::vector<PointVectorPair> pVP)
Delaunay exportTriWithCnFun(std::vector<PointVectorPair> pVP, const char* ofn)
{
//    Delaunay Dt = triangulationFromFile(ifn);
//    Delaunay Dt = triangulationSimple();
    
    Delaunay Dt(pVP.begin(), pVP.end());

    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();

    // create PLY output file for outputting the triangulation, with point coordinates, color, normals and triangle facets
    std::ofstream fo;
    fo.open (ofn);
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "comment VCGLIB generated" << std::endl;
    fo << "element vertex " << nv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property float nx" << std::endl;
    fo << "property float ny" << std::endl;
    fo << "property float nz" << std::endl;
    fo << "element face " << nf << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(3);

    // give every vertex from the triangulation an index starting at 0
    // and already print the point coordinates, color and normal of the vertex to the PLY file
    std::map<Vertex_handle, int> Vertices;
    int index = 0;
    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        // assign index to vertex handle. this is needed later for finding the index a certain cell is constructed with
        // TODO: maybe this could already be included in the info vector of every point/vertex, and then instead of using the finite_vertices_iterator,
        // one could iterate over the index vector, so that the point order of the file is determined by the info vector
        // and you would do the find operation to find the point in a std::map<Point, idx>
        Vertices[vft] = index;
        // print data to file
        fo << vft->point() << " "                           // coordinates
           << vft->info() << std::endl;                     // normal
//        fo << vft->point() << " "                           // coordinates
//           << int(std::get<1>(vft->info())[0]) << " "       // red
//           << int(std::get<1>(vft->info())[1]) << " "       // green
//           << int(std::get<1>(vft->info())[2]) << " "       // blue
//           << std::get<0>(vft->info()) << std::endl;        // normal
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

            // print the indicies of each cell to the file
            // vertices is a map of all vertices of the triangulation to an index
            fo << Vertices.find(v)->second << ' ';

        }
        // new cell
        fo << std::endl;
    }
    fo.close();

    std::cout << "Delaunay triangulation done and exported to PLY file!" << std::endl;
    return Dt;

}


/////// ray tracing functions
///
///
int traverseCells(Delaunay Dt, std::map<Cell_handle, std::pair<int, int>> all_cells, //std::map<Vertex_handle, int> all_vertices,
                  Ray ray, Cell_handle current_cell, int oppositeVertex)
{

    if(!Dt.is_infinite(current_cell)){

        //int oppositeVertex = current_cell->index(face_vertex);

        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){


            int idx = (oppositeVertex+i)%4;

            Triangle tri = Dt.triangle(current_cell, idx);
            Facet fac = std::make_pair(current_cell, idx);


            // intersection from here: https://doc.cgal.org/latest/Kernel_23/group__intersection__linear__grp.html
            // get the intersection of the ray and the triangle
            CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
              result = intersection(tri, ray);

            // check if there is an intersection between the current ray and current triangle
            if (result) {

//                std::cout << "in second cell" << std::endl;
//                std::cout << "Face number: " << i << std::endl;


                // first of all I need to locate the current cell in the global context of the triangulation,
                // so I can mark that it is crossed by a ray
                // 1. mark the current cell as "positivelly" traversed, i.e. add one to count
                // add some point this needs to be weighted by the distance from the original point
                (all_cells.find(current_cell)->second.second)++;
                //std::cout << all_cells.find(current_cell)->second.first << ": " << all_cells.find(current_cell)->second.second << std::endl;

                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                Facet mirror_fac = Dt.mirror_facet(fac);
                Triangle mirror_tri = Dt.triangle(mirror_fac);
                // now from this new cell that I am in (get it from mirror_fac), iterate over all the triangles that are not the mirror triangle
                // and check if there is an intersection
                // this should be entered again at if(!Dt.is_infinite(current_cell)), since like this I can check if the cell is not already the infinite cell
                // so start from there to put this into a function

                Cell_handle newCell = mirror_fac.first;
                int newIdx = mirror_fac.second;

                // go on to the next cell, via mirror facet

                // check if ray triangle intersection is a point (probably in most cases)
                // or a line segment (if ray lies inside the triangle)
                // if result is a point

                // for every vertex of the cell, check if there is an intersection
                // because that will generate an infinite loop, because it means that one cell has multiple triangles that are intersected by
                // the ray (since a vertex is shared by two facets)
                // so it will always re-enter the same cell
                // that's why I am - for now - just returning from that cell if it happens
                // what could be done is make this intersection point the new source of the ray
                // simply say ray source = this point, and ray target = new point in the opposite direction of the previous point
                // and than start from RayTracingFun again, because there it will just look for the opposite facet of the intersected vertex
                for(int i=0; i<4; i++){

                    Point pt = newCell->vertex(i)->point();
                    CGAL::cpp11::result_of<Intersect(Point, Ray)>::type
                      point_intersection = intersection(pt, ray);

                    if(point_intersection){
                        const Point* p = boost::get<Point>(&*point_intersection);
                        std::cout << "intersection with a vertex of the cell: " << p << std::endl;
                        return 0;
                    }
                }

                if (const Point* p = boost::get<Point>(&*result))
                {//std::cout << "point of ray-triangle-intersection :  " << *p << std::endl;


                }
                else{
                    const Segment* s = boost::get<Segment>(&*result);
                    std::cout << "segment 3:  " << *s << std::endl;

                    // get the three edges of the current triangle
                    // check how they intersect with the current ray
                    // since the ray passes straight through the triangle in this case

    //                Point end_point = s->target();
    //                int vertexIndexOfEdgeCrossedRay = all_vertices.find(end_point)->second;

    //                std::cout << vertexIndexOfEdgeCrossedRay << std::cout;

                    // for now just return in this case, until it is solved
                    return 0;
                }
                //traverseCells(Dt, all_cells, all_vertices, ray, newCell, newIdx);
                traverseCells(Dt, all_cells, ray, newCell, newIdx);
            }
        }
    }
    return 0;
}




int rayTracingFun(Delaunay Dt){


    // assign index to vertex handle. this is needed later for finding the index a certain cell is constructed with
    // TODO: maybe this could already be included in the info vector of every point/vertex, and then instead of using the finite_vertices_iterator,
    // one could iterate over the index vector, so that the point order of the file is determined by the info vector
    // and you would do the find operation to find the point in a std::map<Point, idx>
    std::map<Vertex_handle, int> all_vertices;
    int vIndex = 0;
    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        all_vertices[vft] = vIndex;
        vIndex++;
    }

    // from Efficient volumetric fusion paper, it follows that I need to sum over all rays to get the score for a tetrahedron
    // every ray means iterating over every vertice
    // "get the score for a tetrahedron" also means I need to keep track of all the tetrahedrons in the Triangulation, thus:
    // give every cell from the triangulation an index starting at 0, and a count, which counts how often it is traversed by rays, also starting at 0
    std::map<Cell_handle, std::pair<int, int>> all_cells;
    int cIndex = 0;
    Delaunay::All_cells_iterator cft;
    for (cft = Dt.all_cells_begin() ; cft != Dt.all_cells_end() ; cft++){
        // assign index to vertex handle. this is needed later for finding the index a certain cell is constructed with
        // TODO: maybe this could already be included in the info vector of every point/vertex, and then instead of using the finite_vertices_iterator,
        // one could iterate over the index vector, so that the point order of the file is determined by the info vector
        // and you would do the find operation to find the point in a std::map<Point, idx>

        // make the map where the key is the cell and the value is a pair of <index, traversal_count>
        // actually, an index is probably not even needed
        all_cells[cft] = std::make_pair(cIndex,0);
        cIndex++;
        //std::cout << std::to_string(cIndex) << std::endl;
    }
    std::cout << "all cells enumerated..." << std::endl;

    // iterate over every vertices = iterate over every ray
    Delaunay::Finite_vertices_iterator vft2;
    int counter = 0;
    for(vft2 = Dt.finite_vertices_begin() ; vft2 != Dt.finite_vertices_end() ; vft2++){

        // ray constructed from point origin to (end of) normal
        Ray ray(vft2->point(), vft2->info());
//        std::cout << ray << std::endl;

        // vector of incident cells to the vertex
        std::vector<Cell_handle> inc_cells;


        // get all incident cells of a vertex: https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html
        Dt.incident_cells(vft2, std::back_inserter(inc_cells));


        // for every cell of incident cells, check if facet(cell, vertex) intersects with vertex normal
        // so this is checking in all directions of a vertex, but we will only have an intersection in (maximum) one direction
        for(std::size_t i=0; i < inc_cells.size(); i++){

            Cell_handle current_cell = inc_cells[i];


            if(!Dt.is_infinite(current_cell))
            {

                int cellBasedVertexIndex = current_cell->index(vft2);
                Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
                Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);



                // intersection from here: https://doc.cgal.org/latest/Kernel_23/group__intersection__linear__grp.html
                // get the intersection of the ray and the triangle
                CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
                  result = intersection(tri, ray);


                // check if there is an intersection between the current ray and current triangle
                if (result){

                    //std::cout << "found first cell in direction of normal / ray" << std::endl;

                    // first of all I need to locate the current cell in the global context of the triangulation,
                    // so I can mark that it is crossed by a ray
                    // 1. mark the current cell as "positivelly" traversed, i.e. add one to count
                    // add some point this needs to be weighted by the distance from the original point
                    (all_cells.find(current_cell)->second.second)++;
                    //std::cout << ++(all_cells.find(current_cell)->second.second) << std::endl;

                    // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                    Facet mirror_fac = Dt.mirror_facet(fac);
                    Triangle mirror_tri = Dt.triangle(mirror_fac);
                    // now from this new cell that I am in (get it from mirror_fac), iterate over all the triangles that are not the mirror triangle
                    // and check if there is an intersection
                    // this should be entered again at if(!Dt.is_infinite(current_cell)), since like this I can check if the cell is not already the infinite cell
                    // so start from there to put this into a function

                    Cell_handle newCell = mirror_fac.first;
                    int newIdx = mirror_fac.second;

                    // if result is a point
                    // check if ray triangle intersection is a point (probably in most cases)
                    // or a line segment (if ray lies inside the triangle)
                    if (const Point* p = boost::get<Point>(&*result)){
                    //std::cout << "point of ray-triangle-intersection :  " << *p << std::endl;
                    }
                    else{
                        const Segment* s = boost::get<Segment>(&*result);
                        //std::cout << "segment 3:  " << *s << std::endl;
                        // for now just return in this case, until it is solved
                        return 0;
                    }

                    //traverseCells(Dt, all_cells, all_vertices, ray, newCell, newIdx);
                    traverseCells(Dt, all_cells, ray, newCell, newIdx);


                }

                else {
                //std::cout << "no intersection" << std::endl;
                // pack this intersection thing into a function that gets called throughout the traversel and returns once there is no more intersection or once you hit the infinite cell

                }

            }

        }
        std::cout << "ray " << std::to_string(++counter) << " done" << std::endl;

    }
    return 0;
}











std::vector<Point> readPlyFun(const char* fname)
{
    // Reads a .ply point set file with normal vectors and colors
    std::vector<Point> points; // store points
    std::ifstream in(fname);
    CGAL::read_ply_points(in, std::back_inserter (points));

    // Display points read
//    for (auto pt = points.begin(); pt != points.end(); pt++)
//    {
//        const Point& p = *pt;
//        std::cout << p << std::endl;
//    }

    std::cout << "PLY file read!" << std::endl;
    return points;
}


int main()
{



//    const char* ifn = "/home/raphael/PhD_local/data/museeZoologic/ALS_TLS_clipped.ply";
//    const char* ofn = "/home/raphael/PhD_local/data/museeZoologic/ALS_TLS_CGAL_meshed.ply";
    const char* ifn = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/2cube_50sampled.ply";
    const char* ofn = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/2cube_CGAL_normals.ply";


//    exportTriWithCnFun(ifn, ofn_test);
    //Delaunay Dt = triangulationSimple();

    auto pts = readPlyFun(ifn);

    auto pVP = estimateNormalsFun(pts);

    Delaunay Dt = exportTriWithCnFun(pVP, ofn);


    //for(int i = 1; i<4; i++){std::cout << i << std::endl;}


    rayTracingFun(Dt);

    return 0;

}






