#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/intersections.h>

#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_ply_points.h>

#include <CGAL/Euclidean_distance.h>

// for GCoptimization
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "GCoptimization.h"


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
typedef Kernel::Ray_3                                               Ray;
typedef Kernel::Triangle_3                                          Triangle;
typedef Kernel::Intersect_3                                         Intersect;
typedef Kernel::Segment_3                                           Segment;

// map cell of the Dt, to an index, the outside score, the inside score, and the final label
typedef std::map<Cell_handle, std::tuple<int, float, float, int>>            Cell_map;
//typedef std::map<Cell_handle, std::pair<double, double>>            Cell_map;


typedef std::pair<Point, Vector> PointVectorPair;

typedef CGAL::cpp11::tuple<Point, Vector> PN;
typedef CGAL::Nth_of_tuple_property_map<0, PN> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PN> Normal_map;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

////////////////////////////////////////////////////////////
/////////////////// preprocessing functions ////////////////
////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
float cellScore(float dist2, bool inside){

    float sigma;
    if(inside){
        sigma = 0.05;
        // truncate inside ray
        if(dist2 > 3*sigma)
            dist2=0;
    }
    else {
        sigma = 0.1;
    }

    float S = 1 - exp(-dist2/(2*sigma*sigma));
    return S;
}

int traverseCells(const Delaunay& Dt, Cell_map& all_cells, Ray ray, Cell_handle current_cell, int oppositeVertex, Point source, bool inside)
{
    if(!Dt.is_infinite(current_cell)){
        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int idx = (oppositeVertex+i)%4;

            Triangle tri = Dt.triangle(current_cell, idx);
            Facet fac = std::make_pair(current_cell, idx);

            // btw, here I don't have the problem of ray intersecting multiple cells, because I'm only checking in the current cell
            CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
              result = intersection(tri, ray);

            // check if there is an intersection between the current ray and current triangle
            if (result) {
                // check if ray triangle intersection is a point (probably in most cases)
                // or a line segment (if ray lies inside the triangle)
                // if result is a point
                float score;
                if (const Point* p = boost::get<Point>(&*result))
                {//std::cout << "point of ray-triangle-intersection :  " << *p << std::endl;
                    // get the distance of this point to the current source:
//                    double dist = sqrt(CGAL::squared_distance(*p, source));
                    float dist2 = CGAL::squared_distance(*p, source);
                    score = cellScore(dist2, inside);

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
                // now locate the current cell in the global context of the triangulation,
                // so I can mark that it is crossed by a ray
                if(!inside){
                    std::get<1>(all_cells.find(current_cell)->second) += score;
                }
                else {
                    std::get<2>(all_cells.find(current_cell)->second) += score;
                }

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
                        const Point* p = boost::get<Point>(&*point_intersection);
                        std::cout << "intersection with a vertex of the cell: " << *p << std::endl;
                        return 0;
                    }
                }
                traverseCells(Dt, all_cells, ray, newCell, newIdx, source, inside);
            }
        }
    }
    // put outside score of infinite cell very high
    else{
        // set outside score
        std::get<1>(all_cells.find(current_cell)->second)+=1;
        // set inside score
        std::get<2>(all_cells.find(current_cell)->second)+=0;
    }
    return 0;
}

void firstCell(const Delaunay& Dt, Delaunay::Finite_vertices_iterator& vit, Cell_map& all_cells, bool inside){

    // ray constructed from point origin to (end of) normal
    Ray ray(vit->point(), vit->info());

    // make the inside ray
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
            if (result){
                // check if ray triangle intersection is a point (probably in most cases) or a line segment (if ray lies inside the triangle).
                // so far this is not used (and not needed) since I am not handling the unlikely case where the ray goes through a triangle
                // if result is a point
                float score;
                Point source = vit->point();
                if (const Point* p = boost::get<Point>(&*result)){
                //std::cout << "point of ray-triangle-intersection :  " << *p << std::endl;

                    // get the distance of this point to the current source:
//                    double dist = sqrt(CGAL::squared_distance(*p, source));
                    float dist2 = CGAL::squared_distance(*p, source);
                    score = cellScore(dist2, inside);
                }
                // else result is a line
                else{
                    const Segment* s = boost::get<Segment>(&*result);
                    std::cout << "segment 3:  " << *s << std::endl;
                    // for now just return in this case, until it is solved
                    //continue;
                    score = 0;
                }
                // now locate the current cell in the global context of the triangulation,
                // so I can mark that it is crossed by a ray
                if(!inside){
                    std::get<1>(all_cells.find(current_cell)->second) += score;
                }
                else {
                    std::get<2>(all_cells.find(current_cell)->second) += score;
                }
                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                // now from this new cell that I am in (get it from mirror_fac), iterate over all the triangles that are not the mirror triangle
                // and check if there is an intersection
                // this should be entered again at if(!Dt.is_infinite(current_cell)), since like this I can check if the cell is not already the infinite cell
                // so start from there to put this into a function
                Facet mirror_fac = Dt.mirror_facet(fac);
                Cell_handle newCell = mirror_fac.first;
                int newIdx = mirror_fac.second;

                // go to next cell
                traverseCells(Dt, all_cells, ray, newCell, newIdx, source, inside);
            }
        }
        // put outside score of infinite cell very high
        else{
//            std::cout << "infinite cell score set" << std::endl;
            // set outside score
            std::get<1>(all_cells.find(current_cell)->second)+=1.0;
            // set inside score
            std::get<2>(all_cells.find(current_cell)->second)+=0.0;
        }
    }


}

void rayTracingFun(const Delaunay& Dt, Cell_map& all_cells){

    std::cout << "Start tracing rays to every point..." << std::endl;

    // "get the score for a tetrahedron" also means I need to keep track of all the tetrahedrons in the Triangulation, thus:
    // give every cell from the triangulation a count, which counts how often it is traversed by rays starting at 0
    int cindex = 0;
    Delaunay::All_cells_iterator cft;
    for (cft = Dt.all_cells_begin(); cft != Dt.all_cells_end(); cft++){
        // make the map where the key is the cell and the value is:
        // index
        std::get<0>(all_cells[cft]) = cindex;
        // outside votes
        std::get<1>(all_cells[cft]) = 0.0;
        // inside votes
        std::get<2>(all_cells[cft]) = 0.0;
        // final label
        std::get<3>(all_cells[cft]) = 0;
        cindex++;
    }
    // from Efficient volumetric fusion paper, it follows that I need to sum over all rays to get the score for a tetrahedron
    // iterate over every vertices = iterate over every ray
    // TODO: go in here and map the following to a function
    // this function can than be called with the vft2 iterator
    // and I pass this iterator around and also pass it from within the traversal function
    // when I want to restart from a point that is hit
    Delaunay::Finite_vertices_iterator vit;
    for(vit = Dt.finite_vertices_begin() ; vit != Dt.finite_vertices_end() ; vit++){
        // collect outside votes
        firstCell(Dt, vit, all_cells, 0);
        // collect inside votes
        firstCell(Dt, vit, all_cells, 1);
    }
    // now that all rays have been traced, apply the last function to all the cells:
    float gamma = 3.0;
    Cell_map::iterator it;
    for(it = all_cells.begin(); it!=all_cells.end(); it++)
    {
        std::get<1>(it->second) = 1 - exp(-std::get<1>(it->second)/gamma);
        std::get<2>(it->second) = 1 - exp(-std::get<2>(it->second)/gamma);
    }
}

//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
/// //////////////////////////////////////////////////////
std::vector<Point> readPlyFun(const char* fname)
{
    // Reads a .ply point set file
    std::vector<Point> points; // store points
    std::ifstream in(fname);
    CGAL::read_ply_points(in, std::back_inserter (points));

    std::cout << "PLY file read!" << std::endl;
    return points;
}

std::vector<PN> readPlyWithCnFun(const char* fname)
{
    std::vector<PN> points; // store points
    std::ifstream in(fname);

    CGAL::read_ply_points_with_properties
      (in,
       std::back_inserter (points),
       CGAL::make_ply_point_reader (Point_map()),
       CGAL::make_ply_normal_reader (Normal_map()));

    std::cout << "PLY file read!" << std::endl;
    return points;
}

// generate a Delaunay triangulation from a PLY file
Delaunay triangulationFromFile(const char* ifn)
{
   // get data as vector of tuples(point, normal, color, intensity)
   auto ply = readPlyWithCnFun(ifn);

   std::vector<Point> points;
   std::vector<Vector> infos;
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
void exportSoup(const Delaunay& Dt, Cell_map all_cells, const char* ofn)
{
    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();

    // calculate how many faces to print
    Delaunay::Finite_facets_iterator fft;
    int deletedFaceCount = 0;
    for(fft = Dt.finite_facets_begin() ; fft != Dt.finite_facets_end() ; fft++){
        Cell_handle c = fft->first;
        int clabel = std::get<3>(all_cells.find(c)->second);
        Cell_handle m = Dt.mirror_facet(*fft).first;
        int mlabel = std::get<3>(all_cells.find(m)->second);
        if(clabel == mlabel){deletedFaceCount++;}
    }
    int sub = nf - deletedFaceCount;

    // create PLY output file for outputting the triangulation, with point coordinates, color, normals and triangle facets
    std::fstream fo;
    fo.open(ofn, std::fstream::out);

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
    fo << "element face " << sub << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(3);

    // give every vertex from the triangulation an index starting at 0
    // and already print the point coordinates, color and normal of the vertex to the PLY file
    std::map<Vertex_handle, int> Vertices;
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
        // with only outside labelling:
        // cell of current facet
        float cscore = std::get<1>(all_cells.find(c)->second);
        // cell of mirror facet and see if the labels are different
        Cell_handle c1 = Dt.mirror_facet(*fft).first;
        float mscore = std::get<1>(all_cells.find(c1)->second);

//        if((cscore !=0.0 && mscore !=0.0) || (cscore != 0.0 && Dt.is_infinite(c1)) || (mscore !=0.0 && Dt.is_infinite(c)) ){
//            deletedFaceCount++;
//            continue;
//        }

        // with inside labelling:
//        double coutside = all_cells.find(c)->second.first;
//        double cinside = all_cells.find(c)->second.second;
//        double moutside = all_cells.find(c1)->second.first;
//        double minside = all_cells.find(c1)->second.second;

//        double fscore = sqrt(pow(coutside-moutside,2)+pow(cinside-minside,2));
//        std::cout << fscore << std::endl;

//        if(!(fscore > 0.1)){
//            deletedFaceCount++;
//            continue;
//        }

//        // with GC labelling:
        int clabel = std::get<3>(all_cells.find(c)->second);

        Cell_handle m = Dt.mirror_facet(*fft).first;
        int mlabel = std::get<3>(all_cells.find(m)->second);

        if(clabel == mlabel){
//            deletedFaceCount++;
            continue;}

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

    fo.close();

    std::cout << "before face count: " << nf << std::endl;
    std::cout << "remaining faces: " << sub << std::endl;
    std::cout << "Delaunay triangulation done and exported to PLY file!" << std::endl;
}

//// in this version, set data and smoothness terms using arrays
//// grid neighborhood is set up "manually". Uses spatially varying terms. Namely
//// V(p1,p2,l1,l2) = w_{p1,p2}*[min((l1-l2)*(l1-l2),4)], with
//// w_{p1,p2} = p1+p2 if |p1-p2| == 1 and w_{p1,p2} = p1*p2 if |p1-p2| is not 1
//std::pair<std::map<Cell_handle, int>, std::vector<int>> GeneralGraph_DArraySArraySpatVarying(std::pair<Delaunay&, Cell_map&> dt_cells, std::map<Cell_handle, int>& cell_indexMap, std::vector<int> result, int num_iterations)
void GeneralGraph_DArraySArraySpatVarying(const Delaunay& Dt, Cell_map& all_cells, int num_iterations)
{

    std::cout << "Starting Optimization..." << std::endl;

    int num_cells = all_cells.size();
    int num_labels = 2;
    GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_cells,num_labels);


    // first set up the array for data costs and set the initial label in the same loop
    float *data = new float[num_cells*num_labels];
    int idx;
    Cell_map::iterator it;
    // iterate over the all_cells map
    for(it = all_cells.begin(); it!=all_cells.end(); it++)
    {
        // use the idx of the all_cells map
        idx = std::get<0>(it->second);
        // I am initializing my label s.t. if the outside vote is bigger than the inside vote, than it should be labelled 1 (outside) - and vice versa (0 = inside)
        // that means the cost/penalty for labelling an cell that initially has label 1 with the opposite label, is the opposite vote (so the inside vote)
        // so labelling 0 costs outside votes
        data[idx*2+0] = std::get<1>(it->second);
        // and labelling 1 costs inside votes
        data[idx*2+1] = std::get<2>(it->second);

        // set an initial label
        if(std::get<1>(it->second) > std::get<2>(it->second))
            {gc->setLabel(idx, 1);}
        else
            {gc->setLabel(idx, 0);}
    }

    // next set up the array for smooth costs
    float *smooth = new float[num_labels*num_labels];
    for ( int l1 = 0; l1 < num_labels; l1++ )
        for (int l2 = 0; l2 < num_labels; l2++ )
            if(l1 == l2){smooth[l1+l2*num_labels] = 0.0;}
            else{smooth[l1+l2*num_labels] = 1.0;}

    try{
        gc->setDataCost(data);
        gc->setSmoothCost(smooth);

        // set neighborhood:
        int current_index;
        Cell_map::iterator it2;
        // iterate over the all_cells map
        for(it2 = all_cells.begin(); it2!=all_cells.end(); it2++)
        {

            current_index = std::get<0>(it2->second);
            Cell_handle current_cell = it2->first;
            for(int i = 0; i < 4; i++){

                Cell_handle neighbour_cell = current_cell->neighbor(i);
                int neighbour_index = std::get<0>(all_cells.find(neighbour_cell)->second);

                // if current_cell and neighbour_cell are BOTH infinite, then continue
                if(Dt.is_infinite(current_cell) && Dt.is_infinite(neighbour_cell)){
                    continue;
                }

                // prevent to call setNeighbour(s2,s1) if setNeighbour(s1,s2)was already called
                if(neighbour_index < current_index)
                    continue;

                // since i is giving me the cell that is opposite of vertex i, as well as the facet that is opposite of vertex i, I can just use that same index
                Triangle tri = Dt.triangle(current_cell, i);
                float area = sqrt(tri.squared_area());

                // call the neighbourhood function
                float area_weight = 1.5;    // this clearly shows that minimization does something, since energy changes when weight is changed
                gc->setNeighbors(current_index, neighbour_index, area_weight*area);
            }
        }

        // Optimization
        std::cout << "Before optimization data energy is " << gc->giveDataEnergy() << std::endl;
//        for(int i = 0; i < num_cells; i++){
//            std::cout << "cell " << i << " label 0: " << data[i*2+0] << std::endl;
//            std::cout << "cell " << i << " label 1: " << data[i*2+1] << std::endl;
//        }
        std::cout << "Before optimization smoothness energy is " << gc->giveSmoothEnergy() << std::endl;

        std::cout << "Before optimization energy is " << gc->compute_energy() << std::endl;
        // use swap because it is good when you have two labels
        gc->swap(num_iterations);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
//        gc->expansion(num_iterations);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
        std::cout << "After optimization energy is " << gc->compute_energy() << std::endl;

        // save the results in the all_cells Cell_map
        int idx3;
        Cell_map::iterator it3;
        // iterate over the all_cells map
        for(it3 = all_cells.begin(); it3!=all_cells.end(); it3++)
        {
            idx3 = std::get<0>(it3->second);
            std::get<3>(it3->second) = gc->whatLabel(idx3);
        }

        delete gc;
    }
    catch (GCException e){
        e.Report();
    }
}

void energyMin(const Delaunay& Dt, Cell_map& all_cells, int num_iterations)
{

    int num_cells = all_cells.size();
    int num_labels = 2;

    int* label = new int[num_cells];

    // first set up the array for data costs and set the initial label in the same loop

    int current_idx;
    Cell_map::iterator it;
    // iterate over the all_cells map
    for(it = all_cells.begin(); it!=all_cells.end(); it++)
    {
        Cell_handle current_cell = it->first;

        // data term
        current_idx = std::get<0>(it->second);
//        data = std::get<1>(it->second) + std::get<2>(it->second);

        // smoothness term
        if(std::get<1>(it->second) > std::get<2>(it->second))
            {label[current_idx] = 1;}
        else
            {label[current_idx] = 0;}
    }


    float data_energy = 0.0;
    float smoothness_energy = 0.0;
    float total_energy = 0.0;
    // iterate over the all_cells map
    for(it = all_cells.begin(); it!=all_cells.end(); it++)
    {
        Cell_handle current_cell = it->first;
        current_idx = std::get<0>(it->second);

        // data term
        // so I am initializing my label s.t. if the outside vote is bigger than the inside vote, than it should be labelled 1 - and vice versa
        // that means the cost for labelling an cell that has label 1 with label 0, is the opposite vote (so the inside vote)
        float data;
        if(label[current_idx] == 0)
            data = std::get<1>(it->second);
        else
            data = std::get<2>(it->second);
        data_energy+=data;

        for(int i = 0; i < 4; i++){

            Cell_handle neighbour_cell = current_cell->neighbor(i);
            int neighbour_idx = std::get<0>(all_cells.find(neighbour_cell)->second);

            // if current_cell and neighbour_cell are BOTH infinite, then continue
            if(Dt.is_infinite(current_cell) && Dt.is_infinite(neighbour_cell)){
                continue;
            }

            // prevent to call setNeighbour(s2,s1) if setNeighbour(s1,s2)was already called
            if(neighbour_idx < current_idx)
                continue;

            // since i is giving me the cell that is opposite of vertex i, as well as the facet that is opposite of vertex i, I can just use that same index
            Triangle tri = Dt.triangle(current_cell, i);
            float area = sqrt(tri.squared_area());

            // call the neighbourhood function
            float area_weight = 1.5;    // this clearly shows that minimization does something, since energy changes when weight is changed

            float smooth;
            if(label[current_idx]!=label[neighbour_idx])
                smooth = area_weight*area;
            else
                smooth = 0.0;
            smoothness_energy+=smooth;
        }
    }

    total_energy = data_energy + smoothness_energy;

    std::cout << "Calculated total=data+smoothness energy is: " << total_energy << "=" << data_energy << "+" << smoothness_energy << std::endl;


}




//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int main()
{
//    const char* ifn = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/3cube_10000sampled_messyNormals.ply";
//    const char* ofn = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/3cube_CGAL_pruned.ply";
    const char* ifn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/3cube_10000sampled_messyNormals.ply";
    const char* ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/3cube_CGAL_pruned.ply";


    Delaunay Dt = triangulationFromFile(ifn);

    Cell_map all_cells;

    rayTracingFun(Dt, all_cells);

    energyMin(Dt, all_cells, 2);

    GeneralGraph_DArraySArraySpatVarying(Dt, all_cells, -1);

    exportSoup(Dt, all_cells, ofn);

    return 0;
}






