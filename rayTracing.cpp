#include <cgal_typedefs.h>
#include <fileIO.h>

// TODO: check if I can replace the whole thing with the nearest_vertex(const Point& p, Cell_handle start) function,
// which can be found in the Delaunay_triangulation_3.h file in /usr/lib/CGAL

namespace rayTracing{



////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
std::pair<float, float> cellScore(float dist2, double eig3, bool inside){

    float score_inside;
    float score_outside;
    // noise
    float sigma_d = eig3;
    // scene thickness // good for fontaine dataset is 0.1
    float sigma_o = 1.0;
    // scale of the outside area?? // good for fontaine dataset is 1.0
    float sigma_e = 1.0;
    // not to be confused with the following, which means if I am walking inside/outside
    if(inside){
        //        sigma = 0.05;
        score_inside = (1 - 0.5*exp(-pow((sqrt(dist2)/sigma_d),2)))*exp(-pow((sqrt(dist2)/sigma_o),2));
        score_outside = 0.5*exp(-pow((sqrt(dist2)/sigma_d),2));
    }
    else {
        score_outside = (1 - 0.5*exp(-pow((sqrt(dist2)/sigma_d),2)))*exp(-pow((sqrt(dist2)/sigma_e),2));
        score_inside = 0.5*exp(-pow((sqrt(dist2)/sigma_d),2));
    }
//    if(score_inside != score_inside || score_outside != score_outside)
//        std::cout << score_outside << "  " << score_inside << std::endl;

    // first element is the outside score, second the inside score
    std::pair<float,float> sigmas(score_outside, score_inside);
    return sigmas;
}

//float cellScore(float dist2, double eig3, bool inside){
      // implementation according to Efficient volumetric fusion paper

//    float sigma;
//    if(inside){
//        //        sigma = 0.05;
//        sigma = eig3;
//        // truncate inside ray
//        // in efficient volumetric fusion paper they also introduce outside limit of 3*sigma, simply for having "a shorter walk in the 3DT".
////        if(dist2 > 8*sigma)
////            dist2=0;
//    }
//    else {
////        sigma = 0.25;
//        // good with normals
////        sigma = eig3*5;
//        sigma = eig3;
//    }

//    float S = 1 - exp(-dist2/(2*sigma*sigma));
//    return S;
//}

// TODO: why can the Delaunay be const here? I'm changing the cell scores that are saved inside the Delaunay!
int traverseCells(const Delaunay& Dt, double sigma, Ray ray, Cell_handle current_cell, int oppositeVertex, Point source, bool inside)
{

    // input::
    // &Delaunay            Reference to a Delaunay triangulation
    // ray                  Current ray
    // source               The "Delaunay" point of the current ray
    // current_cell         The cell that has just been entered (in the global context)
    // oppositeVertex       The opposite vertex of the facet where the current_cell was entered
    // inside               Bool that says if the current cell is before or after the point

    if(!Dt.is_infinite(current_cell)){
        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int idx = (oppositeVertex+i)%4;

            Triangle tri = Dt.triangle(current_cell, idx);
//            double tri_area = tri.squared_area();
//            if(tri_area < 1e-100){
//                std::cout << "triangle with " << tri_area << " skipped." << std::endl;
////                return 0;
//            }
//            if(tri.is_degenerate() && !CGAL::is_valid(tri)){
//                std::cout << "triangle is degenerate or not valid " << std::endl;
//                return 0;
//            }
//            if(!do_intersect(tri, ray))
//                return 0;

            typedef CGAL::Cartesian_converter<EPICK,EPECK>                         IK_to_EK;
            typedef CGAL::Cartesian_converter<EPECK,EPICK>                         EK_to_IK;
            IK_to_EK to_exact;
            EK_to_IK to_inexact;

            // btw, here I don't have the problem of ray intersecting multiple cells, because I'm only checking in the current cell
            CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
              result;
            try{
                result = intersection(tri, ray);
            }
            catch(...){
                std::cout << "ray-triangle intersection failed!" << std::endl;
                return 0;
            }


            // check if there is an intersection between the current ray and current triangle
            if(result){
                // check if ray triangle intersection is a point (probably in most cases)
                // or a line segment (if ray lies inside the triangle)
                // if result is a point
                std::pair<float,float> score;
                if (const Point* p = boost::get<Point>(&*result))
                {//std::cout << "point of ray-triangle-intersection :  " << *p << std::endl;
                    // get the distance of this point to the current source:
//                    double dist = sqrt(CGAL::squared_distance(*p, source));
                    float dist2 = CGAL::squared_distance(*p, source);
                    //std::cout << dist2 << std::endl;
                    // dist2 = squared distance from intersection to the point; sigma = noise of the point; inside = bool
                    score = cellScore(dist2, sigma, inside);
                }
                else{
                    const Segment* s = boost::get<Segment>(&*result);
                    std::cout << "segment 3 intersection behind the first cell:  " << *s << std::endl;
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
                current_cell->info().outside_score += score.first;
                current_cell->info().inside_score += score.second;

                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                Facet fac = std::make_pair(current_cell, idx);
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
//                        const Point* p = boost::get<Point>(&*point_intersection);
//                        std::cout << "intersection with a vertex of the cell: " << *p << std::endl;
                        return 0;
                    }
                }
                traverseCells(Dt, sigma, ray, newCell, newIdx, source, inside);
            }
        }
    }
    // put outside score of infinite cell very high
    else{
        current_cell->info().outside_score+=1;
        current_cell->info().inside_score+=0;

    }
    return 0;
}

void firstCell(const Delaunay& Dt, Delaunay::Finite_vertices_iterator& vit, bool inside, bool one_cell){


    // get sigma of the current vertex
    double sigma = vit->info().sigma;


    // ray constructed from point origin to (end of) normal
    // introduces a ray r with source p and with a direction given by v.
    Ray ray(vit->point(), vit->info().sensor_vec);

    // make the inside ray
    // in fact MicMac saves the camera normals pointing away from the camera,
    // so I take the opposite ray for outside traversal and the normal ray for inside
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
//            if(tri_area < 1e-100){
//                std::cout << "triangle with " << tri_area << " skipped." << std::endl;
////                return 0;
//            }
//            if(tri.is_degenerate() && !CGAL::is_valid(tri)){
//                std::cout << "triangle is degenerate or not valid " << std::endl;
//                continue;
//            }
//            if(!do_intersect(tri, ray))
//                continue;

            // intersection from here: https://doc.cgal.org/latest/Kernel_23/group__intersection__linear__grp.html
            // get the intersection of the ray and the triangle
            CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
              result;
            try{
                result = intersection(tri, ray);
            }
            catch(...){
                std::cout << "ray-triangle intersection failed!" << std::endl;
                continue;
            }



            // check if there is an intersection between the current ray and current triangle
            if(result){
                // check if ray triangle intersection is a point (probably in most cases) or a line segment (if ray lies inside the triangle).
                // so far this is not used (and not needed) since I am not handling the unlikely case where the ray goes through a triangle
                // if result is a point
                std::pair<float,float> score;
                Point source = vit->point();
                if(const Point* p = boost::get<Point>(&*result)){

                //std::cout << "point of ray-triangle-intersection :  " << *p << std::endl;

                    // get the distance of this point to the current source:
//                    double dist = sqrt(CGAL::squared_distance(*p, source));
                    float dist2 = CGAL::squared_distance(*p, source);
                    // dist2 = squared distance from intersection to the point; sigma = noise of the point; inside = bool
                    score = cellScore(dist2, sigma, inside);
                }
                // else result is a line
                else{
                    const Segment* s = boost::get<Segment>(&*result);
                    std::cout << "segment 3 intersection in first cell:  " << *s << std::endl;
                    // for now just return in this case, until it is solved
                    //continue;
                    // TOOD: simply calculate the distance to this edge, and then I can also get a score from cellScore()
                    score = std::make_pair(0.0, 0.0);
                }
                // now locate the current cell in the global context of the triangulation,
                // so I can mark that it is crossed by a ray
                current_cell->info().outside_score += score.first;
                current_cell->info().inside_score += score.second;

                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
                // now from this new cell that I am in (get it from mirror_fac), iterate over all the triangles that are not the mirror triangle
                // and check if there is an intersection
                // this should be entered again at if(!Dt.is_infinite(current_cell)), since like this I can check if the cell is not already the infinite cell
                // so start from there to put this into a function
                if(!inside){
                    Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
                    Facet mirror_fac = Dt.mirror_facet(fac);
                    Cell_handle newCell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    // go to next cell
                    traverseCells(Dt, sigma, ray, newCell, newIdx, source, inside);
                }
            // if there was a match, break the loop around this vertex, so it can go to the next one
            // this is only done for speed reasons, it shouldn't have any influence on the label
                // because a ray can only hit more than one facet of a cell if it hits another point of the cell
                // in which case it goes THROUGH a facet of the cell, in which case the intersection is not a point
                // but an edge.
                // it does however make a difference if this is turn on or not. why??
            break;
            }
        }
        // put outside score of infinite cell very high
        else{
//            std::cout << "infinite cell score set" << std::endl;
            current_cell->info().outside_score+=1;
            current_cell->info().inside_score+=0;
        }
    }


}

void rayTracingFun(const Delaunay& Dt, bool one_cell){

    std::cout << "Start tracing rays to every point..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // from Efficient volumetric fusion paper, it follows that I need to sum over all rays to get the score for a tetrahedron
    // iterate over every vertices = iterate over every ray
    // TODO: go in here and map the following to a function
    // this function can than be called with the vft2 iterator
    // and I pass this iterator around and also pass it from within the traversal function
    // when I want to restart from a point that is hit
    Delaunay::Finite_vertices_iterator vit;
    for(vit = Dt.finite_vertices_begin() ; vit != Dt.finite_vertices_end() ; vit++){

        // collect outside votes
        firstCell(Dt, vit, 0, one_cell);    // one_cell currently not used in the correct way
        // collect inside votes
        firstCell(Dt, vit, 1, one_cell);    // one_cell currently not used in the correct way
    }
    // now that all rays have been traced, apply the last function to all the cells:
//    float gamma = 2.0;
//    Cell_map::iterator it;
//    for(it = all_cells.begin(); it!=all_cells.end(); it++)
//    {
//        std::get<1>(it->second) = 1 - exp(-std::get<1>(it->second)/gamma);
//        std::get<2>(it->second) = 1 - exp(-std::get<2>(it->second)/gamma);
//    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Ray tracing done in " << full_duration.count() << "s" << std::endl;
}

// end of namespace rayTracing
}






