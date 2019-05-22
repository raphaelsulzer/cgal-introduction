#include "fileIO.cpp"

typedef Kernel::Ray_3                                               Ray;
typedef Kernel::Triangle_3                                          Triangle;
typedef Kernel::Intersect_3                                         Intersect;
typedef Kernel::Segment_3                                           Segment;

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
float cellScore(float dist2, bool inside){

    float sigma;
    if(inside){
        sigma = 0.2;
        // truncate inside ray
        if(dist2 > 3*sigma)
            dist2=0;
    }
    else {
        sigma = 0.3;
    }

    float S = 1 - exp(-dist2/(2*sigma*sigma));
    return S;
}

int traverseCells(const Delaunay& Dt, Cell_map& all_cells, Ray ray, Cell_handle current_cell, int oppositeVertex, Point source, bool inside)
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

    Vector camera;
    if(vit->info() == 0){
        camera = Vector(-0.360035117847216257, -0.0243440832633011854, 0.284447917373549908);
    }
    else if(vit->info() == 1){
        camera = Vector(-0.144933279455665032, -10.4598329635251179, -6.34409732148353278);
    }
    else{
        camera = Vector(-0.229706673957515983, 9.05508818222588552, -9.21427702085086331);
    }

    // ray constructed from point origin to (end of) normal
//    Ray ray(vit->point(), vit->info());
    Ray ray(vit->point(), camera);

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
    float gamma = 1.0;
    Cell_map::iterator it;
    for(it = all_cells.begin(); it!=all_cells.end(); it++)
    {
        std::get<1>(it->second) = 1 - exp(-std::get<1>(it->second)/gamma);
        std::get<2>(it->second) = 1 - exp(-std::get<2>(it->second)/gamma);
    }
}


