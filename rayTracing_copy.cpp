#include <cgal_typedefs.h>
#include <fileIO.h>

// TODO: check if I can replace the whole thing with the nearest_vertex(const Point& p, Cell_handle start) function,
// which can be found in the Delaunay_triangulation_3.h file in /usr/lib/CGAL

namespace rayTracingCopy{

std::vector<double> crossProduct(std::vector<double> a, std::vector<double> b){
    std::vector<double> cp(3);
    cp[0] = a[1] * b[2] - a[2] * b[1];
    cp[1] = a[0] * b[2] - a[2] * b[0];
    cp[2] = a[0] * b[1] - a[1] * b[0];

    return cp;
}
double dotProduct(std::vector<double> a, std::vector<double> b){
    return  a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


bool rayTriangleIntersection(Point& rayOrigin,
                           Vector& rayVector,
                           Triangle& inTriangle,
                           Point& outIntersectionPoint){
    const double EPSILON = 0.0000001;

    std::vector<double> rayO{rayOrigin.x(), rayOrigin.y(), rayOrigin.z()};
    std::vector<double> rayV{rayVector.x(), rayVector.y(), rayVector.z()};

    std::vector<double> vertex0{inTriangle.vertex(0).x(), inTriangle.vertex(0).y(), inTriangle.vertex(0).z()};
    std::vector<double> vertex1{inTriangle.vertex(1).x(), inTriangle.vertex(1).y(), inTriangle.vertex(1).z()};
    std::vector<double> vertex2{inTriangle.vertex(2).x(), inTriangle.vertex(2).y(), inTriangle.vertex(2).z()};

    double a,f,u,v;

    std::vector<double> edge1{  vertex1[0] - vertex0[0],
                                vertex1[1] - vertex0[1],
                                vertex1[2] - vertex0[2]};

    std::vector<double> edge2{  vertex2[0] - vertex0[0],
                                vertex2[1] - vertex0[1],
                                vertex2[2] - vertex0[2]};

    std::vector<double> h = crossProduct(rayV,edge2);
    a = dotProduct(edge1,h);
    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.
    f = 1.0/a;
    std::vector<double> s{  rayO[0] - vertex0[0],
                            rayO[1] - vertex0[1],
                            rayO[2] - vertex0[2]};
    u = f * dotProduct(s,h);
    if (u < 0.0 || u > 1.0)
        return false;
    std::vector<double> q = crossProduct(s, edge1);
    v = f * dotProduct(rayV, q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * dotProduct(edge2, q);
    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = rayOrigin + rayVector * t;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}

////////////////////////////////////////////////////////////
/////////////////// ray tracing functions //////////////////
////////////////////////////////////////////////////////////
// TODO: why can the Delaunay be const here? I'm changing the cell scores that are saved inside the Delaunay!
int traverseCells(Delaunay& Dt,
                  double sigma, Ray& ray,
                  Cell_handle& current_cell, int oppositeVertex, Point& source,
                  bool inside, std::unordered_set<Cell_handle>& processed_cells)
{

    // input::
    // &Delaunay            Reference to a Delaunay triangulation
    // ray                  Current ray
    // source               The "Delaunay" point of the current ray
    // current_cell         The cell that has just been entered (in the global context)
    // oppositeVertex       The opposite vertex of the facet where the current_cell was entered
    // inside               Bool that says if the current cell is before or after the point

    if(processed_cells.find(current_cell) != processed_cells.end())
        return 0;

    if(!Dt.is_infinite(current_cell)){
        // iterate over the faces of the current cell
        for(int i=1; i<4; i++){
            // I'm starting here at the oppositeVertex+1 of the facet, so it will not go back to the same cell it came from
            int idx = (oppositeVertex+i)%4;

            Facet fac = std::make_pair(current_cell, idx);
            if(Dt.is_infinite(fac))
                return 0;

//            Triangle tri = Dt.triangle(current_cell, idx);
            Triangle tri = Dt.triangle(fac);
            if(tri.is_degenerate()){
                std::cout<<"tri is degenerate" << std::endl;
                return 0;
            }
            Point intersectionPoint;
            Vector rayV = ray.to_vector();

            std::cout << tri[0][0] << std::endl;
//            std::cout << tri[1][1] << std::endl;
//            std::cout << tri[2][1] << std::endl;
            bool result = rayTriangleIntersection(source, rayV, tri, intersectionPoint);
            std::pair<float,float> score;
            // check if there is an intersection between the current ray and current triangle
            if(result){

                processed_cells.insert(current_cell);
                Triangle tri = Dt.triangle(fac);
                if(tri.is_degenerate()){
                    std::cout<<"tri is degenerate" << std::endl;
                    return 0;
                }
                std::cout << tri[0][0] << std::endl;
//                std::cout << tri[1][1] << std::endl;
//                std::cout << tri[2][1] << std::endl;

                Facet mirror_fac = Dt.mirror_facet(fac);
                if(Dt.is_infinite(mirror_fac)){
                    std::cout << "infinite facet hit" << std::endl;
                    return 0;
                }


                Triangle mirror_tri = Dt.triangle(mirror_fac);
                if(mirror_tri.is_degenerate()){
                    std::cout<<"tri is degenerate" << std::endl;
                    return 0;
                }
                std::cout << mirror_tri[0][0] << std::endl;
//                std::cout << mirror_tri[1][1] << std::endl;
//                std::cout << mirror_tri[2][1] << std::endl;

                Cell_handle newCell = mirror_fac.first;
                int newIdx = mirror_fac.second;
                for(int i=0; i<4; i++){
                    Point pt = newCell->vertex(i)->point();
                    CGAL::cpp11::result_of<Intersect(Point, Ray)>::type
                      point_intersection = intersection(pt, ray);
                    if(point_intersection){
    //                        const Point* p = boost::get<Point>(&*point_intersection);
                        std::cout << "intersection with a vertex of the cell: " << *point_intersection << std::endl;
                        return 0;
                    }
                }
                traverseCells(Dt, sigma, ray, newCell, newIdx, source, inside, processed_cells);
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

void firstCell(Delaunay& Dt, Delaunay::Finite_vertices_iterator& vit, bool inside){

    std::unordered_set<Cell_handle> processed_cells;
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
//            Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
            Triangle tri = Dt.triangle(current_cell, cellBasedVertexIndex);
            std::cout << tri[0][0] << std::endl;
//            std::cout << tri[1] << std::endl;
//            std::cout << tri[2] << std::endl;

            Point intersectionPoint;
            Point source = vit->point();
            Vector rayV = ray.to_vector();
            bool result = rayTriangleIntersection(source, rayV, tri, intersectionPoint);
            std::pair<float,float> score;
            // check if there is an intersection between the current ray and current triangle
            if(result){
                float dist2 = CGAL::squared_distance(intersectionPoint, source);

                current_cell->info().outside_score += score.first;
                current_cell->info().inside_score += score.second;
                processed_cells.insert(current_cell);

                if(!inside){
                    Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
                    Facet mirror_fac = Dt.mirror_facet(fac);
                    Cell_handle newCell = mirror_fac.first;
                    int newIdx = mirror_fac.second;
                    // go to next cell
                    traverseCells(Dt, sigma, ray, newCell, newIdx, source, inside, processed_cells);
                }
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

void iterateOverAllFacets(Delaunay& Dt, Delaunay::Finite_vertices_iterator vit){

    // vector of incident cells to the current vertex (from vertex iterator vit)
    std::vector<Cell_handle> inc_cells;
    Dt.incident_cells(vit, std::back_inserter(inc_cells));
    for(std::size_t i=0; i < inc_cells.size(); i++){
        Cell_handle current_cell = inc_cells[i];
        for(int cellBasedVertexIndex=0; cellBasedVertexIndex < 3; cellBasedVertexIndex++){
            Facet fac = std::make_pair(current_cell, cellBasedVertexIndex);
            if(Dt.is_infinite(fac)){
                continue;
            }
            Triangle tri = Dt.triangle(fac);
            std::cout << tri[0] << std::endl;
            std::cout << tri[1] << std::endl;
            std::cout << tri[2] << std::endl;
            Facet mirror_fac = Dt.mirror_facet(fac);
            Triangle mirror_tri = Dt.triangle(mirror_fac);
            std::cout << mirror_tri[0] << std::endl;
            std::cout << mirror_tri[1] << std::endl;
            std::cout << mirror_tri[2] << std::endl;
        }
    }
}

void iterateOverFacets(Delaunay& Dt){
    Delaunay::Finite_facets_iterator vit;
    int count = 0;
    for(vit = Dt.finite_facets_begin() ; vit != Dt.finite_facets_end() ; vit++){
        Triangle tri = Dt.triangle(*vit);
        std::cout << count++ << std::endl;
        std::cout << tri[0][0] << std::endl;
//        std::cout << tri[1][0] << std::endl;
//        std::cout << tri[2][0] << std::endl;
        Facet mfac = Dt.mirror_facet(*vit);
        if(Dt.is_infinite(mfac)){
            continue;
        }
        Triangle mtri = Dt.triangle(mfac);
        std::cout << mtri[0][0] << std::endl;
//        std::cout << mtri[1][0] << std::endl;
//        std::cout << mtri[2][0] << std::endl;

    }
}


void rayTracingFun(Delaunay& Dt){

    std::cout << "Start tracing rays to every point..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

//    iterateOverFacets(Dt);


    int count=0;
    Delaunay::Finite_vertices_iterator vit;
    for(vit = Dt.finite_vertices_begin() ; vit != Dt.finite_vertices_end() ; vit++){

        std::cout << count++ << std::endl;
        // collect outside votes
        firstCell(Dt, vit, 0);    // one_cell currently not used in the correct way
        // collect inside votes
        firstCell(Dt, vit, 1);    // one_cell currently not used in the correct way
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Ray tracing done in " << full_duration.count() << "s" << std::endl;
}

// end of namespace rayTracing
}






