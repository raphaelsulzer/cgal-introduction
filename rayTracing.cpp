#include <cgal_typedefs.h>
#include <fileIO.h>

// TODO: check if I can replace the whole thing with the nearest_vertex(const Point& p, Cell_handle start) function,
// which can be found in the Delaunay_triangulation_3.h file in /usr/lib/CGAL

namespace rayTracing{


//// cross product of two vector array.
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

//bool rayTriangleIntersection(Point& rayOrigin,
//                           Vector& rayVector,
//                           Triangle& inTriangle,
//                           Point& outIntersectionPoint){

//    const float EPSILON = 0.0000001;
//    Point vertex0 = inTriangle.vertex(0);
//    Point vertex1 = inTriangle.vertex(1);
//    Point vertex2 = inTriangle.vertex(2);
//    Vector edge1, edge2, h, s, q;
//    float a,f,u,v;
//    edge1 = Vector(vertex1.x() - vertex0.x(),
//                   vertex1.y() - vertex0.y(),
//                   vertex1.z() - vertex0.z());
//    edge2 = Vector(vertex2.x() - vertex0.x(),
//                   vertex2.y() - vertex0.y(),
//                   vertex2.z() - vertex0.z());
//    h = crossProduct(rayVector,edge2);
//    a = dotProduct(edge1,h);
//    if (a > -EPSILON && a < EPSILON)
//        return false;    // This ray is parallel to this triangle.
//    f = 1.0/a;
//    s = rayOrigin - vertex0;
//    u = f * dotProduct(s,h);
//    if (u < 0.0 || u > 1.0)
//        return false;
//    q = crossProduct(s,edge1);
//    v = f * dotProduct(rayVector,q);
//    if (v < 0.0 || u + v > 1.0)
//        return false;
//    // At this stage we can compute t to find out where the intersection point is on the line.
//    float t = f * dotProduct(edge2,q);
//    if (t > EPSILON) // ray intersection
//    {
//        outIntersectionPoint = rayOrigin + rayVector * t;
//        return true;
//    }
//    else // This means that there is a line intersection but not a ray intersection.
//        return false;
//}

bool rayTriangleIntersection(Point& rayOrigin,
                           Vector& rayVector,
                           Triangle& inTriangle,
                           Point& outIntersectionPoint){
    using namespace Eigen;
    const double EPSILON = 0.0000001;

    // init "Eigen" vectors
//    Vector3d rayO, rayV;
    Vector3d rayO(rayOrigin.x(), rayOrigin.y(), rayOrigin.z());
    Vector3d rayV(rayVector.x(), rayVector.y(), rayVector.z());

//    Vector3d vertex0, vertex1, vertex2;
//    vertex0 << inTriangle.vertex(0).x(), inTriangle.vertex(0).y(), inTriangle.vertex(0).z();
//    vertex1 << inTriangle.vertex(1).x(), inTriangle.vertex(1).y(), inTriangle.vertex(1).z();
//    vertex2 << inTriangle.vertex(2).x(), inTriangle.vertex(2).y(), inTriangle.vertex(2).z();
    Vector3d vertex0(inTriangle.vertex(0).x(), inTriangle.vertex(0).y(), inTriangle.vertex(0).z());
    Vector3d vertex1(inTriangle.vertex(1).x(), inTriangle.vertex(1).y(), inTriangle.vertex(1).z());
    Vector3d vertex2(inTriangle.vertex(2).x(), inTriangle.vertex(2).y(), inTriangle.vertex(2).z());


    Vector3d edge1, edge2, h, s, q;
    double a,f,u,v;
//    edge1 = vertex1 - vertex0;
//    edge2 = vertex2 - vertex0;
    edge1(0) = vertex1(0) - vertex0(0);
    edge1(1) = vertex1(1) - vertex0(1);
    edge1(2) = vertex1(2) - vertex0(2);
    edge2(0) = vertex2(0) - vertex0(0);
    edge2(1) = vertex2(1) - vertex0(1);
    edge2(2) = vertex2(2) - vertex0(2);
    h = rayV.cross(edge2);
    a = edge1.dot(h);
    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.
    f = 1.0/a;
    s = rayO - vertex0;
    u = f * s.dot(h);
    if (u < 0.0 || u > 1.0)
        return false;
    q = s.cross(edge1);
    v = f * rayV.dot(q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * edge2.dot(q);
    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = rayOrigin + rayVector * t;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}

//bool rayTriangleIntersection(Point rayOrigin,
//                           Vector rayVector,
//                           Triangle inTriangle,
//                           Point& outIntersectionPoint){
//    const double EPSILON = 0.0001;

//    std::vector<double> rayO{rayOrigin.x(), rayOrigin.y(), rayOrigin.z()};
//    std::vector<double> rayV{rayVector.x(), rayVector.y(), rayVector.z()};

//    std::vector<double> vertex0{inTriangle.vertex(0).x(), inTriangle.vertex(0).y(), inTriangle.vertex(0).z()};
//    std::vector<double> vertex1{inTriangle.vertex(1).x(), inTriangle.vertex(1).y(), inTriangle.vertex(1).z()};
//    std::vector<double> vertex2{inTriangle.vertex(2).x(), inTriangle.vertex(2).y(), inTriangle.vertex(2).z()};

//    double a,f,u,v;

//    std::vector<double> edge1{  vertex1[0] - vertex0[0],
//                                vertex1[1] - vertex0[1],
//                                vertex1[2] - vertex0[2]};

//    std::vector<double> edge2{  vertex2[0] - vertex0[0],
//                                vertex2[1] - vertex0[1],
//                                vertex2[2] - vertex0[2]};

//    std::vector<double> h = crossProduct(rayV,edge2);
//    a = dotProduct(edge1,h);
//    if (a > -EPSILON && a < EPSILON)
//        return false;    // This ray is parallel to this triangle.
//    f = 1.0/a;
//    std::vector<double> s{  rayO[0] - vertex0[0],
//                            rayO[1] - vertex0[1],
//                            rayO[2] - vertex0[2]};
//    u = f * dotProduct(s,h);
//    if (u < 0.0 || u > 1.0)
//        return false;
//    std::vector<double> q = crossProduct(s, edge1);
//    v = f * dotProduct(rayV, q);
//    if (v < 0.0 || u + v > 1.0)
//        return false;
//    // At this stage we can compute t to find out where the intersection point is on the line.
//    double t = f * dotProduct(edge2, q);
//    if (t > EPSILON) // ray intersection
//    {
//        outIntersectionPoint = rayOrigin + rayVector * t;
//        return true;
//    }
//    else // This means that there is a line intersection but not a ray intersection.
//        return false;
//}






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


// TODO: why can the Delaunay be const here? I'm changing the cell scores that are saved inside the Delaunay!
int traverseCells(Delaunay& Dt,
                  Cell_handle& current_cell, std::unordered_set<Cell_handle>& processed_cells,
                  Ray ray, Point source, double sigma, int oppositeVertex,
                  bool inside)
{
    // input::
    // &Delaunay            Reference to a Delaunay triangulation
    // ray                  Current ray
    // source               The "Delaunay" point of the current ray
    // current_cell         The cell that has just been entered (in the global context)
    // oppositeVertex       The opposite vertex of the facet where the current_cell was entered
    // inside               Bool that says if the current cell is before or after the point

    // check if the current cell was not already processed
    // this will only happen if this function is called by itself (not from the firstCell function)
    // and if the intersection point is a point of the cell
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

            Triangle tri = Dt.triangle(fac);
            Point intersectionPoint;
            Point rayO = ray.source();
            Vector rayV = ray.to_vector();
            bool result = rayTriangleIntersection(rayO, rayV, tri, intersectionPoint);
            std::pair<float,float> score;
            // check if there is an intersection between the current ray and current triangle
            if(result){
                // mark the current cell as processed
                processed_cells.insert(current_cell);
                // get the distance between the source of the ray and the intersection point with the current cell
                float dist2 = CGAL::squared_distance(intersectionPoint, rayO);
                // calculate the score for the current cell based on the distance
                score = cellScore(dist2, sigma, inside);
                // now locate the current cell in the global context of the triangulation,
                // so I can set the score
                current_cell->info().outside_score += score.first;
                current_cell->info().inside_score += score.second;

                // 2. get the neighbouring cell of the current triangle and check for ray triangle intersections in that cell
//                Facet fac = std::make_pair(current_cell, idx);
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
//                for(int i=0; i<4; i++){
//                    Point pt = newCell->vertex(i)->point();
//                    CGAL::cpp11::result_of<Intersect(Point, Ray)>::type
//                      point_intersection = intersection(pt, ray);


//                    if(point_intersection){
////                        const Point* p = boost::get<Point>(&*point_intersection);
////                        std::cout << "intersection with a vertex of the cell: " << *p << std::endl;
//                        return 0;
//                    }
//                }
                traverseCells(Dt, newCell, processed_cells, ray, source, sigma, newIdx, inside);
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

void firstCell(Delaunay& Dt, Delaunay::Finite_vertices_iterator& vit,
               bool inside,
               int& intersection_count){


    // get sigma of the current vertex
    double sigma = vit->info().sigma;

    // init a set with the already processed cells (for this ray!)
    std::unordered_set<Cell_handle> processed_cells;

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

            Point intersectionPoint;
            Point source = vit->point();
            Vector rayV = ray.to_vector();
            // check if there is an intersection between the current ray and current triangle
            bool oresult = rayTriangleIntersection(source, rayV, tri, intersectionPoint);
            std::pair<float,float> score;

            CGAL::cpp11::result_of<Intersect(Triangle, Ray)>::type
              result = intersection(tri, ray);


            if(oresult){
                intersection_count++;
//                if(const Point* p = boost::get<Point>(&*result)){
                    // get the distance between the source of the ray and the intersection point with the current cell
                    float dist2 = CGAL::squared_distance(intersectionPoint, source);
                    // calculate the score for the current cell based on the distance
                    score = cellScore(dist2, sigma, inside);
                    // now locate the current cell in the global context of the triangulation,
                    // so I can set the score
                    current_cell->info().outside_score += score.first;
                    current_cell->info().inside_score += score.second;
                    // mark the current cell as processed
                    processed_cells.insert(current_cell);

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
                        traverseCells(Dt, newCell, processed_cells, ray, source, sigma, newIdx, inside);
                    }
                // if there was a match, break the loop around this vertex, so it can go to the next one
                // this is only done for speed reasons, it shouldn't have any influence on the label
                    // because a ray can only hit more than one facet of a cell if it hits another point of the cell
                    // in which case it goes THROUGH a facet of the cell, in which case the intersection is not a point
                    // but an edge.
                    // it does however make a difference if this is turn on or not. why??
                    break;
                // here the if(boost point) ends...
//                }

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

void rayTracingFun(Delaunay& Dt){

    std::cout << "Start tracing rays to every point..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // from Efficient volumetric fusion paper, it follows that I need to sum over all rays to get the score for a tetrahedron
    // iterate over every vertices = iterate over every ray
    // TODO: go in here and map the following to a function
    // this function can than be called with the vft2 iterator
    // and I pass this iterator around and also pass it from within the traversal function
    // when I want to restart from a point that is hit
    int intersection_count = 0;
    Delaunay::Finite_vertices_iterator vit;
    for(vit = Dt.finite_vertices_begin() ; vit != Dt.finite_vertices_end() ; vit++){

        // collect outside votes
        firstCell(Dt, vit, 0, intersection_count);    // one_cell currently not used in the correct way
        // collect inside votes
        firstCell(Dt, vit, 1, intersection_count);    // one_cell currently not used in the correct way
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Ray tracing done in " << full_duration.count() << "s" << std::endl;
    std::cout << "intersection count " << intersection_count << std::endl;

}

// end of namespace rayTracing
}






