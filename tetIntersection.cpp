#include <cgal_typedefs.h>
#include <fileio.h>


void tetIntersectionFun(Polyhedron& a, Polyhedron& b){

    Polyhedron::Facet_iterator sfi;
    std::vector<Triangle> a_triangles;
    for(sfi = a.facets_begin(); sfi != a.facets_end(); sfi++){
        std::cout << "new face" << std::endl;
        Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
        std::vector<Point> triangle_points;
        do{
            std::cout << circ->vertex()->point() << std::endl;
            triangle_points.push_back(circ->vertex()->point());
        }
        while (++circ != sfi->facet_begin());
        Triangle tri(triangle_points[0], triangle_points[1], triangle_points[2]);
        a_triangles.push_back(tri);
    }

    std::vector<Triangle> b_triangles;
    for(sfi = b.facets_begin(); sfi != b.facets_end(); sfi++){
        std::cout << "new face" << std::endl;
        Polyhedron::Halfedge_around_facet_circulator circ = sfi->facet_begin();
        std::vector<Point> triangle_points;
        do{
            std::cout << circ->vertex()->point() << std::endl;
            triangle_points.push_back(circ->vertex()->point());
        }
        while (++circ != sfi->facet_begin());
        Triangle tri(triangle_points[0], triangle_points[1], triangle_points[2]);
        b_triangles.push_back(tri);
    }

    // double loop over the triangles of each polyhedron
    std::set<Point> intersection_points;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){

            CGAL::cpp11::result_of<Intersect(Triangle, Triangle)>::type
              result = intersection(a_triangles[i], b_triangles[j]);

            if (result){
                if (const Segment* seg = boost::get<Segment>(&*result)){
                    intersection_points.insert(seg->point(0));
                    intersection_points.insert(seg->point(1));
                }
            }
        }
    }

    double vol;
    Polyhedron intersection_poly;
    if(intersection_points.size() > 3){
        CGAL::convex_hull_3(intersection_points.begin(), intersection_points.end(), intersection_poly);
        vol = CGAL::Polygon_mesh_processing::volume(intersection_poly);

    }
    else{
        vol = 0.0;
    }
    std::cout << vol << std::endl;

//    exportOFF(a, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/a");
//    exportOFF(b, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/b");
//    exportOFF(c, "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/tetras/c");

}
