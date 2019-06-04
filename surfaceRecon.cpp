#include <cgal_typedefs.h>
#include "fileIO.cpp"
#include "pointSetProcessing.cpp"
#include "rayTracing.cpp"
#include "optimization.cpp"







//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int main()
{

    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";
    std::string ifn = path+"fontaine/fontaine_10000_normals";
    std::string ofn = ifn;
    ifn+=".ply";

//    std::string ifn = path+"2cube_10000sampled_messyNormals.ply";
//    std::string ofn = path+"2cube_10000sampled_messyNormals";


//    std::string ifn = path+"cube5000.ply";
//    std::string ofn = path+"cube5000";

    Delaunay Dt = triangulationFromFile(ifn);

    VNC_map all_vertices;
    pca(Dt, all_vertices);

//    exportSimple(Dt, all_vertices, ofn+"_pca.ply");

    Cell_map all_cells;

    rayTracingFun(Dt, all_cells, all_vertices);

//    checkEnergyTerms(Dt, all_cells, 10.0);

    // Dt, all_cells, file_output, area_weight, iteration
    GeneralGraph_DArraySArraySpatVarying(Dt, all_cells, 15.0, -1);

    Polyhedron out_mesh;
    std::vector<std::vector<int>> polygons;
    std::vector<Point> points;
    //    std::vector<Polygon_3> tris;

    // Dt, all_cells, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    exportSoup(Dt, all_cells, ofn, 1, 1, 1, points, polygons);

    bool oriented = CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);


    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points,polygons,out_mesh);

    Kernel kernel;

    double max_dist = CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(points,polygons,0.35,kernel);


    // TODO: implement a quality check function that measures the distance between the mesh and the point cloud
    // replace the ray tracing - it takes to long.

    // export a Polyhedron surface mesh (as .OFF)
    unsigned int nb_components_to_keep = 2;
    unsigned int erased_components = out_mesh.keep_largest_connected_components(nb_components_to_keep);
    std::cout << "number of erased components: " << erased_components << std::endl;
    std::cout << "max dist to point set: " << max_dist << std::endl;
    std::ofstream out("/home/raphael/Dropbox/Studium/PhD/data/sampleData/tet-oriented1.off");
    out << out_mesh;
    out.close();

    return 0;
}






