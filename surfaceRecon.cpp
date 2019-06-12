#include <cgal_typedefs.h>
#include "fileIO.cpp"
#include "pointSetProcessing.cpp"
#include "rayTracing.cpp"
#include "optimization.cpp"

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/tags.h>


//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int main()
{

    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string ifn = path+"fontaine/fontaine_10000_normals";
    std::string ifn = path+"office/clouds/office_15000";
//    std::string ifn = path+"daratech/daratech25000";
//    std::string ifn = path+"musee/musee";
    std::string ofn = ifn;
    ifn+=".ply";

    // for reading file with normals
//    std::vector<PN> point_with_info;
//    std::vector<Vector> info;
    // for reading file with camera_index
//    std::vector<PC> point_with_info;
//    std::vector<int> info;
    // for reading file with camera_index
//    std::vector<PNC> point_with_info;
//    std::vector<NC> info;

    std::vector<PCN> ply_lines;
    Delaunay Dt = triangulationFromFile(ifn, ply_lines);

    VPS_map all_vertices;
    // populate the VPS_map all vertices with vertex, point and sigma for every vertex of the Dt
    pcaKNN(Dt, all_vertices);

//    exportSimple(Dt, all_vertices, ofn+"_pca.ply");

    Cell_map all_cells;

    // 0 = camera, 1 = normal
    bool ray_construction = 1;

    // ray tracing for Dt for saving initial cell labels in all_cells;
    // parameters: is one_cell traversel only.
    rayTracingFun(Dt, all_cells, all_vertices, 1);

//    checkEnergyTerms(Dt, all_cells, 10.0);

    // Dt, all_cells, file_output, area_weight, iteration
    GeneralGraph_DArraySArraySpatVarying(Dt, all_cells, 15.0, -1);
    // good area weight for fontaine dataset is 15.0, for daratec 0.01,

    // Dt, all_cells, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    exportPLY(Dt, all_cells, ofn, ray_construction, 1, 1);

    // create surface mesh
    Polyhedron out_mesh;
    std::vector<std::vector<int>> polygons;
    std::vector<Point> points;
    // create surface mesh, with orient and nb_of_components_to_keep
    createSurfaceMesh(Dt, all_cells, points, polygons, out_mesh, 1, -1);


//    double max_dist =
//          CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(out_mesh,
//                                                                   points,
//                                                                   4000);
    double max_dist =
          CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(points,
                                                                   out_mesh);

    std::cout << "Max distance to point_set: " << max_dist << std::endl;
    // export surface mesh as OFF
    exportOFF(out_mesh, ofn);


//    std::cout << "max dist to point set: " << max_dist << std::endl;


    return 0;
}






