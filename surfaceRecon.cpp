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

//    std::string ifn = path+"fontaine/fontaine_10000_normals";
//    std::string ifn = path+"office/clouds/office_15000";
    std::string ifn = path+"daratech/daratech25000";
//    std::string ifn = path+"musee/museeAP_05m";
//    std::string ifn2 = path+"musee/TLS_registered";
    std::string ofn = ifn;
    ifn+=".ply";


    std::vector<Point> a_points;
    std::vector<vertex_info> a_infos;
    readPLY(ifn, a_points, a_infos);

//    std::vector<Point> t_points;
//    std::vector<vertex_info> t_infos;
//    readPLY(ifn2, t_points, t_infos);
//    ifn2+=".ply";

    // combine point clouds
//    copyInfo(a_points, a_infos, t_points, t_infos);
//    a_points.insert(a_points.end(), t_points.begin(), t_points.end());
//    a_infos.insert(a_infos.end(), t_infos.begin(), t_infos.end());



//    // init a vector that has the correct elements of the PLY file, so either PNC or PN
//    std::vector<PNC> ply_lines;
//    Delaunay Dt = triangulationFromFile(ifn, ply_lines);


    Delaunay Dt = makeDelaunayWithInfo(a_points, a_infos);

    // calculate noise per point and save it in the vertex_info of the Dt
    pcaKNN(Dt, a_points);
//    pcaDt(Dt);
    // TODO: calculate a sigma = sigmaKNN * sigmaDelaunay

    // 0 = camera, 1 = normal
    bool ray_construction = 1;

    // ray tracing for Dt for saving initial cell labels in cell info;
    // parameters: is one_cell traversel only.
    rayTracingFun(Dt, 1);


//    checkEnergyTerms(Dt, all_cells, 10.0);

    // Dt, area_weight, iteration
    GeneralGraph_DArraySArraySpatVarying(Dt, 0.1, -1);
    // good area weight for fontaine dataset is 15.0, for daratec 0.01,

    // Dt, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    exportPLY(Dt, ofn, ray_construction, 0, 0);
    exportPLY(Dt, ofn, ray_construction, 0, 1);

//    // create surface mesh
//    Polyhedron out_mesh;
//    std::vector<std::vector<int>> polygons;
//    std::vector<Point> points;
//    // create surface mesh, with orient and nb_of_components_to_keep
//    createSurfaceMesh(Dt, points, polygons, out_mesh, 1, -1);
//    // export surface mesh as OFF
//    exportOFF(out_mesh, ofn);

//    // Quality control
////    double max_dist =
////          CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(out_mesh, points, 4000);
//    double max_dist =
//          CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(points, out_mesh);
//    std::cout << "Max distance to point_set: " << max_dist << std::endl;


    // TODO:
    // 1. do the full ray tracing to the outside, but not to the inside
    // simply replace the one_cell if-statement with one_cell && inside
    // 2. get the correct sensor orientation from COLMAP or from seperate depth map from each image from MicMac
    // 3. intersect a sensor topology tetrahedron (formed by 3 pixels next to each other, or LiDAR points next to each other and their (almost common -> barycenter) ray source
    // use this for outside vote of the Delaunay tetrahedra, and keep the ray for inside votes for now


    return 0;
}






