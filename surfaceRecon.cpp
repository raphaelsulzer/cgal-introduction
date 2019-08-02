#include <cgal_typedefs.h>
#include <fileIO.h>
#include <rayTracing.h>
#include <tetIntersection.h>

#include "meshProcessing.cpp"
#include "pointSetProcessing.cpp"
#include "tetTracing.cpp"
#include "optimization.cpp"

//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void surfaceReconstruction()
{
    auto start = std::chrono::high_resolution_clock::now();


//    std::string path1 = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
    std::string path1 = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";


    std::string ifn1 = path1+"musee/TLS/Est1.mesh_cut2";
    std::string ifn2 = path1+"musee/AP/fused_fixedSensor_cut_alligned";     // there might be a problem with this file since it was exported as an ASCII from the CC

//    std::string ifn1 = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/fused";
    std::string ofn = ifn1;
//    std::string ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/fused_mesh";

    ifn1+=".ply";
    ifn2+=".ply";

//     read ASCII PLY with normal
    std::vector<Point> t_points;
    std::vector<vertex_info> t_infos;
    std::vector<std::vector<int>> t_polys;
    readBinaryPLY(ifn1, t_points, t_infos, t_polys, 0);
//    std::vector<Point> a_points;
//    std::vector<vertex_info> a_infos;
//    readASCIIPLY(ifn2, a_points, a_infos);
//    concatenateData(a_points, a_infos, t_points, t_infos, 0);

    Delaunay Dt = makeDelaunayWithInfo(t_points, t_infos);

    // calculate noise per point and save it in the vertex_info of the Dt
    pcaKNN(Dt, t_points);
//    pcaDt(Dt);
    // TODO: calculate a sigma = sigmaKNN * sigmaDelaunay


    // ray tracing for Dt for saving initial cell labels in cell info;
    // parameters: is one_cell traversel only.
    rayTracing::rayTracingFun(Dt);
    tetTracing::firstCell(Dt, t_polys);


    // Dt, area_weight, iteration
    GeneralGraph_DArraySArraySpatVarying(Dt, 0.01, -1);
    // good area weight for fontaine dataset is 15.0, for daratec 0.01,

    // Dt, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    // 0 = camera, 1 = normal
    bool ray_construction = 1;
    exportPLY(Dt, ofn, ray_construction, 1, 0);
    exportPLY(Dt, ofn, ray_construction, 1, 1);

    exportCellCenter(ofn, Dt);


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

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Surface reconstruction done in " << full_duration.count() << "s" << std::endl;


}






