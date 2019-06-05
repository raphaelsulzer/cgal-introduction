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
    std::string ifn = path+"daratech/daratech25000";
    std::string ofn = ifn;
    ifn+=".ply";



    Delaunay Dt = triangulationFromFile(ifn);

    VPS_map all_vertices;
    // populate the VPS_map all vertices with vertex, point and sigma for every vertex of the Dt
    pca(Dt, all_vertices);

//    exportSimple(Dt, all_vertices, ofn+"_pca.ply");

    Cell_map all_cells;

    // ray tracing for Dt for saving initial cell labels in all_cells; last parameter is one_cell traversel only
    rayTracingFun(Dt, all_cells, all_vertices, 0);

//    checkEnergyTerms(Dt, all_cells, 10.0);

    // Dt, all_cells, file_output, area_weight, iteration
    GeneralGraph_DArraySArraySpatVarying(Dt, all_cells, 0.1, -1);
    // good area weight for fontaine dataset is 15.0

    // Dt, all_cells, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    exportPLY(Dt, all_cells, ofn, 1, 0, 0);

//    // create surface mesh
//    Polyhedron out_mesh;
//    std::vector<std::vector<int>> polygons;
//    std::vector<Point> points;
//    // create surface mesh, with orient and nb_of_components_to_keep
//    createSurfaceMesh(Dt, all_cells, points, polygons, out_mesh, 1, -1);

//    // export surface mesh as OFF
//    exportOFF(out_mesh, ofn);


//    std::cout << "max dist to point set: " << max_dist << std::endl;


    return 0;
}






