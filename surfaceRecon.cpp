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
    std::string ifn = path+"fontaine/fontaine_200000_normals";
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

    // Dt, all_cells, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    exportSoup(Dt, all_cells, ofn, 1, 1, 0);

    // TODO: implement a quality check function that measures the distance between the mesh and the point cloud

    return 0;
}






