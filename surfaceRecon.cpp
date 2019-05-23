#include <cgal_typedefs.h>
#include "fileIO.cpp"
#include "rayTracing.cpp"
#include "optimization.cpp"

//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int main()
{
//    const char* ifn = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/3cube_10000sampled_messyNormals.ply";
//    const char* ofn = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/3cube_CGAL_pruned.ply";
//    const char* ifn = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/fontaine_10000.ply";
//    const char* ofn = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/fontaine_pruned.ply";

//    const char* ifn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/3cube_10000sampled_messyNormals.ply";
//    const char* ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/3cube_CGAL_pruned.ply";
//    const char* ifn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/fontaine_10000.ply";
//    const char* ifn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/fontaine_10000_messyNormals.ply";
//    const char* ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/fontaine_normals_pruned.ply";
//    const char* ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/fontaine_cameras_pruned.ply";
//    const char* ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/fontaine_initial_colored.ply";

    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";
    std::string ifn = path+"fontaine_10000_messyNormals.ply";
    std::string ofn = path+"fontaine_10000_messyNormals";

    Delaunay Dt = triangulationFromFile(ifn);

    Cell_map all_cells;

    rayTracingFun(Dt, all_cells);

    // Dt, all_cells, file_output, area_weight, iteration
    GeneralGraph_DArraySArraySpatVarying(Dt, all_cells, 10.0, -1);

    // Dt, all_cells, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    exportSoup(Dt, all_cells, ofn, 1, 1, 1);

    return 0;
}






