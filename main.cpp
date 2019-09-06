#include <cgal_typedefs.h>
#include <fileIO.h>
#include <tetIntersection.h>

#include "surfaceRecon.cpp"

#include <boost/dynamic_bitset.hpp>

//#include <mvs/workspace.h>

#include "importDepthMap.cpp"

int main(int argc, char const *argv[]){

    std::string file_number;
    double regularization_weight;
    if(argc<3){
        file_number = "4";
        regularization_weight = 0.1;
    }
    else{
        file_number = argv[1];
        regularization_weight = atof(argv[2]);
    }
    std::cout << "Surface reconstruction of file number: " << file_number << " with regularization weight " << regularization_weight << std::endl;
//    surfaceReconstruction(file_number, regularization_weight);

    std::string tiff_path = "/home/raphael/PhD_local/data/fountain/micmac/images/PIMs-QuickMac/Nuage-Depth-0000.png_Prof.tif";
    readTiff(tiff_path);




    // process colmap PLY file
//    readColmapPLY();

//    tetIntersectionTest();

//    int hits = 84223;
//    for(int i = 0; i < hits; i++){
//        std::cout<<i<<std::endl;
//        newTetIntersectionTest();
//    }

//    auto cw_options = colmap::mvs::Workspace::Options();
//    cw_options.workspace_path = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/reconstruction/sfm_result";

//    auto cw_cached_image = colmap::mvs::Workspace::CachedImage();
//    auto cw = colmap::mvs::Workspace(cw_options);



    return 0;

}
