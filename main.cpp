#include <cgal_typedefs.h>
#include <fileIO.h>
#include <tetIntersection.h>

#include "surfaceRecon.cpp"

#include <boost/dynamic_bitset.hpp>

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
    surfaceReconstruction(file_number, regularization_weight);

    // process colmap PLY file
//    readColmapPLY();

//    tetIntersectionTest();

//    int hits = 84223;
//    for(int i = 0; i < hits; i++){
//        std::cout<<i<<std::endl;
//        newTetIntersectionTest();
//    }

    return 0;

}
