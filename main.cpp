#include <cgal_typedefs.h>
#include <fileIO.h>
#include <tetIntersection.h>

#include "surfaceRecon.cpp"

int main(){

    //std::string path="", double regularization_weight=1
    surfaceReconstruction(0.1);

//    tetIntersectionTest();

//    int hits = 84223;
//    for(int i = 0; i < hits; i++){
//        std::cout<<i<<std::endl;
//        newTetIntersectionTest();
//    }

    return 0;

}
