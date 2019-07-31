#include <cgal_typedefs.h>
#include <fileIO.h>
#include <tetIntersection.h>

#include "surfaceRecon.cpp"

int main(){

    surfaceReconstruction();

//    tetIntersectionTest();

//    int hits = 84223;
//    for(int i = 0; i < hits; i++){
//        std::cout<<i<<std::endl;
//        newTetIntersectionTest();
//    }

    return 0;

}
