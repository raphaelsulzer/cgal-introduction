#include <cgal_typedefs.h>
#include <fileIO.h>

#include "surfaceRecon.cpp"
#include "newTetIntersection.cpp"

int main(){

    surfaceReconstruction();

//    int hits = 84223;
//    for(int i = 0; i < hits; i++){
//        std::cout<<i<<std::endl;
//        newTetIntersectionTest();
//    }

    return 0;

}
