#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <tuple>
#include <utility>
//#include "exportTri.cpp"
//#include "exportTriWithCn.cpp"
//#include "info_insert_with_zip_iterator.cpp"
#include "rayTriIntersection.cpp"


int main()
{


    const char* ifn = "/home/raphael/PhD_local/data/tanksAndTemples/Barn_COLMAP_subsampled.ply";
    const char* ofn = "/home/raphael/PhD_local/data/tanksAndTemples/Barn_COLMAP_ss_triangulated.ply";
    const char* ofn_test = "/home/raphael/PhD_local/data/tanksAndTemples/test.ply";
//    int result = exportTriangulationFun(ifn, ofn);


//    exportTriWithCnFun(ifn, ofn);
    rayTriIntersectionFun(ofn_test);

    return 0;




}
