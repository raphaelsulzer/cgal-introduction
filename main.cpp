#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <tuple>
#include <utility>
#include "exportPLY.cpp"
#include "readPlyWithCn.cpp"


int main()
{


    const char* ifn = "/home/raphael/PhD_local/data/tanksAndTemples/Barn_COLMAP_subsampled.ply";
    const char* ofn = "/home/raphael/PhD_local/data/tanksAndTemples/Barn_COLMAP_ss_triangulated.ply";
//    int result = exportTriangulationFun(ifn, ofn);

    auto ply = readPlyWithCnFun(ifn);

    auto points = ply.size();

    std::cout << points << std::endl;


    return 0;


}
