#ifndef SURFACERECON_H
#define SURFACERECON_H


#include <cgal_typedefs.h>
#include <fileIO.h>
#include <rayTracing.h>
#include <tetIntersection.h>

#include "meshProcessing.cpp"
#include "pointSetProcessing.cpp"
#include "tetTracing.cpp"
//#include "tetTracingBB.cpp"
#include "optimization.cpp"


//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void surfaceReconstruction(std::string path, double regularization_weight);










#endif // SURFACERECON_H
