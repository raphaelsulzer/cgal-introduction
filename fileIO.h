#ifndef FILEIO_H
#define FILEIO_H

#include <cgal_typedefs.h>
//#include "plyDefinition.cpp"
//#include "colmapPLY.cpp"


//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
//////////////////////////////////////////////////////////
// generate a Delaunay triangulation from a PLY file
Delaunay triangulationFromFile();

void readBinaryPLY(std::string ifn,
                   std::vector<Point>& points, std::vector<vertex_info>& infos,
                   bool colmap);

void readBinaryPLY(std::string ifn,
                   std::vector<Point>& points, std::vector<vertex_info>& infos,
                   std::vector<std::vector<int>>& sensor_triangle,
                   bool colmap);

void readASCIIPLY(std::string ifn,
                  std::vector<Point>& points, std::vector<vertex_info>& infos);

void concatenateData(std::vector<Point>& a_points, std::vector<vertex_info>& a_info,
                     std::vector<Point>& t_points, std::vector<vertex_info>& t_info,
                     bool copyInfo);

/////////////////////////////////////////////////////////////////////
/////////////////////////////// OUTPUT //////////////////////////////
/////////////////////////////////////////////////////////////////////
void exportEdges();

void exportPoints();

void exportOFF(Polyhedron& out_mesh, std::string path);
void exportOFF(Polyhedron_Exact& out_mesh, std::string path);
void exportOFF(Tetrahedron& in_tet, std::string path);

void importOff(std::string path, Polyhedron& import_poly);
void importOff(std::string path, Tetrahedron& import_tet);
void importOff(std::string path, std::vector<Point>& points);
void importOff(std::string path, std::vector<Plane>& planes);


void fixSensorCenter();

void exportCellCenter(std::string path,
                      const Delaunay& Dt);

void exportPLY(const Delaunay& Dt,
               std::string path,
               bool normals, bool optimized, bool prune_or_color);



#endif // FILEIO_H
