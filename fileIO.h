#ifndef FILEIO_H
#define FILEIO_H

#include <cgal_typedefs.h>



//////////////////////////////////////////////////////////
///////////////////// FILE I/O ///////////////////////////
//////////////////////////////////////////////////////////
// generate a Delaunay triangulation from a PLY file
Delaunay triangulationFromFile();

void readColmapPLY();

void readAP(std::string ifn,
                   std::vector<Point>& points, std::vector<vertex_info>& infos);

void readTLS(std::string ifn,
                   std::vector<Point>& points, std::vector<vertex_info>& infos,
                   std::vector<std::vector<int>>& sensor_triangle);

void readASCIIPLY(std::string ifn,
                  std::vector<Point>& points, std::vector<vertex_info>& infos);

void concatenateData(std::vector<Point>& a_points, std::vector<vertex_info>& a_info,
                     std::vector<Point>& t_points, std::vector<vertex_info>& t_info,
                     bool copyInfo);

/////////////////////////////////////////////////////////////////////
/////////////////////////////// OUTPUT //////////////////////////////
/////////////////////////////////////////////////////////////////////


void printPLYHeader(std::fstream& fo,
                    int nv, int nf,
                    bool normals, bool color, bool sensor, bool cam_index, bool score,
                    bool fcolor,
                    int precision);

void exportEdges();

void exportPoints(std::string path, std::vector<Point>& points, std::vector<vertex_info>& infos);

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

void exportColoredFacetsPLY(const Delaunay& Dt,
                            std::string path,
                            bool optimized);

void exportSurfacePLY(const Delaunay& Dt,
                std::vector<Point>& remaining_points,
                std::vector<std::vector<int>>& remaining_polygons,
                std::string path,
                bool optimized,
                bool fixedManifold);

#endif // FILEIO_H
