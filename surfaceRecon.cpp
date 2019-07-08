#include <cgal_typedefs.h>

#include "fileIO.cpp"
#include "meshProcessing.cpp"
#include "rayTracing.cpp"
#include "optimization.cpp"

#include <CGAL/IO/PLY_reader.h>

#include <stdio.h>
#include "rply.h"
//#include "rply.c"

static int vertex_cb(p_ply_argument argument) {
    long eol;
    ply_get_argument_user_data(argument, NULL, &eol);
    printf("%g", ply_get_argument_value(argument));
    if (eol) printf("\n");
    else printf(" ");
    return 1;
}

static int face_cb(p_ply_argument argument) {
    long length, value_index;
    ply_get_argument_property(argument, NULL, &length, &value_index);
    switch (value_index) {
        case 0:
        case 1:
            printf("%g ", ply_get_argument_value(argument));
            break;
        case 2:
            printf("%g\n", ply_get_argument_value(argument));
            break;
        default:
            break;
    }
    return 1;
}


//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int main()
{

    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";

//    std::string ifn = path+"fontaine/fontaine_10000_normals";
//    std::string ifn = path+"office/clouds/office_15000";
//    std::string ifn = path+"daratech/daratech25000";
//    std::string ifn = path+"musee/museeAP_05m";
    std::string ifn = path+"musee/Est1.mesh_cut";
//    std::string ifn2 = path+"musee/TLS_registered";
    std::string ofn = ifn;
    ifn+=".ply";

//    Mesh_ply aMesh;
//    Import_PLY("/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/Est1.mesh_cut.ply", &aMesh);

//    std::string input = "/home/raphael/PhD_local/data/museeZoologic/TLS/EMSMesh/Est1.mesh.off";
//    std::string output = "/home/raphael/PhD_local/data/museeZoologic/TLS/EMSMesh/Est1.mesh_decimated.off";

//    decimateSurfaceMesh(input, output);

//    std::vector<Point> a_points;
//    std::vector<vertex_info> a_infos;
//    readPLYWithSensor(ifn, a_points, a_infos);

//    fixSensorCenter(ofn, a_points, a_infos);

    long nvertices, ntriangles;
    p_ply ply = ply_open("/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/Est1.mesh_cut.ply", 0);
    ply_read_header(ply);
    nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);
    ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 0);
    ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 1);
    ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", face_cb, NULL, 0);
    printf("%ld\n%ld\n", nvertices, ntriangles);
    ply_read(ply);
    ply_close(ply);





//    // else try polygon soup
//    std::vector<Kernel::Point_3> points;
//    std::vector<std::vector<std::size_t> > polygons;
//    std::vector<CGAL::Color> fcolors;
//    std::vector<CGAL::Color> vcolors;

//    std::ifstream in(ifn);
//    CGAL::read_PLY(in, points, polygons, fcolors, vcolors);

//    std::vector<std::vector<std::size_t> > polygons_sample;
//    for(std::size_t i = 0; i < 10; i++)
//    {
//        std::cout << polygons[i][0] << std::endl;

//    }

//    CGAL::read_PLY(in, points, polygons);




    // TODO: save the sensor mesh as an off file, e.g. with Meshlab and read it into a surface mesh with CGAL
    // then do the tetrahedron - tetrahedron intersection
    // problem is that Meshlab and CC both mess up the OFF file.
    // another problem is that Meshlab cannot export the additional information (namely sensor center) and CC even messes up PLY files
    // this means I cannot read in the original PLY file, neither with Meshlab nor with CC and also not with CGAL
    // not with CGAL because the original PLY loader does not allow to read in the mesh structure, and in order to write one myself, the
    // file would have to be ASCII not Binary


////    std::vector<Point> t_points;
////    std::vector<vertex_info> t_infos;
////    readPLY(ifn2, t_points, t_infos);
////    ifn2+=".ply";

//    // combine point clouds
////    copyInfo(a_points, a_infos, t_points, t_infos);
////    a_points.insert(a_points.end(), t_points.begin(), t_points.end());
////    a_infos.insert(a_infos.end(), t_infos.begin(), t_infos.end());



////    // init a vector that has the correct elements of the PLY file, so either PNC or PN
////    std::vector<PNC> ply_lines;
////    Delaunay Dt = triangulationFromFile(ifn, ply_lines);


//    Delaunay Dt = makeDelaunayWithInfo(a_points, a_infos);

//    // calculate noise per point and save it in the vertex_info of the Dt
//    pcaKNN(Dt, a_points);
////    pcaDt(Dt);
//    // TODO: calculate a sigma = sigmaKNN * sigmaDelaunay

//    // 0 = camera, 1 = normal
//    bool ray_construction = 1;

//    // ray tracing for Dt for saving initial cell labels in cell info;
//    // parameters: is one_cell traversel only.
//    rayTracingFun(Dt, 1);


////    checkEnergyTerms(Dt, all_cells, 10.0);

//    // Dt, area_weight, iteration
//    GeneralGraph_DArraySArraySpatVarying(Dt, 0.1, -1);
//    // good area weight for fontaine dataset is 15.0, for daratec 0.01,

//    // Dt, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
//    exportPLY(Dt, ofn, ray_construction, 0, 0);
//    exportPLY(Dt, ofn, ray_construction, 0, 1);

////    // create surface mesh
////    Polyhedron out_mesh;
////    std::vector<std::vector<int>> polygons;
////    std::vector<Point> points;
////    // create surface mesh, with orient and nb_of_components_to_keep
////    createSurfaceMesh(Dt, points, polygons, out_mesh, 1, -1);
////    // export surface mesh as OFF
////    exportOFF(out_mesh, ofn);

////    // Quality control
//////    double max_dist =
//////          CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(out_mesh, points, 4000);
////    double max_dist =
////          CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(points, out_mesh);
////    std::cout << "Max distance to point_set: " << max_dist << std::endl;


//    // TODO:
//    // 1. do the full ray tracing to the outside, but not to the inside
//    // simply replace the one_cell if-statement with one_cell && inside
//    // 2. get the correct sensor orientation from COLMAP or from seperate depth map from each image from MicMac
//    // 3. intersect a sensor topology tetrahedron (formed by 3 pixels next to each other, or LiDAR points next to each other and their (almost common -> barycenter) ray source
//    // use this for outside vote of the Delaunay tetrahedra, and keep the ray for inside votes for now


    return 0;
}






