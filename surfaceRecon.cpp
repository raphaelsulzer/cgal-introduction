#include <cgal_typedefs.h>

#include "fileIO.cpp"
#include "meshProcessing.cpp"
#include "pointSetProcessing.cpp"
#include "rayTracing.cpp"
#include "optimization.cpp"
#include "plyDefinition.cpp"


void readSensorMesh(std::string ofn, Mesh_ply aMesh){


    Polyhedron out_mesh;

    std::vector<Point> points;
    std::map<Point, std::pair<int, Vector>> pvm;
    for(int i = 0; i < aMesh.mVertices.size(); i++){
        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
        points.push_back(pt);
        Vector vec(aMesh.mvCapture[i].x - pt.x(), aMesh.mvCapture[i].y - pt.y(), aMesh.mvCapture[i].z - pt.z());
        pvm[pt] = std::make_pair(i, vec);
    }
    std::vector<std::vector<int>> polygons;
    for(int i = 0; i < aMesh.mIndices.size()/3; i++){
        std::vector<int> poly(3);
        poly[0] = aMesh.mIndices[(i*3)%3+0];
        poly[1] = aMesh.mIndices[(i*3)%3+1];
        poly[2] = aMesh.mIndices[(i*3)%3+2];
        polygons.push_back(poly);
    }

    std::cout<< "orientation: " << CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons) << std::endl;
//    std::cout << CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons) << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, out_mesh);


    //    decimateSurfaceMesh(input, output);


    std::fstream fo;
    fo.open(ofn, std::fstream::out);
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "element vertex " << out_mesh.size_of_vertices() << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property float nx" << std::endl;
    fo << "property float ny" << std::endl;
    fo << "property float nz" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(8);

    Polyhedron::Point_iterator pit;
//    int i = 0;
    for(pit = out_mesh.points_begin(); pit != out_mesh.points_end(); pit++){
        fo << *pit << " " << pvm.find(*pit)->second.second << std::endl;
//        std::cout << i << std::endl;
//        std::cout << *pit << " " << pvm.find(*pit)->second.second << std::endl;
//        i++;
    }

    fo.close();

    Polyhedron::Facet_iterator fit;
    for(fit = out_mesh.facets_begin(); fit < out_mesh.facets_end(); fit++){

        fit;

    }

    // TODO: save the sensor mesh as an off file, e.g. with Meshlab and read it into a surface mesh with CGAL
    // then do the tetrahedron - tetrahedron intersection
    // problem is that Meshlab and CC both mess up the OFF file.
    // another problem is that Meshlab cannot export the additional information (namely sensor center) and CC even messes up PLY files
    // this means I cannot read in the original PLY file, neither with Meshlab nor with CC and also not with CGAL
    // not with CGAL because the original PLY loader does not allow to read in the mesh structure, and in order to write one myself, the
    // file would have to be ASCII not Binary

}


//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int main()
{

    std::string path = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string path = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";

    std::string ifn = path+"musee/Est1.mesh_cut2";
//    std::string ifn = path+"daratech/daratech25000";
    std::string ofn = ifn;
    ifn+=".ply";


    // read ASCII PLY with normal
    std::vector<Point> a_points;
    std::vector<vertex_info> a_infos;
//    readPLY(ifn, a_points, a_infos);

    // read Binary PLY with sensor
    Mesh_ply aMesh;
    Import_PLY(ifn.c_str(), &aMesh);

    for(int i = 0; i < aMesh.mVertices.size(); i++){
        Point pt(aMesh.mVertices[i].x, aMesh.mVertices[i].y, aMesh.mVertices[i].z);
        a_points.push_back(pt);
        Vector vec(aMesh.mvCapture[i].x - pt.x(), aMesh.mvCapture[i].y - pt.y(), aMesh.mvCapture[i].z - pt.z());
        vertex_info vec_inf;
        vec_inf.sensor = vec;
        a_infos.push_back(vec_inf);
    }



    std::vector<std::vector<int>> sensor_triangle;
    for(int i = 0; i < aMesh.mIndices.size()/3; i++){
        std::vector<int> poly(3);
        poly[0] = aMesh.mIndices[(i*3)%3+0];
        poly[1] = aMesh.mIndices[(i*3)%3+1];
        poly[2] = aMesh.mIndices[(i*3)%3+2];
        sensor_triangle.push_back(poly);
    }


    // combine point clouds
//    std::vector<Point> t_points;
//    std::vector<vertex_info> t_infos;
//    readPLY(ifn2, t_points, t_infos);
//    ifn2+=".ply";

//    copyInfo(a_points, a_infos, t_points, t_infos);
//    a_points.insert(a_points.end(), t_points.begin(), t_points.end());
//    a_infos.insert(a_infos.end(), t_infos.begin(), t_infos.end());


    Delaunay Dt = makeDelaunayWithInfo(a_points, a_infos);

    iterateOverTetras(Dt, a_points, a_infos, sensor_triangle);

    // calculate noise per point and save it in the vertex_info of the Dt
    pcaKNN(Dt, a_points);
//    pcaDt(Dt);
    // TODO: calculate a sigma = sigmaKNN * sigmaDelaunay

    // 0 = camera, 1 = normal
    bool ray_construction = 1;

    // ray tracing for Dt for saving initial cell labels in cell info;
    // parameters: is one_cell traversel only.
    rayTracingFun(Dt, 1);

    // Dt, area_weight, iteration
    GeneralGraph_DArraySArraySpatVarying(Dt, 0.1, -1);
    // good area weight for fontaine dataset is 15.0, for daratec 0.01,

    // Dt, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    exportPLY(Dt, ofn, ray_construction, 1, 0);
    exportPLY(Dt, ofn, ray_construction, 1, 1);

//    // create surface mesh
//    Polyhedron out_mesh;
//    std::vector<std::vector<int>> polygons;
//    std::vector<Point> points;
//    // create surface mesh, with orient and nb_of_components_to_keep
//    createSurfaceMesh(Dt, points, polygons, out_mesh, 1, -1);
//    // export surface mesh as OFF
//    exportOFF(out_mesh, ofn);

//    // Quality control
////    double max_dist =
////          CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(out_mesh, points, 4000);
//    double max_dist =
//          CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(points, out_mesh);
//    std::cout << "Max distance to point_set: " << max_dist << std::endl;


    // TODO:
    // 1. do the full ray tracing to the outside, but not to the inside
    // simply replace the one_cell if-statement with one_cell && inside
    // 2. get the correct sensor orientation from COLMAP or from seperate depth map from each image from MicMac
    // 3. intersect a sensor topology tetrahedron (formed by 3 pixels next to each other, or LiDAR points next to each other and their (almost common -> barycenter) ray source
    // use this for outside vote of the Delaunay tetrahedra, and keep the ray for inside votes for now


    return 0;
}






