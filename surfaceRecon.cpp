//#include <cgal_typedefs.h>
//#include <fileIO.h>
//#include <rayTracing.h>
//#include <tetIntersection.h>

//#include "meshProcessing.cpp"
//#include "pointSetProcessing.cpp"
//#include "tetTracing.cpp"
////#include "tetTracingBB.cpp"
//#include "optimization.cpp"

#include <surfaceRecon.h>





//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void surfaceReconstruction(std::string file_number, double regularization_weight)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::string path1 = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string path1 = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";

    std::string ifn1 = path1+"musee/TLS/Est1.mesh_cut"+file_number;
    std::string ifn2 = path1+"musee/AP/AP_alligned_cut1";

//    std::string ifn1 = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/fused";
    std::string ofn = ifn1;
//    std::string ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/fused_mesh";

    ifn1+=".ply";
    ifn2+=".ply";

    // TLS
    std::vector<Point> t_points;
    std::vector<vertex_info> t_infos;
    std::vector<std::vector<int>> t_polys;
    readTLS(ifn1, t_points, t_infos, t_polys);
    // AP
//    std::vector<Point> a_points;
//    std::vector<vertex_info> a_infos;
//    readAP(ifn2, a_points, a_infos);
//    concatenateData(a_points, a_infos, t_points, t_infos, 1);
//    t_points = a_points;
//    t_infos = a_infos;

    Delaunay Dt = makeDelaunayWithInfo(t_points, t_infos);

    // calculate noise per point and save it in the vertex_info of the Dt
    double medianNoise = pcaKNN(Dt, t_points);
//    pcaDt(Dt);
    // TODO: calculate a sigma = sigmaKNN * sigmaDelaunay


    rayTracing::rayTracingFun(Dt, medianNoise);

//    int outside_weight = 2.2;
//    tetTracing::firstCell(Dt, t_polys, outside_weight);
//    tetTracingBB::treeIntersection(Dt, t_polys);

    /////// feature normalization
//    standardizeScores(Dt);
//    normalizeScores(Dt);
//    logScore(Dt);
//    softmax(Dt);


    // Dt, area_weight, iteration
    std::cout << "regularization weight: " << regularization_weight << std::endl;
    GeneralGraph_DArraySArraySpatVarying(Dt, regularization_weight, -1);


//    Delaunay::Finite_cells_iterator cit;
//    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
//        std::cout << cit->info().final_label << std::endl;
//    }
    // good area weight for fontaine dataset is 15.0, for daratec 0.01,

    // Dt, output_file, optimized, (pruned=1 or colored=0)
    std::vector<Point> remaining_points;
    std::vector<std::vector<int>> remaining_facets;
    exportSurfacePLY(Dt, remaining_points, remaining_facets, ofn, 0);
    exportSurfacePLY(Dt, remaining_points, remaining_facets, ofn, 1);
//    exportSurfacePLY(Dt, ofn, 0, 1);
//    exportSurfacePLY(Dt, ofn, 0, 0);
//    exportSurfacePLY(Dt, ofn, 1, 1);

    exportCellCenter(ofn, Dt);


//    // create surface mesh
//    Polyhedron out_mesh;
//    std::vector<std::vector<int>> polygons;
//    std::vector<Point> points;
//    // create surface mesh, with orient and nb_of_components_to_keep
//    createSurfaceMesh(Dt, points, polygons, out_mesh, 1, 1);
//    // export surface mesh as OFF
//    exportOFF(out_mesh, ofn);

//    // Quality control
////    double max_dist =
////          CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(out_mesh, points, 4000);
//    double max_dist =
//          CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(points, out_mesh);
//    std::cout << "Max distance to point_set: " << max_dist << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto full_duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Surface reconstruction done in " << full_duration.count() << "s" << std::endl;

}

// comments:
// unlabelled tets, i.e. that are missed by the ray are not such a big problem, as they can be fixed by the
// optimization.
// bigger problem are wrongly labelled tetrahedrons, i.e. mostly wrong inside labeles.
// optimization can also get rid of them, but at the cost of producing wholes, by pruning away valid surface triangles.
// the wholes are then filled in the back from the unlabelled tets in the back.
// problem now is that with a tetrahedron tracing I can get rid of the


// TODO:
// add a post processing step where I find non-manifold triangles and fix the region by applying face selection in the way that
// PolyFit does it. Use their Gurobi energy minimization solver for it (their energy cannot be minimized with a graph-cut).
// Problem will be that it might be to slow, and that it is not clear how to extract region with problem.

