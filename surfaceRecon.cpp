#include <cgal_typedefs.h>
#include <fileIO.h>
#include <rayTracing.h>
#include <tetIntersection.h>

#include "meshProcessing.cpp"
#include "pointSetProcessing.cpp"
#include "tetTracing.cpp"
//#include "tetTracingBB.cpp"
#include "optimization.cpp"


void standardizeScores(Delaunay& Dt){

    Delaunay::Finite_cells_iterator cit;
    std::vector<double> inside_scores;
    std::vector<double> outside_scores;
    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        inside_scores.push_back(cit->info().inside_score);
        outside_scores.push_back(cit->info().outside_score);
    }

//    double inside_scale = 255/(*std::max_element(inside_scores.begin(), inside_scores.end()));
//    double outside_scale = 255/(*std::max_element(outside_scores.begin(), outside_scores.end()));
    double inside_mean = std::accumulate(inside_scores.begin(), inside_scores.end(), 0.0)/inside_scores.size();
    double outside_mean = std::accumulate(outside_scores.begin(), outside_scores.end(), 0.0)/outside_scores.size();



    std::vector<double> diff1(inside_scores.size());
    std::transform(inside_scores.begin(),
                   inside_scores.end(),
                   diff1.begin(), [inside_mean](double x) { return x - inside_mean; });
    double sq_sum1 = std::inner_product(diff1.begin(), diff1.end(), diff1.begin(), 0.0);
    double inside_stdev = std::sqrt(sq_sum1 / inside_scores.size());

    std::vector<double> diff2(outside_scores.size());
    std::transform(outside_scores.begin(),
                   outside_scores.end(),
                   diff2.begin(), [outside_mean](double x) { return x - outside_mean; });
    double sq_sum2 = std::inner_product(diff2.begin(), diff2.end(), diff2.begin(), 0.0);
    double outside_stdev = std::sqrt(sq_sum2 / outside_scores.size());

    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
//        cit->info().inside_score = (cit->info().inside_score - inside_mean) / inside_stdev + 2;
//        cit->info().outside_score = (cit->info().outside_score - outside_mean) / outside_stdev + 2;
        cit->info().inside_score;
        cit->info().outside_score;
    }

}

void normalizeScores(Delaunay& Dt){

    Delaunay::Finite_cells_iterator cit;
    std::vector<double> inside_scores;
    std::vector<double> outside_scores;
    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        inside_scores.push_back(cit->info().inside_score);
        outside_scores.push_back(cit->info().outside_score);
    }

    double inside_min = *std::min_element(inside_scores.begin(), inside_scores.end());
    double outside_min = *std::min_element(outside_scores.begin(), outside_scores.end());
    double inside_max = *std::max_element(inside_scores.begin(), inside_scores.end());
    double outside_max = *std::max_element(outside_scores.begin(), outside_scores.end());



    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){
        cit->info().inside_score = (cit->info().inside_score - inside_min) / (inside_max - inside_min);
        cit->info().outside_score = (cit->info().outside_score - outside_min) / (outside_max - outside_min);
//        cit->info().inside_score;
//        cit->info().outside_score;
    }

}




//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void surfaceReconstruction()
{
    auto start = std::chrono::high_resolution_clock::now();

    std::string path1 = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string path1 = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";

    std::string ifn1 = path1+"musee/TLS/Est1.mesh_cut4";
    std::string ifn2 = path1+"musee/AP/fused_fixedSensor_cut_alligned";     // there might be a problem with this file since it was exported as an ASCII from the CC

//    std::string ifn1 = "/home/raphael/PhD_local/data/museeZoologic/aerial_images/BIOM-EMS/colmap/results/fused";
    std::string ofn = ifn1;
//    std::string ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/fused_mesh";

    ifn1+=".ply";
    ifn2+=".ply";

//     read ASCII PLY with normal
    std::vector<Point> t_points;
    std::vector<vertex_info> t_infos;
    std::vector<std::vector<int>> t_polys;
    readBinaryPLY(ifn1, t_points, t_infos, t_polys, 0);
//    std::vector<Point> a_points;
//    std::vector<vertex_info> a_infos;
//    readASCIIPLY(ifn2, a_points, a_infos);
//    concatenateData(a_points, a_infos, t_points, t_infos, 0);

    Delaunay Dt = makeDelaunayWithInfo(t_points, t_infos);

    // calculate noise per point and save it in the vertex_info of the Dt
//    pcaKNN(Dt, t_points);
//    pcaDt(Dt);
    // TODO: calculate a sigma = sigmaKNN * sigmaDelaunay


    // ray tracing for Dt for saving initial cell labels in cell info;
    // parameters: is one_cell traversel only.
//    bool full=true;
    rayTracing::rayTracingFun(Dt);
    tetTracing::firstCell(Dt, t_polys);
//    tetTracingBB::treeIntersection(Dt, t_polys);

    standardizeScores(Dt);
//    normalizeScores(Dt);

    // Dt, area_weight, iteration
    GeneralGraph_DArraySArraySpatVarying(Dt, 150, -1);
    // good area weight for fontaine dataset is 15.0, for daratec 0.01,

    // Dt, file_output, (normals=1 or cam_index=0), optimized, (pruned=1 or colored=0)
    exportPLY(Dt, ofn, 0, 1);
    exportPLY(Dt, ofn, 1, 1);

    exportCellCenter(ofn, Dt);


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




