#include <surfaceRecon.h>

//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// SURFACE RECONSTRUCTION ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void surfaceReconstruction(std::string file_number, double regularization_weight)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::string path1 = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
//    std::string path1 = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";

    std::string ifn1 = path1+"musee/TLS/Est1.mesh_cut"+file_number;
    std::string ifn2 = path1+"musee/AP/AP_alligned_cut1";

    std::string ofn = ifn1;
//    std::string ofn = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/musee/fused_mesh";

    ifn1+=".ply";
    ifn2+=".ply";

    // TLS
    std::vector<Point> t_points;
    std::vector<vertex_info> t_infos;
    std::vector<std::vector<int>> t_polys;
    readTLS(ifn1, t_points, t_infos, t_polys);
//    // AP
//    std::vector<Point> a_points;
//    std::vector<vertex_info> a_infos;
//    readAP(ifn2, a_points, a_infos);
//    concatenateData(a_points, a_infos, t_points, t_infos, 0);
//    t_points = a_points;
//    t_infos = a_infos;

    Delaunay Dt = makeDelaunayWithInfo(t_points, t_infos);


    ///// calculate noise per point and save it in the vertex_info of the Dt
//    double medianNoise = pcaKNN(Dt, t_points);
    double medianNoise = 0;
//    pcaDt(Dt);
    // TODO: calculate a sigma = sigmaKNN * sigmaDelaunay


    double infiniteScore = 10e3;

//    ///// Tet tracing
//    int outside_weight = 2.2;
//    tetTracing::firstCell(Dt, t_polys, outside_weight, infiniteScore);
////    tetTracingBB::treeIntersection(Dt, t_polys);
//    bool tetTracing = true;
    bool tetTracing = false;


    //// Ray tracing
    rayTracing::rayTracingFun(Dt, tetTracing, infiniteScore);


    /////// feature normalization
//    standardizeScores(Dt);
//    normalizeScores(Dt);
//    logScore(Dt);
//    softmax(Dt);

    ///// optimization: Dt, area_weight, iteration (-1 means until convergence)
    GeneralGraph_DArraySArraySpatVarying(Dt, regularization_weight, -1);

    ///// manifoldness operations
    std::vector<Point> remaining_points;
    std::vector<std::vector<int>> remaining_facets;
    exportSurfacePLY(Dt, remaining_points, remaining_facets, ofn, 1, 0);
    fixNonManifoldEdges(Dt, regularization_weight);
    fixNonManifoldPoints(Dt, regularization_weight);


    ////// Dt, ..., ..., output_file, optimized, manifold

//    std::vector<Point> remaining_points;
//    std::vector<std::vector<int>> remaining_facets;
//    exportSurfacePLY(Dt, remaining_points, remaining_facets, ofn, 0, 0);
    exportSurfacePLY(Dt, remaining_points, remaining_facets, ofn, 1, 1);

    exportColoredFacetsPLY(Dt, ofn, 1);
    exportCellCenter(ofn, Dt);

    /////// create surface mesh

    Polyhedron out_mesh;
    bool oriented = CGAL::Polygon_mesh_processing::orient_polygon_soup(remaining_points, remaining_facets);
    std::cout << "Created duplicate vertices to generated polygon surface? " << oriented << std::endl;
    std::cout << "... 0 meaning additional vertices were added for orientation (and probably for ensuring manifoldness)." << std::endl;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(remaining_points, remaining_facets, out_mesh);
    int erased_components = out_mesh.keep_largest_connected_components(1);
    std::cout << "Number of erased components from surface mesh: " << erased_components << std::endl;

    exportOFF(out_mesh, ofn);

    ///// Quality control
    double max_dist_to_point_set =
          CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set(out_mesh, t_points, 4000);
    std::cout << "Max distance to point set (precision): " << max_dist_to_point_set << std::endl;
    double max_dist_to_triangle_mesh =
          CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(t_points, out_mesh);
    std::cout << "Max distance to tiangle mesh (recall): " << max_dist_to_triangle_mesh << std::endl;

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
// 1.
// add a post processing step where I find non-manifold triangles and fix the region by applying face selection in the way that
// PolyFit does it. Use their Gurobi energy minimization solver for it (their energy cannot be minimized with a graph-cut).
// Problem will be that it might be to slow, and that it is not clear how to extract region with problem.
// problem: need to incrementally build the polygon surfaces for that to get a "controlled" facet ID with which I can query the neighbouring
// cells to get their label


// 2.
// make a synthetic dataset with sensor information
// label tets
// add noise to the dataset
// label tets again, and see what problems come with noise and how they can be fixed.




