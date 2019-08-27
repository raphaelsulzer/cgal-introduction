#include <surfaceRecon.h>

typedef std::pair<std::pair<Cell_handle, int>, double> Facet_score;

void nonManifoldEdges(Delaunay& Dt, std::vector<std::vector<Facet_score>>& problematic_facets_per_edge){


    Delaunay::Finite_edges_iterator fei;
//    std::vector<std::vector<std::tuple<Cell_handle, Cell_handle>>> problematic_facets_per_edge;
    for(fei = Dt.finite_edges_begin(); fei != Dt.finite_edges_end(); fei++){

        Delaunay::Facet_circulator fac = Dt.incident_facets(*fei);
        Facet f1 = *fac;
//        std::vector<std::pair<Cell_handle,Cell_handle>> problematic_facets;
        std::vector<Facet_score> problematic_facets;
        do{
            Cell_handle c1 = fac->first;
            Facet mf = Dt.mirror_facet(*fac);
            Cell_handle c2 = mf.first;
            if(c1->info().final_label != c2->info().final_label){
                double score = sqrt(pow(c1->info().inside_score - c2->info().outside_score,2) +
                                    pow(c2->info().inside_score - c1->info().outside_score,2));
//                problematic_facets.push_back(std::make_tuple(c1,c2));
                problematic_facets.push_back(std::make_pair(std::make_pair(c1,fac->second),score));
            }
            fac++;
        }
        while(*fac != f1);
        if(problematic_facets.size()>2)
            problematic_facets_per_edge.push_back(problematic_facets);

    }
    std::cout << "number of non manifold edges: " << problematic_facets_per_edge.size() << std::endl;

}


bool sortbysec(const Facet_score &a,
              const Facet_score &b)
{
    return (a.second < b.second);
}



void selectFacets(Delaunay& Dt,
                  std::vector<std::vector<Facet_score>>& problematic_facets_per_edge,
                  std::string path){


    std::unordered_set<Vertex_handle> points;
    std::vector<Facet> remaining_facets;

    for(int i = 0; i < problematic_facets_per_edge.size(); i++){

        std::vector<Facet_score> problematic_facets;
        problematic_facets = problematic_facets_per_edge[i];
        std::sort(problematic_facets.begin(), problematic_facets.end(), sortbysec);

        // save the first two facets with the lowest score:
        Cell_handle c1 = problematic_facets[0].first.first;
        int vidx = problematic_facets[0].first.second;
        for(int j = vidx + 1; j <= vidx + 3; j++){
            // so c->vertex() gives me the global vertex handle from the Dt
            points.insert(c1->vertex(j%4));
        }
        remaining_facets.push_back(problematic_facets[0].first);

        c1 = problematic_facets[1].first.first;
        vidx = problematic_facets[1].first.second;
        for(int j = vidx + 1; j <= vidx + 3; j++){
            // so c->vertex() gives me the global vertex handle from the Dt
            points.insert(c1->vertex(j%4));
        }
        remaining_facets.push_back(problematic_facets[1].first);

//        Facet f2 = problematic_facets[1].first;
//        points.insert(f2.first->vertex((f2.second+1)%3));
//        points.insert(f2.first->vertex((f2.second+2)%3));
//        points.insert(f2.first->vertex((f2.second+3)%3));
//        remaining_facets.push_back(f2);

    }

    std::fstream fo;
    fo.open(path+"_selected_facets.ply", std::fstream::out);
    int nv = points.size();
    int nf = remaining_facets.size();

    printPLYHeader(fo,
                   nv, nf,
                   true, true, false, false, true,
                   15);

    std::unordered_set<Vertex_handle>::iterator vit;
    int vidx = 0;
    for(vit = points.begin(); vit != points.end(); vit++){
        // reset the vertex index here, because I need to know the order of exactly this loop here
        // for the indexing of the facets in the PLY file
        (*vit)->info().idx = vidx++;
        // also create a remaining_points vector
        // print data to file
        // coordinates
        fo  << (*vit)->point() << " "
        // color
            << int((*vit)->info().color[0]) <<  " " << int((*vit)->info().color[1]) <<  " " << int((*vit)->info().color[2]) <<  " "
        // normal
            << (*vit)->info().sensor_vec
        // endl
            << std::endl;

    }

    std::vector<Facet>::iterator fit;
    for(fit = remaining_facets.begin(); fit != remaining_facets.end(); fit++){

        // start printed facet line with a 3
        fo << 3 << ' ';
        Cell_handle c = fit->first;
        int vidx = fit->second;
        for(int j = vidx + 1; j <= vidx + 3; j++){
            // so c->vertex() gives me the global vertex handle from the Dt
            fo << c->vertex(j%4)->info().idx << ' ';
        }
//        for(int j = 1; j < 4; j++){
//            // print the indicies of each cell to the file
//            fo << fit->first->vertex((fit->second+j)%3)->info().idx << ' ';
//        }
        // put a color on the face, so that in Meshlab I can activate the color by facet mode, to compare with the "colored facet file"
        fo << "0 200 0";
        fo << std::endl;
    }
    fo.close();


}




//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// SURFACE RECONSTRUCTION ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void surfaceReconstruction(std::string file_number, double regularization_weight)
{
    auto start = std::chrono::high_resolution_clock::now();

//    std::string path1 = "/home/raphael/Dropbox/Studium/PhD/data/sampleData/";
    std::string path1 = "/Users/Raphael/Dropbox/Studium/PhD/data/sampleData/";

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
//    double medianNoise = pcaKNN(Dt, t_points);
    double medianNoise = 0;
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
    exportColoredFacetsPLY(Dt, ofn, 1);


    exportCellCenter(ofn, Dt);

    std::vector<std::vector<Facet_score>> problematic_facets_per_edge;
    nonManifoldEdges(Dt, problematic_facets_per_edge);
    selectFacets(Dt, problematic_facets_per_edge, ofn);



    // create surface mesh
    Polyhedron out_mesh;

    bool oriented = CGAL::Polygon_mesh_processing::orient_polygon_soup(remaining_points, remaining_facets);
    std::cout << "Created duplicate vertices to generated polygon surface? " << oriented << std::endl;
    std::cout << "... 0 meaning additional vertices were added for orientation (and probably for ensuring manifoldness)." << std::endl;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(remaining_points, remaining_facets, out_mesh);
//    int erased_components = out_mesh.keep_largest_connected_components(1);
//    std::cout << "Number of erased components from surface mesh: " << erased_components << std::endl;

//    std::vector<std::vector<int>> polygons;
//    std::vector<Point> points;
//    // create surface mesh, with orient and nb_of_components_to_keep
//    createSurfaceMesh(Dt, points, polygons, out_mesh, 1, 1);

    // export surface mesh as OFF
    exportOFF(out_mesh, ofn);

//    // Quality control
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




