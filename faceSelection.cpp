#include <cgal_typedefs.h>
#include <fileIO.h>
#include <rayTracing.h>
#include <tetIntersection.h>

typedef std::pair<std::pair<Cell_handle, int>, double> Facet_score;

// get non manifold edges
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


void exportProblematicFacets(Delaunay& Dt,
                  std::vector<std::vector<Facet_score>>& problematic_facets_per_edge,
                  std::string path){

    std::fstream fo;
    fo.open(path+"_facets_center.ply", std::fstream::out);
    int nf = problematic_facets_per_edge.size()*4;

    printPLYHeader(fo,
                   nf, 0,
                   0, true, false, false, true,
                   false,
                   15);

    for(int i = 0; i < problematic_facets_per_edge.size(); i++){

        // get facets per edge and sort them according to their weight
        std::vector<Facet_score> problematic_facets = problematic_facets_per_edge[i];

        for(int j = 0; j < problematic_facets.size(); j++){

            Facet fac = problematic_facets[j].first;
            if(Dt.is_infinite(fac))
                continue;
            Triangle tri = Dt.triangle(fac);
            Point p = CGAL::centroid(tri.vertex(0),tri.vertex(1),tri.vertex(2));
            fo  << p << " " << problematic_facets[j].second << std::endl;
        }

    }

    fo.close();
}


void exportProblematicFacets2(Delaunay& Dt,
                  std::vector<std::vector<Facet_score>>& problematic_facets_per_edge,
                  std::string path){


    std::unordered_set<Vertex_handle> points;
    std::vector<Facet_score> remaining_facets;

    for(int i = 0; i < problematic_facets_per_edge.size(); i++){

        // get facets per edge and sort them according to their weight
        std::vector<Facet_score> problematic_facets = problematic_facets_per_edge[i];
        std::sort(problematic_facets.begin(), problematic_facets.end(), sortbysec);

        double min_score = problematic_facets[0].second;
        double max_score = problematic_facets.back().second;

        for(int j = 0; j < problematic_facets.size(); j++){

            double scaled_score = 255*((problematic_facets[j].second-min_score)/(max_score-min_score));

            // save the first two facets with the lowest score:
            Cell_handle c1 = problematic_facets[j].first.first;
            if(Dt.is_infinite(c1))
                continue;
            int vidx = problematic_facets[j].first.second;
            for(int j = vidx + 1; j <= vidx + 3; j++){
                // so c->vertex() gives me the global vertex handle from the Dt
                points.insert(c1->vertex(j%4));
            }
            remaining_facets.push_back(std::make_pair(problematic_facets[j].first, scaled_score));


        }
    }

    std::fstream fo;
    fo.open(path+"_problematic_facets.ply", std::fstream::out);
    int nv = points.size();
    int nf = remaining_facets.size();

    printPLYHeader(fo,
                   nv, nf,
                   true, true, false, false, false,
                   true,
                   15);

    std::unordered_set<Vertex_handle>::iterator vit;
    int vidx = 0;
    for(vit = points.begin(); vit != points.end(); vit++){
        // reset the vertex index here, because I need to know the order of exactly this loop here
        // for the indexing of the facets in the PLY file
        (*vit)->info().idx = vidx++;
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

    std::vector<Facet_score>::iterator fit;
    for(fit = remaining_facets.begin(); fit != remaining_facets.end(); fit++){
        // start printed facet line with a 3
        fo << 3 << ' ';
        Cell_handle c = fit->first.first;
        int vidx = fit->first.second;
        for(int j = vidx + 1; j <= vidx + 3; j++){
            // so c->vertex() gives me the global vertex handle from the Dt
            fo << c->vertex(j%4)->info().idx << ' ';
        }
        // put a color on the face, so that in Meshlab I can activate the color by facet mode, to compare with the "colored facet file"
        fo << "0 " << int(fit->second) << " 0";
        fo << std::endl;
    }
    fo.close();
}

void exportSelectedFacets(Delaunay& Dt,
                  std::vector<std::vector<Facet_score>>& problematic_facets_per_edge,
                  std::string path){


    std::unordered_set<Vertex_handle> points;
    std::vector<Facet> remaining_facets;

    for(int i = 0; i < problematic_facets_per_edge.size(); i++){

        // get facets per edge and sort them according to their weight
        std::vector<Facet_score> problematic_facets = problematic_facets_per_edge[i];
        std::sort(problematic_facets.begin(), problematic_facets.end(), sortbysec);

        // save the first two facets with the lowest score:
        Cell_handle c1 = problematic_facets[0].first.first;
        if(Dt.is_infinite(c1))
            continue;
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

    }

    std::fstream fo;
    fo.open(path+"_selected_facets.ply", std::fstream::out);
    int nv = points.size();
    int nf = remaining_facets.size();

    printPLYHeader(fo,
                   nv, nf,
                   true, true, false, false, false,
                   true,
                   15);

    std::unordered_set<Vertex_handle>::iterator vit;
    int vidx = 0;
    for(vit = points.begin(); vit != points.end(); vit++){
        // reset the vertex index here, because I need to know the order of exactly this loop here
        // for the indexing of the facets in the PLY file
        (*vit)->info().idx = vidx++;
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
        // put a color on the face, so that in Meshlab I can activate the color by facet mode, to compare with the "colored facet file"
        fo << "0 200 0";
        fo << std::endl;
    }
    fo.close();
}

