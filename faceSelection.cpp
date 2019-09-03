#include <cgal_typedefs.h>
#include <fileIO.h>
#include <rayTracing.h>
#include <tetIntersection.h>

typedef std::pair<std::pair<Cell_handle, int>, double> Facet_score;

// TODO: ask Laurent how to use QUBO solver in order to optimize clique with manifold constrained.


// idea: take a clique, starting with a non manifold edge, ALL its sourounding facets and their facet neighbours.
// now optimize a binary linear program by enforcing that an edge of this clique either has 2 or 0 facet neighbours.
// the idea is that the neighbouring facets of the neighbours of the non manifold edge, which a lot of them did not get selected, now
// contribute with their data attachment (low difference between inside and outside label -> take the inverse)
// to the non selection of the obvious case where the
// 2 facet around the non-manifold edge with the highest difference of inside outside get selected.



// get non manifold edges
int isManifoldEdge(Delaunay& Dt, Delaunay::Finite_edges_iterator& e){

    Delaunay::Cell_circulator cc = Dt.incident_cells(*e);
    Cell_handle first_cell = cc;
    int borders = 0;
    do{
        // first cell
        int cc_label;
        if(cc->info().manifold_label < 2)
            cc_label = cc->info().manifold_label;
        else
            cc_label = cc->info().gc_label;
        // go to next cell
        int nc_label;
        cc++;
        if(cc->info().manifold_label < 2)
            nc_label = cc->info().manifold_label;
        else
            nc_label = cc->info().gc_label;
        // check if there is a border
        if(cc_label != nc_label){borders+=1;}
    }
    while(first_cell != cc && borders < 3);

    if(borders > 2)
        return 0;
    else
        return 1;
}

// get non manifold edges
int isManifoldEdge(Delaunay& Dt, Edge e){

    Delaunay::Cell_circulator cc = Dt.incident_cells(e);
    Cell_handle first_cell = cc;
    int borders = 0;
    do{
        // first cell
        int cc_label;
        if(cc->info().manifold_label < 2)
            cc_label = cc->info().manifold_label;
        else
            cc_label = cc->info().gc_label;
        // go to next cell
        int nc_label;
        cc++;
        if(cc->info().manifold_label < 2)
            nc_label = cc->info().manifold_label;
        else
            nc_label = cc->info().gc_label;
        // check if there is a border
        if(cc_label != nc_label){borders+=1;}
    }
    while(first_cell != cc && borders < 3);

    if(borders > 2)
        return 0;
    else
        return 1;
}


int isManifoldClique(Delaunay& Dt, Delaunay::Finite_edges_iterator& e){

    if(!isManifoldEdge(Dt,e))
        return 0;

    Delaunay::Cell_circulator cc = Dt.incident_cells(*e, e->first);
    Cell_handle first_cell = cc;
    Edge current_edge = *e;
    do{

        int i = current_edge.second;
        int j = current_edge.third;
        int k = Dt.next_around_edge(i,j);

        // check first edge
        Edge e1 = CGAL::Triple<Cell_handle,int,int>(cc,i,k);
        if(!isManifoldEdge(Dt,e1))
            return 0;
        // check second edge
        Edge e2 = CGAL::Triple<Cell_handle,int,int>(cc,j,k);
        if(!isManifoldEdge(Dt,e2))
            return 0;
        // check third edge
        Edge e3 = CGAL::Triple<Cell_handle,int,int>(cc,k,6-(i+j+k));
        if(!isManifoldEdge(Dt,e3))
            return 0;

        // go to next cell with its corresponding edge
        cc++;
        int new_i;
        switch(i){
        case 0: new_i = 1;
            break;
        case 1: new_i = 0;
            break;
        case 2: new_i = 2;
            break;
        case 3: new_i = 3;
        }
        int new_j;
        switch(j){
        case 0: new_j = 1;
            break;
        case 1: new_j = 0;
            break;
        case 2: new_j = 2;
            break;
        case 3: new_j = 3;
        }
        current_edge = CGAL::Triple<Cell_handle, int, int>(cc, new_j, new_i);
    }
    while(first_cell != cc);

    return 1;

}



int getCellLabel(Cell_handle& c, bool manifold=true){

    if(c->info().manifold_label < 2 && manifold)
        return c->info().manifold_label;
    else
        return c->info().gc_label;
}


double nonManifoldCliqueEnergy(const Delaunay& Dt, Delaunay::Finite_edges_iterator& e, double reg_weight){

    Delaunay::Cell_circulator cc = Dt.incident_cells(*e, e->first);
    Cell_handle first_cell = cc;
    double unary = 0;
    double binary = 0;
    do{
        // first cell
        Cell_handle current_cell = cc;
        int cc_label = getCellLabel(current_cell);

        // add the unary term of this cell
        if(cc_label==1){unary+=cc->info().inside_score;} // if label = 1 = outside, penalize with inside score
        else{unary+=cc->info().outside_score;} // if label = 0 = inside, penalize with outside score

        // add the binary term of the other two facets that are not connected to this edge
        Facet f1 = std::make_pair(cc,e->second);
        Cell_handle mc1 = Dt.mirror_facet(f1).first;
        int mc1_label = getCellLabel(mc1);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(cc) && Dt.is_infinite(mc1)) && mc1_label != cc_label){
            Triangle tri = Dt.triangle(f1);
            binary+=sqrt(tri.squared_area());
        }
        Facet f2 = std::make_pair(cc,e->third);
        Cell_handle mc2 = Dt.mirror_facet(f2).first;
        int mc2_label = getCellLabel(mc2);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(cc) && Dt.is_infinite(mc2)) && mc2_label != cc_label){
            Triangle tri = Dt.triangle(f2);
            binary+=sqrt(tri.squared_area());
        }

        // go to next cell and add the binary term between previous cell and next cell
        cc++;
        int nidx = current_cell->index(cc); // the index of the new cell seen from the old cell = facet between previous and new cell
        Facet f3 = std::make_pair(current_cell,nidx);
        int mc3_label = getCellLabel(current_cell);
        //if they are not both infinite and have a different label, add their area to the binary term
        if(!(Dt.is_infinite(current_cell) && Dt.is_infinite(cc)) && mc3_label != cc_label){
            Triangle tri = Dt.triangle(f3);
            binary+=sqrt(tri.squared_area());
        }
    }
    while(first_cell != cc);
    return unary+reg_weight*binary;
}

typedef std::pair<double, boost::dynamic_bitset<>> Combination_score;
int makeCombinations(int n, std::vector<Combination_score>& c){
    const int x = pow(2,n);
    for(int i = 0; i < x; i++){
        const boost::dynamic_bitset<> b(n, i);
        c.push_back(std::make_pair(0.0,b));
    }
    return x;
}

bool sortCombinations(const Combination_score &a,
              const Combination_score &b)
{
    return (a.first < b.first);
}

void fixNonManifoldEdges(Delaunay& Dt, double regularization_weight){

    // iterate over all edges of the Dt
    Delaunay::Finite_edges_iterator fei;
    for(fei = Dt.finite_edges_begin(); fei != Dt.finite_edges_end(); fei++){

        // if the current edge is (already) manifold, continue
        if(isManifoldEdge(Dt, fei))
            continue;

        // make a container with all incident cells
        std::vector<Cell_handle> cells_around_nmedge;
        Delaunay::Cell_circulator cc = Dt.incident_cells(*fei);
        Cell_handle first_cell = cc;
        do{cells_around_nmedge.push_back(cc++);}
        while(first_cell != cc);

        // make the bitset with according size
        std::vector<Combination_score> combinations;
        int number_of_cells = cells_around_nmedge.size();
        int number_of_possible_combinations = makeCombinations(number_of_cells, combinations);

        // iterate over the container, relabel the cells with the bitset, and get their corresponding energy and save it in a vector
        for(int c = 0; c < number_of_possible_combinations; c++){
            for(int v = 0; v < number_of_cells; v++)
                cells_around_nmedge[v]->info().manifold_label = combinations[c].second[v];
            // calculate the energy of the current configuration
            combinations[c].first = nonManifoldCliqueEnergy(Dt, fei, regularization_weight);
        }

        // sort the container while keeping track of its original index, which gives you the corresponding configuration from the bitset index
        std::sort(combinations.begin(), combinations.end(), sortCombinations);

        // take the lowest energy and check if it is manifold
        int solution_found = 0;
        for(int sc = 0; sc < number_of_possible_combinations; sc++){

            // relabel to the correct combination
            for(int v = 0; v < number_of_cells; v++){
                cells_around_nmedge[v]->info().manifold_label = combinations[sc].second[v];
            }
            if(isManifoldClique(Dt, fei)){
                for(int v = 0; v < number_of_cells; v++){
                    cells_around_nmedge[v]->info().gc_label = combinations[sc].second[v];
                }
                solution_found = 1;
                break;
            }
        } // end of trying all possible solutions
        if(!solution_found){
            std::cout << "no manifold combination found for edge ("
                      << fei->first->info().idx << ", "
                      << fei->second << ", "
                      << fei->third << ")"
                      << std::endl;
        }
    } // end of iteration over all edges of the Delaunay
}






////////////////////////////////////
//////// Point manifoldness ////////
////////////////////////////////////


// TODO: could also implement this with kind of a recursive function that always looks in neighbourhood
// so start with a seed cell, go to next cell and see if it has a different label and if it is connected to the previous cell
// ... if I have more than two regions, point is not manifold

// get non manifold points
int isManifoldPoint(Delaunay& Dt, Vertex_handle p){

    std::vector<Cell_handle> incident_cells;
    Dt.incident_cells(p, std::back_inserter(incident_cells));
    std::vector<Cell_handle>::iterator cc;
    // sort cells in inside and outside cells
    std::vector<Cell_handle> inside_cells;
    std::vector<Cell_handle> outside_cells;
    for(cc = incident_cells.begin(); cc != incident_cells.end(); cc++){
        if(getCellLabel(*cc) == 0)
            inside_cells.push_back(*cc);
        else
            outside_cells.push_back(*cc);
    }
    // if one of them is empty, point is manifold
    int is = inside_cells.size();
    int os = outside_cells.size();
    std::vector<Cell_handle> region;
    std::vector<Cell_handle> rest;
    std::vector<Cell_handle> stack;
    if(is < 2 || os < 2){
        return 1;
    }
    else if(is <= os){
        region.push_back(inside_cells[0]);
        for(int i = 1; i < is; i++)
            rest.push_back(inside_cells[i]);
    }
    else if(is > os){
        region.push_back(outside_cells[0]);
        for(int i = 1; i < os; i++)
            rest.push_back(outside_cells[i]);
    }
    int new_stack = 1000;
    int old_stack = 1001;
    do{
        Cell_handle first_rest_cell = rest[0];
        bool first_rest_cell_added = 0;
        for(auto r = region.begin(); r != region.end(); r++){
            if(first_rest_cell->has_neighbor(*r)){
                region.push_back(first_rest_cell);
                rest.erase(rest.begin());
                first_rest_cell_added = 1;
                break;
            }
        }
        if(!first_rest_cell_added){
            stack.push_back(first_rest_cell);
            rest.erase(rest.begin());
        }
        if(rest.size() == 0){
            if(stack.size()==0)
                return 1;
            rest = stack;
            old_stack=new_stack;
            new_stack=stack.size();
            stack.clear();
        }
    }while(old_stack > new_stack);
    if(rest.size()==0 && stack.size()==0)
        return 1;
    else
        return 0;
}

// get non manifold points
int isManifoldPoint(Delaunay& Dt, Delaunay::Finite_vertices_iterator& p, std::vector<Cell_handle>& incident_cells){

    Dt.incident_cells(p, std::back_inserter(incident_cells));
    std::vector<Cell_handle>::iterator cc;
    // sort cells in inside and outside cells
    std::vector<Cell_handle> inside_cells;
    std::vector<Cell_handle> outside_cells;
    for(cc = incident_cells.begin(); cc != incident_cells.end(); cc++){
        if(getCellLabel(*cc) == 0)
            inside_cells.push_back(*cc);
        else
            outside_cells.push_back(*cc);
    }
    // if one of them is empty, point is manifold
    int is = inside_cells.size();
    int os = outside_cells.size();
    std::vector<Cell_handle> region;
    std::vector<Cell_handle> rest;
    std::vector<Cell_handle> stack;
    if(is < 2 || os < 2){
        return 1;
    }
    else if(is <= os){
        region.push_back(inside_cells[0]);
        for(int i = 1; i < is; i++)
            rest.push_back(inside_cells[i]);
    }
    else if(is > os){
        region.push_back(outside_cells[0]);
        for(int i = 1; i < os; i++)
            rest.push_back(outside_cells[i]);
    }
    int new_stack = 1000;
    int old_stack = 1001;
    do{
        Cell_handle first_rest_cell = rest[0];
        bool first_rest_cell_added = 0;
        for(auto r = region.begin(); r != region.end(); r++){
            if(first_rest_cell->has_neighbor(*r)){
                region.push_back(first_rest_cell);
                rest.erase(rest.begin());
                first_rest_cell_added = 1;
                break;
            }
        }
        if(!first_rest_cell_added){
            stack.push_back(first_rest_cell);
            rest.erase(rest.begin());
        }
        if(rest.size() == 0){
            if(stack.size()==0)
                return 1;
            rest = stack;
            old_stack=new_stack;
            new_stack=stack.size();
            stack.clear();
        }
    }while(old_stack > new_stack);
    if(rest.size()==0 && stack.size()==0)
        return 1;
    else
        return 0;
}



// get non manifold edges
int isPointManifoldClique2(Delaunay& Dt, Delaunay::Finite_vertices_iterator& p){


    // first check if input point is manifold
    if(!isManifoldPoint(Dt, p))
        return 0;

    // now assemble all other points in a set
    std::vector<Cell_handle> incident_cells;
    Dt.incident_cells(p, std::back_inserter(incident_cells));
    std::vector<Cell_handle>::iterator cc;
    std::unordered_set<Vertex_handle> vertices;
    for(cc = incident_cells.begin(); cc != incident_cells.end(); cc++){

        // add its 4th face neighbour
        int idx = (*cc)->index(p);

        // switch case
        // if idx 0, other vertices to add are 1,2,3, and so on
        std::vector<int> other_idx(3);
        switch(idx){
        case 0: other_idx[0] = 1; other_idx[1] = 2; other_idx[2] = 3;
            break;
        case 1: other_idx[0] = 0; other_idx[1] = 2; other_idx[2] = 3;
            break;
        case 2: other_idx[0] = 0; other_idx[1] = 1; other_idx[2] = 3;
            break;
        case 3: other_idx[0] = 0; other_idx[1] = 1; other_idx[2] = 2;
        }
        for(int i = 0; i < 3; i++)
            vertices.insert((*cc)->vertex(other_idx[i]));
    }
    // now check all of them
    std::unordered_set<Vertex_handle>::iterator vs;
    for (vs = vertices.begin(); vs != vertices.end(); vs++){
        if(!isManifoldPoint(Dt, *vs))
            return 0;
    }
    return 1;
}

// get non manifold edges
int isPointManifoldClique(Delaunay& Dt, Delaunay::Finite_vertices_iterator& p){

    std::vector<Cell_handle> incident_cells;
    Dt.incident_cells(p, std::back_inserter(incident_cells));
    // TODO: also add the cells around the points around it,
    // rest should work exactly the same

    std::vector<Cell_handle>::iterator cc;
    // sort cells in inside and outside cells
    std::unordered_set<Cell_handle> inside_cells;
    std::unordered_set<Cell_handle> outside_cells;
    for(cc = incident_cells.begin(); cc != incident_cells.end(); cc++){

        // add the current cell
        if(getCellLabel(*cc, false) == 0)
            inside_cells.insert(*cc);
        else
            outside_cells.insert(*cc);

        // add its 4th face neighbour
        Cell_handle nei = (*cc)->neighbor((*cc)->index(p));
        if(getCellLabel(nei, false) == 0)
            inside_cells.insert(nei);
        else
            outside_cells.insert(nei);

    }


    // if one of them is empty, point is manifold
    int is = inside_cells.size();
    int os = outside_cells.size();
    std::vector<Cell_handle> region;
    std::vector<Cell_handle> rest;
    std::vector<Cell_handle> stack;
    if(is < 2 || os < 2){
        return 1;
    }
    else if(is <= os){
        std::unordered_set<Cell_handle>::iterator cit = inside_cells.begin();
        region.push_back(*cit);
        while(cit != inside_cells.end()){
            cit++;
            rest.push_back(*cit);
        }
    }
    else if(is > os){
        std::unordered_set<Cell_handle>::iterator cit = outside_cells.begin();
        region.push_back(*cit);
        while(cit != outside_cells.end()){
            cit++;
            rest.push_back(*cit);
        }
    }
    int new_stack = 1000;
    int old_stack = 1001;
    do{
        Cell_handle first_rest_cell = rest[0];
        bool first_rest_cell_added = 0;
        for(auto r = region.begin(); r != region.end(); r++){
            if(first_rest_cell->has_neighbor(*r)){
                region.push_back(first_rest_cell);
                rest.erase(rest.begin());
                first_rest_cell_added = 1;
                break;
            }
        }
        if(!first_rest_cell_added){
            stack.push_back(first_rest_cell);
            rest.erase(rest.begin());
        }
        if(rest.size() == 0){
            if(stack.size()==0)
                return 1;
            rest = stack;
            old_stack=new_stack;
            new_stack=stack.size();
            stack.clear();
        }
    }while(old_stack > new_stack);
    if(rest.size()==0 && stack.size()==0)
        return 1;
    else
        return 0;
}


double nonManifoldCliqueEnergy(const Delaunay& Dt, std::vector<Cell_handle> all_cells, double reg_weight){

    std::vector<Cell_handle>::iterator cit;
    double unary = 0;
    double binary = 0;
    for(cit = all_cells.begin(); cit != all_cells.end(); cit++){

        // first cell
        int cc_label = getCellLabel(*cit);

        // add the unary term of this cell
        if(cc_label==1){unary+=(*cit)->info().inside_score;} // if label = 1 = outside, penalize with inside score
        else{unary+=(*cit)->info().outside_score;} // if label = 0 = inside, penalize with outside score

        int current_index = (*cit)->info().idx;
        for(int i = 0; i < 4; i++){

            Cell_handle neighbour_cell = (*cit)->neighbor(i);
            int neighbour_index = neighbour_cell->info().idx;

            // if current_cell and neighbour_cell are BOTH infinite, then continue
            if(Dt.is_infinite(*cit) && Dt.is_infinite(neighbour_cell)){
                continue;
            }

            // prevent to call setNeighbour(s2,s1) if setNeighbour(s1,s2)was already called
            if(neighbour_index < current_index)
                continue;

            // since i is giving me the cell that is opposite of vertex i, as well as the facet that is opposite of vertex i, I can just use that same index
            Triangle tri = Dt.triangle(*cit, i);
            binary+=sqrt(tri.squared_area());

        }

    }
    return unary+reg_weight*binary;
}


void getNonManifoldPoints(Delaunay& Dt, double regularization_weight, std::string path){

    auto start = std::chrono::high_resolution_clock::now();

    path+="_nmPoints.ply";

    std::fstream fo;
    fo.open(path, std::fstream::out);

    int nv = Dt.number_of_vertices();

    printPLYHeader(fo,
                   nv, 0,
                   false, true, false, false, false, false, 15);

    // iterate over all edges of the Dt
    Delaunay::Finite_vertices_iterator fvi;
    int count = 0;
    for(fvi = Dt.finite_vertices_begin(); fvi != Dt.finite_vertices_end(); fvi++){

        if(!isManifoldPoint(Dt, fvi)){
            fo << fvi->point() << " " << "0 255 0" << std::endl;
            std::array<unsigned char, 3> col = {0,255,0};
            fvi->info().color = col;
            count++;
        }
        else{
            fo << fvi->point() << " " << "0 0 0" << std::endl;
        }
    }
    std::cout << "number of non-manifold points: " << count << std::endl;
    fo.close();

}


void fixNonManifoldPoints(Delaunay& Dt, double regularization_weight){

    auto start = std::chrono::high_resolution_clock::now();

    // iterate over all edges of the Dt
    Delaunay::Finite_vertices_iterator fvi;
    for(fvi = Dt.finite_vertices_begin(); fvi != Dt.finite_vertices_end(); fvi++){

        std::vector<Cell_handle> incident_cells;
        if(!isManifoldPoint(Dt, fvi, incident_cells))
            continue;


        // make the bitset with according size
        std::vector<Combination_score> combinations;
        int number_of_cells = incident_cells.size();
        int number_of_possible_combinations = makeCombinations(number_of_cells, combinations);

        // iterate over the container, relabel the cells with the bitset, and get their corresponding energy and save it in a vector
        for(int c = 0; c < number_of_possible_combinations; c++){
            for(int v = 0; v < number_of_cells; v++)
                incident_cells[v]->info().manifold_label = combinations[c].second[v];
            // calculate the energy of the current configuration
            combinations[c].first = nonManifoldCliqueEnergy(Dt, incident_cells, regularization_weight);
        }

        // sort the container while keeping track of its original index, which gives you the corresponding configuration from the bitset index
        std::sort(combinations.begin(), combinations.end(), sortCombinations);

        // take the lowest energy and check if it is manifold
        int solution_found = 0;
        for(int sc = 0; sc < number_of_possible_combinations; sc++){

            // relabel to the correct combination
            for(int v = 0; v < number_of_cells; v++){
                incident_cells[v]->info().manifold_label = combinations[sc].second[v];
            }
            if(isPointManifoldClique(Dt, fvi)){
                for(int v = 0; v < number_of_cells; v++){
                    incident_cells[v]->info().gc_label = combinations[sc].second[v];
                }
                solution_found = 1;
                break;
            }
        } // end of trying all possible solutions
        if(!solution_found){
            std::cout << "no manifold combination found for edge ("
                      << fvi->info().idx
                      << std::endl;
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Fixed non-manifold points in " << duration.count() << "s" << std::endl;
}



















































////////////////////////////////////////////////////////////////////
/////////////////////// old facet based stuff //////////////////////
////////////////////////////////////////////////////////////////////

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
            if(c1->info().gc_label != c2->info().gc_label){
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


double facetScore(Cell_handle& c1, Cell_handle& c2){
    return sqrt(pow(c1->info().inside_score - c2->info().inside_score,2) +
                pow(c1->info().outside_score - c2->info().outside_score,2));
}
bool sortFacets(const Facet_score &a,
              const Facet_score &b)
{
    return (a.second < b.second);
}

// get cliques as a list of facets including their score
void getNonManifoldCliques(Delaunay& Dt, std::vector<std::vector<Facet_score>>& problematic_facets_per_edge){

    // iterate over all edges of the Delaunay to find the non-manifold ones
    Delaunay::Finite_edges_iterator fei;
    int edge_count = 0;
    for(fei = Dt.finite_edges_begin(); fei != Dt.finite_edges_end(); fei++){
        std::cout << edge_count++ << std::endl;
        if(Dt.is_infinite(*fei)){
            std::cout << "...is infinite edge" << std::endl;
            continue;
        }
        // circulate over all the facet of the current edge
        Delaunay::Facet_circulator fac = Dt.incident_facets(*fei);
        Facet nm_f1 = *fac; // save the first facet to now when to stop again
        std::vector<Facet_score> problematic_facets;         // init a vector of facets for the current edge
        do{
            // get the neighbouring of the current facet fac
            Cell_handle c1 = fac->first;
            Cell_handle c2 = Dt.mirror_facet(*fac).first;
            // first add the facets around the non-manifold edge which have different labels (usually 4)
            if(c1->info().gc_label != c2->info().gc_label){
                double score = facetScore(c1,c2);
                problematic_facets.push_back(std::make_pair(std::make_pair(c1,fac->second),score));
            } // end of facet around non manifold edge selection if they have different labels
            fac++;
        // end of facet around current non manifold edge circulation
        }while(*fac != nm_f1);
        // check if we are even talking about a non-manifold edge, and if not, continue to the next edge
        if(problematic_facets.size()<3)
            continue;
        // now add the neighbouring facets of the problematic facets
        std::vector<Facet_score> new_problematic_facets;
        for(int j = 0; j < problematic_facets.size(); j++){
            Facet current_facet = problematic_facets[j].first;
//            std::cout << "new facet" << std::endl;
            // get the vertex index of the current facet, and the current edge
            int facet_vertex = current_facet.second;
            int e1_vertex = fei->second;
            int e2_vertex = fei->third;
            if(facet_vertex == e1_vertex || facet_vertex == e2_vertex){
                current_facet = Dt.mirror_facet(current_facet);
                facet_vertex = current_facet.second;
                if(facet_vertex == e1_vertex || facet_vertex == e2_vertex){
//                    std::cout << "not working after switch" << std::endl;
                    continue;
                }

            }
            // now build the two other edges of the current facet
            for(int i = e1_vertex + 1; i <= e1_vertex + 3; i++){
                int vertex = i%4;
                // iterate until I am at the 4th vertex of this cell, which is the first vertex of the two remaining edges of the current facet
                if(vertex==facet_vertex || vertex==e1_vertex || vertex==e2_vertex)
                    continue;
                // first edge (edge is a triple of Cell, idx, idx)
                Edge en1 = CGAL::Triple<Cell_handle, int, int>(current_facet.first, vertex, e1_vertex);
                if(Dt.is_infinite(en1)){
                    std::cout << "infinite edge" << std::endl;
                    continue;
                }
                // get all the incident facets to this edge, starting at the current facet
                Delaunay::Facet_circulator fac1 = Dt.incident_facets(en1, current_facet);
                fac1++;
//                std::cout << "current facet (" << current_facet.first->info().idx << ", " << current_facet.second << ")" << std::endl;
                Facet mirror_facet = Dt.mirror_facet(current_facet);
//                std::cout << "mirror facet (" << mirror_facet.first->info().idx << ", " << mirror_facet.second << ")" << std::endl;
                while(*fac1 != current_facet){
                    Cell_handle c1 = fac1->first;
                    Cell_handle c2 = Dt.mirror_facet(*fac1).first;
                    double score = facetScore(c1,c2);
                    new_problematic_facets.push_back(std::make_pair(std::make_pair(c1,fac1->second),score));
//                    std::cout << "before facet (" << fac1->first->info().idx << ", " << fac1->second << ")" << std::endl;
                    fac1++;
                    if(Dt.mirror_facet(*fac1)==current_facet)
                        break;
//                    std::cout << "after facet (" << fac1->first->info().idx << ", " << fac1->second << ")" << std::endl;
                }
                // second edge
                Edge en2 = CGAL::Triple<Cell_handle, int, int>(current_facet.first, e2_vertex, vertex);
                if(Dt.is_infinite(en2)){
                    std::cout << "infinite edge" << std::endl;
                    continue;
                }
                // get all the incident facets to this edge, starting at the current facet
                Delaunay::Facet_circulator fac2 = Dt.incident_facets(en2, current_facet);
                fac2++;
                while(*fac2 != current_facet){
                    Cell_handle c1 = fac2->first;
                    Cell_handle c2 = Dt.mirror_facet(*fac2).first;
                    double score = facetScore(c1,c2);
                    new_problematic_facets.push_back(std::make_pair(std::make_pair(c1,fac2->second),score));
                    fac2++;
                    if(Dt.mirror_facet(*fac2)==current_facet)
                        break;
                }
            }
        }
        // add the new problematic facets to the original facets around the originial non-manifold edge
        problematic_facets.insert(problematic_facets.end(), new_problematic_facets.begin(),
                           new_problematic_facets.end());
        problematic_facets_per_edge.push_back(problematic_facets);
    }
    std::cout << "number of non manifold edges: " << problematic_facets_per_edge.size() << std::endl;
}




void exportProblematicFacets(Delaunay& Dt,
                  std::vector<std::vector<Facet_score>>& problematic_facets_per_edge,
                  std::string path){


    std::unordered_set<Vertex_handle> points;
    std::vector<Facet_score> remaining_facets;

    for(int i = 0; i < problematic_facets_per_edge.size(); i++){

        // get facets per edge and sort them according to their weight
        std::vector<Facet_score> problematic_facets = problematic_facets_per_edge[i];
        std::sort(problematic_facets.begin(), problematic_facets.end(), sortFacets);

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
        std::sort(problematic_facets.begin(), problematic_facets.end(), sortFacets);

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

