#include <cgal_typedefs.h>

// for GCoptimization
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "GCoptimization.h"


////////////////////////////////////////////////////////////
////////////////////// Feature scaling /////////////////////
////////////////////////////////////////////////////////////
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
        cit->info().inside_score = (cit->info().inside_score - inside_mean) / inside_stdev;
        cit->info().outside_score = (cit->info().outside_score - outside_mean) / outside_stdev;
//        cit->info().inside_score;
//        cit->info().outside_score;
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

void logScore(Delaunay& Dt){

    Delaunay::Finite_cells_iterator cit;
    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){

        double inside_score = cit->info().inside_score;
        double outside_score = cit->info().outside_score;
        cit->info().inside_score = log(inside_score+0.00000001) + log(1 - outside_score);
        cit->info().outside_score = log(1-inside_score) + log(outside_score+0.00000001);
        std::cout << "inside score before: " << inside_score <<
                     "  inside score after: " << cit->info().inside_score << std::endl;
        std::cout << "outside score before: " << outside_score <<
                     "  outside score after: " << cit->info().outside_score << std::endl;
    }

}

void softmax(Delaunay& Dt){

    Delaunay::Finite_cells_iterator cit;
    for(cit = Dt.finite_cells_begin(); cit != Dt.finite_cells_end(); cit++){

        double inside_score = cit->info().inside_score;
        double outside_score = cit->info().outside_score;
        cit->info().inside_score = log(exp(inside_score) / (exp(inside_score) + exp(outside_score)));
        cit->info().outside_score = log(1 - (exp(inside_score) / (exp(inside_score) + exp(outside_score))));
//        std::cout << "inside score before: " << inside_score <<
//                     "  inside score after: " << cit->info().inside_score << std::endl;
//        std::cout << "outside score before: " << outside_score <<
//                     "  outside score after: " << cit->info().outside_score << std::endl;
    }
}






////////////////////////////////////////////////////////////
//////////////////////// Optimization //////////////////////
////////////////////////////////////////////////////////////
//// in this version, set data and smoothness terms using arrays
//// grid neighborhood is set up "manually". Uses spatially varying terms. Namely
//// V(p1,p2,l1,l2) = w_{p1,p2}*[min((l1-l2)*(l1-l2),4)], with
//// w_{p1,p2} = p1+p2 if |p1-p2| == 1 and w_{p1,p2} = p1*p2 if |p1-p2| is not 1
//std::pair<std::map<Cell_handle, int>, std::vector<int>> GeneralGraph_DArraySArraySpatVarying(std::pair<Delaunay&, Cell_map&> dt_cells, std::map<Cell_handle, int>& cell_indexMap, std::vector<int> result, int num_iterations)
void GeneralGraph_DArraySArraySpatVarying(Delaunay& Dt, float area_weight, int num_iterations)
{
    std::cout << "Starting Optimization..." << std::endl;

    int num_cells = Dt.number_of_cells();
    int num_labels = 2;
    GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(num_cells,num_labels);

    // first set up the array for data costs and set the initial label in the same loop
    float *data = new float[num_cells*num_labels];
    int idx = 0;
    Delaunay::All_cells_iterator cft;
    // iterate over the all_cells map
    for(cft = Dt.all_cells_begin(); cft!=Dt.all_cells_end(); cft++)
    {
        // set an index for each cell
        cft->info().idx = idx;
        // I am initializing my label s.t. if the outside vote is bigger than the inside vote, than it should be labelled 1 (outside) - and vice versa (0 = inside)
        // that means the cost/penalty for labelling a cell that initially has label 1 with the opposite label, is the opposite vote (so the inside vote)
        // so labelling 0 costs outside votes
        data[idx*2+0] = cft->info().outside_score;
        // and labelling 1 costs inside votes
        data[idx*2+1] = cft->info().inside_score;

        // set an initial label
        if(cft->info().outside_score > cft->info().inside_score)
            {gc->setLabel(idx, 1);}
        else
            {gc->setLabel(idx, 0);}
        // increase index for next cell
        idx++;
    }

    // next set up the array for smooth costs
    float *smooth = new float[num_labels*num_labels];
    for ( int l1 = 0; l1 < num_labels; l1++ )
        for (int l2 = 0; l2 < num_labels; l2++ )
            if(l1 == l2){smooth[l1+l2*num_labels] = 0.0;}
            else{smooth[l1+l2*num_labels] = 1.0;}

    try{
        gc->setDataCost(data);
        gc->setSmoothCost(smooth);

        // set neighborhood:
        int current_index;
        int neighbour_index;
        // iterate over the all_cells map
        for(cft = Dt.all_cells_begin(); cft!=Dt.all_cells_end(); cft++)
        {
            Cell_handle current_cell = cft;
            current_index = current_cell->info().idx;
            for(int i = 0; i < 4; i++){

                Cell_handle neighbour_cell = current_cell->neighbor(i);
                neighbour_index = neighbour_cell->info().idx;

                // if current_cell and neighbour_cell are BOTH infinite, then continue
                if(Dt.is_infinite(current_cell) && Dt.is_infinite(neighbour_cell)){
                    continue;
                }

                // prevent to call setNeighbour(s2,s1) if setNeighbour(s1,s2)was already called
                if(neighbour_index < current_index)
                    continue;

                // since i is giving me the cell that is opposite of vertex i, as well as the facet that is opposite of vertex i, I can just use that same index
                Triangle tri = Dt.triangle(current_cell, i);
                // TODO: think about replacing the area weight with the beta-skeletion weight from
                // Robust and efficient surface reconstruction from range data
                float area = sqrt(tri.squared_area());

                // call the neighbourhood function
                gc->setNeighbors(current_index, neighbour_index, area_weight*area);

                // TODO: what would be nice if I could make like a second grade neighborhood function, in order
                // to penalize facets (=different labels between neighboring cells) even higher if facets from the same cell
                // also already have different labels
            }
        }

        // Optimization
        std::cout << "Before optimization data energy is " << gc->giveDataEnergy() << std::endl;
//        for(int i = 0; i < num_cells; i++){
//            std::cout << "cell " << i << " label 0: " << data[i*2+0] << std::endl;
//            std::cout << "cell " << i << " label 1: " << data[i*2+1] << std::endl;
//        }
        std::cout << "Before optimization smoothness energy is " << gc->giveSmoothEnergy() << std::endl;

        std::cout << "Before optimization energy is " << gc->compute_energy() << std::endl;
        // use swap because it is good when you have two labels
        gc->swap(num_iterations);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
//        gc->expansion(num_iterations);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
        std::cout << "After optimization energy is " << gc->compute_energy() << std::endl;

        // save the results in the all_cells Cell_map
        int idx3;
        // iterate over the all_cells map
        for(cft = Dt.all_cells_begin(); cft!=Dt.all_cells_end(); cft++)
        {
            idx3 = cft->info().idx;
//            std::get<3>(it3->second) = gc->whatLabel(idx3);
            cft->info().final_label = gc->whatLabel(idx3);

        }


        delete gc;
    }
    catch (GCException e){
        e.Report();
    }
}


// manually calculate energy terms for checking for correctness
void checkEnergyTerms(const Delaunay& Dt, Cell_map& all_cells, float area_weight)
{

    int num_cells = Dt.number_of_cells();
    int num_labels = 2;

    int* label = new int[num_cells];

    // first set up the array for data costs and set the initial label in the same loop

    int current_idx;
    Cell_map::iterator it;
    // iterate over the all_cells map
    for(it = all_cells.begin(); it!=all_cells.end(); it++)
    {
        Cell_handle current_cell = it->first;

        // data term
        current_idx = std::get<0>(it->second);
//        data = std::get<1>(it->second) + std::get<2>(it->second);

        // smoothness term
        if(std::get<1>(it->second) > std::get<2>(it->second))
            {label[current_idx] = 1;}
        else
            {label[current_idx] = 0;}
    }

    // calculate energy
    float data_energy = 0.0;
    float smoothness_energy = 0.0;
    float total_energy = 0.0;
    // iterate over the all_cells map
    for(it = all_cells.begin(); it!=all_cells.end(); it++)
    {
        Cell_handle current_cell = it->first;
        current_idx = std::get<0>(it->second);

        // data term
        // so I am initializing my label s.t. if the outside vote is bigger than the inside vote, than it should be labelled 1 - and vice versa
        // that means the cost for labelling a cell that has label 1 with label 0, is the opposite vote (so the inside vote)
        float data;
        if(label[current_idx] == 0)
            data = std::get<1>(it->second);
        else
            data = std::get<2>(it->second);
        data_energy+=data;

        for(int i = 0; i < 4; i++){

            Cell_handle neighbour_cell = current_cell->neighbor(i);
            int neighbour_idx = std::get<0>(all_cells.find(neighbour_cell)->second);

            // if current_cell and neighbour_cell are BOTH infinite, then continue
            if(Dt.is_infinite(current_cell) && Dt.is_infinite(neighbour_cell)){
                continue;
            }

            // prevent to call setNeighbour(s2,s1) if setNeighbour(s1,s2)was already called
            if(neighbour_idx < current_idx)
                continue;

            // since i is giving me the cell that is opposite of vertex i, as well as the facet that is opposite of vertex i, I can just use that same index
            Triangle tri = Dt.triangle(current_cell, i);
            float area = sqrt(tri.squared_area());

            // smooth term
            float smooth;
            if(label[current_idx]!=label[neighbour_idx])
                smooth = area_weight*area;
            else
                smooth = 0.0;
            smoothness_energy+=smooth;
        }
    }

    total_energy = data_energy + smoothness_energy;

    std::cout << "Calculated total=data+smoothness energy is: " << total_energy << "=" << data_energy << "+" << smoothness_energy << std::endl;


}
