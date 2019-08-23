#include <cgal_typedefs.h>


////////////////////////////////////////////////////////////
///////////////////////////// PCA //////////////////////////
////////////////////////////////////////////////////////////
// PCA with kNN neighborhood
double pcaKNN(Delaunay& Dt, std::vector<Point>& all_points){

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Start PCA..." << std::endl;

    // calculate the size of the smallest eigenvalue, which sould serve as a good measurement of noise. and use that as the sigma for the score computation
    // problem: this might have a good effect in noise areas, but it also weakens the few important votes in missing data areas

//    // get the point for every vertex
    Delaunay::Finite_vertices_iterator vft;
//    std::vector<Point> all_points;
//    for(vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
//        Point p = vft->point();
//        all_points.push_back(p);
//    }

    std::size_t NN =  CGAL::estimate_global_k_neighbor_scale(all_points);

    // build a kd-tree with all the points
    Tree tree(all_points.begin(), all_points.end());

    std::vector<double> all_eig3s;
    for(vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){

        // set up the neighborhood search function
        Neighbor_search search(tree, vft->point(), NN);

        std::vector<Vertex_handle> av;
        Dt.adjacent_vertices(vft, std::back_inserter(av));

        std::vector<Point> adj_points;
        int nn = 0;
        for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
            adj_points.push_back(it->first);
            nn++;
        }
        Point centroid = CGAL::centroid(adj_points.begin(), adj_points.end(), CGAL::Dimension_tag<0>());
        Eigen::MatrixXd m(3,nn);
        Vector p;
        for(int i = 0; i < nn; i++){
            p = adj_points[i]-centroid;
            m(0,i) = p.x();
            m(1,i) = p.y();
            m(2,i) = p.z();
        }
        if(nn<1)
            std::cout<<"no neighbours" << std::endl;
        Eigen::MatrixXd A(3,3);
        A = (m*(m.transpose()))/nn;

        // compute eigenvalues of the covariance matrix A:
        // implemented from here: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
        double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
        double eig1;
        double eig2;
        double eig3;
        Eigen::MatrixXd I(3,3);
        I = I.setIdentity();
        if (p1 == 0.0){
           eig1 = A(0,0);
           eig2 = A(1,1);
           eig3 = A(2,2);}
        else{
           double q = A.trace()/3.0;
           double p2 = pow((A(0,0) - q),2.0) + pow((A(1,1) - q),2.0) + pow((A(2,2) - q),2.0) + 2.0 * p1;
           double p = sqrt(p2 / 6.0);
           Eigen::MatrixXd B(3,3);
           B = (1.0 / p) * (A - q * I);
           double r = B.determinant() / 2.0;
           double phi;
           if (r <= -1.01)
              phi = M_PI / 3.0;
           else if (r >= 0.99)
              phi = 0.0;
           else
              phi = acos(r) / 3.0;
           // the eigenvalues satisfy eig3 <= eig2 <= eig1
           eig1 = q + 2.0 * p * cos(phi);
           eig3 = q + 2.0 * p * cos(phi + (2.0*M_PI/3.0));
           eig2 = 3.0 * q - eig1 - eig3;
        }

        // compute eigenvector of the third (smallest) eigenvalue:
        Eigen::MatrixXd EV(3,3);
        EV = (A-eig1*I)*(A-eig2*I);
//        std::cout << "eig1: " << eig1 << "  eig2: " << eig2 << "    eig3: " << eig3 << std::endl;
//        std::cout << "eig3: " << eig3 << std::endl;
//        std::cout << EV << std::endl << std::endl;

        double norm = sqrt(EV(0,0)*EV(0,0)+EV(1,0)*EV(1,0)+EV(2,0)*EV(2,0));

        // save result as a map of Delaunay vertice and eigenvalue
//        all_vertices[vft]=std::make_pair(Point(EV(0,0)/norm,EV(1,0)/norm,EV(2,0)/norm), eig3);
        // save the result directly in the vertex_base of the Delaunay
        vft->info().sigma = eig3;
        all_eig3s.push_back(eig3);
    }

    size_t n = all_eig3s.size() / 2;
    nth_element(all_eig3s.begin(), all_eig3s.begin()+n, all_eig3s.end());

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Calculated noise per point with PCA on " << NN << " neighbors in " << duration.count() << "s" << std::endl;
    std::cout << "Median noise is " << all_eig3s[n] << std::endl;

    return all_eig3s[n];

}

// PCA with Delaunay neighborhood
void pcaDt(Delaunay& Dt){

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Start PCA..." << std::endl;

    Delaunay::Finite_vertices_iterator vft;
    for(vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){

        std::vector<Vertex_handle> av;
        Dt.adjacent_vertices(vft, std::back_inserter(av));

        std::vector<Point> adj_points;
        std::vector<Vertex_handle>::iterator avt;
        int nn = 0;
        for(avt = av.begin(); avt < av.end(); avt++){
            if(!Dt.is_infinite(*avt)){
                adj_points.push_back(Dt.point(*avt));
                nn++;
            }
        }


        Point centroid = CGAL::centroid(adj_points.begin(), adj_points.end(), CGAL::Dimension_tag<0>());
        Eigen::MatrixXd m(3,nn);
        Vector v;
        for(int i = 0; i < nn; i++){
            v = adj_points[i]-centroid;
            m(0,i) = v.x();
            m(1,i) = v.y();
            m(2,i) = v.z();
        }
        if(nn<1)
            std::cout<<"no neighbours" << std::endl;
        Eigen::MatrixXd A(3,3);
        A = (m*(m.transpose()))/nn;

        // compute eigenvalues of the covariance matrix A:
        // implemented from here: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
        double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
        double eig1;
        double eig2;
        double eig3;
        double q,p2,p,r,phi;
        Eigen::MatrixXd B(3,3);
        Eigen::MatrixXd I(3,3);
        I = I.setIdentity();
        if (p1 <= 1e-8){
           eig1 = A(0,0);
           eig2 = A(1,1);
           eig3 = A(2,2);}
        else{
           q = A.trace()/3.0;
           p2 = pow((A(0,0) - q),2.0) + pow((A(1,1) - q),2.0) + pow((A(2,2) - q),2.0) + 2.0 * p1;
           p = sqrt(p2 / 6.0);
           B = (1.0 / p) * (A - q * I);
           r = B.determinant() / 2.0;
           if (r <= -1.01)
              phi = M_PI / 3.0;
           else if (r >= 0.99)
              phi = 0.0;
           else
              phi = acos(r) / 3.0;
           // the eigenvalues satisfy eig3 <= eig2 <= eig1
           eig1 = q + 2.0 * p * cos(phi);
           eig3 = q + 2.0 * p * cos(phi + (2.0*M_PI/3.0));
           eig2 = 3.0 * q - eig1 - eig3;
        }
//        eig1 = eig1/(eig1+eig2+eig3);
//        eig2 = eig2/(eig1+eig2+eig3);
//        eig3 = eig3/(eig1+eig2+eig3);
        // compute eigenvector of the third (smallest) eigenvalue:
        Eigen::MatrixXd EV(3,3);
        EV = (A-eig1*I)*(A-eig2*I);
//        std::cout << "eig1: " << eig1 << "  eig2: " << eig2 << "    eig3: " << eig3 << std::endl;
//        std::cout << "eig3: " << eig3 << std::endl;
//        std::cout << EV << std::endl << std::endl;

        double norm = sqrt(EV(0,0)*EV(0,0)+EV(1,0)*EV(1,0)+EV(2,0)*EV(2,0));

        if(eig3 != eig3){
            std::cout << "eig3 changed from NaN to 0" << std::endl;
            eig3 = 0.0;
        }

//        all_vertices[vft]=std::make_pair(Point(EV(0,0)/norm,EV(1,0)/norm,EV(2,0)/norm), eig3);
//        vft->info().sigma = vft->info().sigma*eig3;
        vft->info().sigma = eig3;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Calculated noise per point with PCA on Delaunay neighborhood in " << duration.count() << "s" << std::endl;
}


////////////////////////////////////////////////////////////
/////////////////// preprocessing functions ////////////////
////////////////////////////////////////////////////////////
// estimate normals of a point set
typedef std::pair<Point, Vector> PointVectorPair;
// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
std::vector<PointVectorPair> estimateNormalsFun(const std::vector<Point>& points)
{

    // initialise a vector of point-vector-pairs
    std::vector<PointVectorPair> pointVectorPairs(points.size());

    // add the points as the first element of the point vector pair
    for(std::size_t i=0; i < points.size(); ++i)
    {
        pointVectorPairs[i].first = points[i];
    }

    // following two blocks are from this CGAL example: https://doc.cgal.org/latest/Point_set_processing_3/Point_set_processing_3_2normals_example_8cpp-example.html
    // Estimates normals direction.
    // Note: pca_estimate_normals() requires a range of points
    // as well as property maps to access each point's position and normal.
    // TODO: estimate the neighbourhood size automatically with: https://doc.cgal.org/latest/Point_set_processing_3/index.html#Point_set_processing_3Scale
    const int nb_neighbors = 4;
    CGAL::pca_estimate_normals<Concurrency_tag>
      (pointVectorPairs, nb_neighbors,
       CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
       normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

    // Orients normals.
    // Note: mst_orient_normals() requires a range of points
    // as well as property maps to access each point's position and normal.
    std::vector<PointVectorPair>::iterator unoriented_points_begin =    // this returns unoriented points that could be removed in the next step
            CGAL::mst_orient_normals(pointVectorPairs, nb_neighbors,
                               CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                               normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

    std::cout << "Normals calculated!" << std::endl;
    return pointVectorPairs;
};
