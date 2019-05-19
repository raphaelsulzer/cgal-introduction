//int exportTriWithCnFun(const char* ifn, const char* ofn, std::vector<PointVectorPair> pVP)
//Delaunay exportTriWithCnFun(std::vector<PointVectorPair> pVP, const char* ofn)
Delaunay exportTriWithCnFun(std::vector<PointVectorPair> pVP, const char* ofn)
{

    Delaunay Dt(pVP.begin(), pVP.end());

    // get number of vertices and triangles of the triangulation
    Delaunay::size_type nv = Dt.number_of_vertices();
    Delaunay::size_type nf = Dt.number_of_finite_facets();

    // create PLY output file for outputting the triangulation, with point coordinates, color, normals and triangle facets
    std::ofstream fo;
    fo.open (ofn);
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "comment VCGLIB generated" << std::endl;
    fo << "element vertex " << nv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property float nx" << std::endl;
    fo << "property float ny" << std::endl;
    fo << "property float nz" << std::endl;
    fo << "element face " << nf << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(3);

    // give every vertex from the triangulation an index starting at 0
    // and already print the point coordinates, color and normal of the vertex to the PLY file
    std::map<Vertex_handle, int> Vertices;
    int index = 0;
    Delaunay::Finite_vertices_iterator vft;
    for (vft = Dt.finite_vertices_begin() ; vft != Dt.finite_vertices_end() ; vft++){
        // assign index to vertex handle. this is needed later for finding the index a certain cell is constructed with
        // TODO: maybe this could already be included in the info vector of every point/vertex, and then instead of using the finite_vertices_iterator,
        // one could iterate over the index vector, so that the point order of the file is determined by the info vector
        // and you would do the find operation to find the point in a std::map<Point, idx>
        // update: don't see anymore what I would gain from that
        Vertices[vft] = index;
        // print data to file
        fo << vft->point() << " "                           // coordinates
           << vft->info() << std::endl;                     // normal
        //        fo << vft->point() << " "                           // coordinates
        //           << int(std::get<1>(vft->info())[0]) << " "       // red
        //           << int(std::get<1>(vft->info())[1]) << " "       // green
        //           << int(std::get<1>(vft->info())[2]) << " "       // blue
        //           << std::get<0>(vft->info()) << std::endl;        // normal
        index++;
    }

    // Save the facets to the PLY file
    //std::cout << "iterate over finite triangles: " << std::endl;
    Delaunay::Finite_facets_iterator fft;
    int vidx;
    // initialise cell and vertex handle
    Cell_handle c;
    Vertex_handle v;
    for(fft = Dt.finite_facets_begin() ; fft != Dt.finite_facets_end() ; fft++){

        // get vertex and cell index that describes the facet:
        // facet fft is represented by std::pair(cell c, int vidx). vidx is the vertex opposite to the cell.
        // even though some of the facets may be described by infinite cells, the facet is still has a neighbouring cell that is finite.
        // see: https://doc.cgal.org/latest/Triangulation_3/index.html
        c = fft->first;         // cell
        vidx = fft->second;     // vertex index

        // start printed facet line with a 3
        fo << 3 << ' ';
        // if opposite vertex vidx is 2, we start at j = vidx + 1 = 3, 3%4 = 3
        // next iteration: j = 4, 4%4 = 0, next iteration: j = 5, 5%4 = 1;
        // so we exactely skip 2 - the opposite vertex.
        for(int j = vidx + 1 ; j <= vidx + 3 ; j++){
            
            v = c->vertex(j%4);
            // print the indicies of each cell to the file
            // vertices is a map of all vertices of the triangulation to an index
            fo << Vertices.find(v)->second << ' ';

        }
        // new cell
        fo << std::endl;
    }
    fo.close();

    std::cout << "Delaunay triangulation done and exported to PLY file!" << std::endl;
    return Dt;

}