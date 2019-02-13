#include <CGAL/Epick_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/predicates_d.h>

int wasure_data::dump_surface(std::vector<Facet> & lft, int nblabs, std::string fileName){
  int dim = 3;
  std::map<typename DTW::Vertex_handle, uint> vertex_map;
  std::vector<Point> point_list;
  int acc = 0;
  Vertex_iterator fvit;

  for(auto fit = lft.begin();  fit != lft.end(); ++fit){
    Facet ft = *fit;
    Full_cell_handle fch = ft.full_cell();
    int idx = ft.index_of_covertex();
    Full_cell_handle fchn = fch->neighbor(idx);

    for(int i = 0; i < dim+1; ++i){
      if(i != ft.index_of_covertex()){
	Vertex_handle v = ft.full_cell()->vertex(i);
	if(vertex_map.find(v) == vertex_map.end()){
	  vertex_map[v] = acc++;
	  point_list.push_back(v->point());
	}
      }
    }
  }

  int nbv = vertex_map.size();
  int nbf = lft.size();


  std::ofstream fo;
  fo.open (fileName.c_str());
  std::cout << "\t Writing " << fileName << "..." << std::endl;
  
  fo << "ply" << std::endl;
  fo << "format ascii 1.0" << std::endl;
  fo << "comment VCGLIB generated" << std::endl;
  fo << "element vertex " << nbv << std::endl;
  fo << "property float x" << std::endl;
  fo << "property float y" << std::endl;
  fo << "property float z" << std::endl;
  fo << "property uchar red" << std::endl;
  fo << "property uchar green" << std::endl;
  fo << "property uchar blue" << std::endl;
  fo << "element face " << nbf << std::endl;
  fo << "property list uchar int vertex_indices" << std::endl;
  fo << "property uchar red" << std::endl;
  fo << "property uchar green" << std::endl;
  fo << "property uchar blue" << std::endl;
  fo << "end_header" << std::endl;
  fo << std::setprecision(12);
  std::cout << "wasure_data:dumping point...." << std::endl;

  for(auto const &pi : point_list) {
    fo << pi[0] << " " <<  pi[1]  << " " <<  pi[2]  << " " <<  255  <<  " " << 255 << " " << 255 << std::endl;
  }

  std::cout << "wasure_data:dumping facets...." << std::endl;
  for(auto fit = lft.begin();  fit != lft.end(); ++fit){
    Facet ft = *fit;
    int acc = 0;
    Full_cell_handle fch = ft.full_cell();
    int idx = ft.index_of_covertex();
    Full_cell_handle fchn = fch->neighbor(idx);
    int lab = fch->data().lab;
    int nlab = fchn->data().lab;
    double plab = ((double)lab/((double)(nblabs-1)));
    double pnlab = ((double)nlab/((double)(nblabs-1)));
    std::vector<int> lid;
    std::vector<Point> lp;
    for(int i = 0 ; i < dim +1; i++){
      if(i != ft.index_of_covertex()){
	lid.push_back(vertex_map[ft.full_cell()->vertex(i)]);
	lp.push_back(ft.full_cell()->vertex(i)->point());
      }
    }
    lid.push_back(ft.index_of_covertex());
    lp.push_back(ft.full_cell()->vertex(ft.index_of_covertex())->point());
    //    assert(vertex_map.find(ft.full_cell()->vertex(ft.index_of_covertex())) != vertex_map.end());
    assert(lid.size() == 4 && lp.size() == 4);
    int cr = 200;
    int cg = fabs(plab - pnlab)*255;
    int cb = fabs(plab - pnlab)*255;

    // Rustine
    bool is_inf = false;
    for(auto ppp : lp){
      std::cout << ppp << std::endl;
      if(ppp.dimension() < dim){
        std::cout << "WARNING, infinit point found" << std::endl;
        is_inf = true;
        break;
      }
    }
    
    int o1 = !is_inf ? CGAL::orientation(lp.begin(),lp.end()) :  0 ;
    bool bl =  ((o1 == -1 && pnlab >= 0.5) || ( o1 == 1 && pnlab < 0.5));
    if(bl)
      fo << "3 " << lid[0] << " " << lid[2] << " " << lid[1] << " " << cr << " " << cg << " " << cb << std::endl;
    else
      fo << "3 " << lid[0] << " " << lid[1] << " " << lid[2] << " " << cr << " " << cg << " " << cb << std::endl;
  }
  std::cout << "dumping ok" << std::endl;
  fo.close();
  return 0;
}


