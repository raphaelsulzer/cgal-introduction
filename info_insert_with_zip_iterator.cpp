#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<std::tuple<int,int>, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;
typedef Delaunay::Point                                             Point;
int triangulationWithZip()
{

    std::vector<std::tuple<int, int>> index1;
    index1.push_back(std::make_tuple(1, 2));
    index1.push_back(std::make_tuple(1, 2));
    index1.push_back(std::make_tuple(1, 2));
    index1.push_back(std::make_tuple(1, 2));
    index1.push_back(std::make_tuple(1, 2));
    index1.push_back(std::make_tuple(1, 2));

//  std::vector<unsigned> index2;
//  index2.push_back(0);
//  index2.push_back(1);
//  index2.push_back(2);
//  index2.push_back(3);
//  index2.push_back(4);
//  index2.push_back(5);

//    std::vector< std::pair<Point,unsigned> > points;
//    points.push_back( std::make_pair(Point(0,0,0),0) );
//    points.push_back( std::make_pair(Point(1,0,0),1) );
//    points.push_back( std::make_pair(Point(0,1,0),2) );
//    points.push_back( std::make_pair(Point(0,0,1),3) );
//    points.push_back( std::make_pair(Point(2,2,2),4) );
//    points.push_back( std::make_pair(Point(-1,0,1),5) );

//    Delaunay T( points.begin(),points.end() );


    std::vector<Point> points;
    points.push_back(Point(0,0,0));
    points.push_back(Point(1,0,0));
    points.push_back(Point(0,1,0));
    points.push_back(Point(0,0,1));
    points.push_back(Point(2,2,2));
    points.push_back(Point(-1,0,1));

    Delaunay T( boost::make_zip_iterator(boost::make_tuple( points.begin(),index1.begin() )),
              boost::make_zip_iterator(boost::make_tuple( points.end(),index1.end() ) )  );




  // check that the info was correctly set.
  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
  {


      std::cout << "coords: " << vit->point() << "   index: " << std::get<0>(vit->info()) << std::endl;
//      std::cout << "coords: " << vit->point() << std::endl;

//      if( points[ vit->info() ] != vit->point() ){
//      std::cerr << "Error different info" << std::endl;


//      exit(EXIT_FAILURE);
//    }
  }



  return 0;
}
