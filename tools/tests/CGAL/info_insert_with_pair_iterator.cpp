// from https://doc.cgal.org/latest/Triangulation_3/Triangulation_3_2info_insert_with_pair_iterator_8cpp-example.html
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/point_generators_3.h>
#include <vector>
#include <iterator>
typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point                                             Point;
int main()
{
  const int NUM_INSERTED_POINTS = 4;
  CGAL::Random_points_in_cube_3<Point> rnd(1.);
  // Construction from a vector of 1,000,000 points
  
  std::vector<std::pair<Point, unsigned>> V;
  V.reserve(NUM_INSERTED_POINTS);
  for (unsigned i = 0; i != 3; ++i)
  {
    V.push_back(std::make_pair(*rnd++, i));
    std::cout << "added point with index: " << i << std::endl;
  }
  std::cout << "manually added point with index: 3"  << std::endl;
  //add first point a second time
  
  double change = 0.0000000000000001;//still accepts 4 points
  /** 
   * point 0 and point 4 are treated as the same point!!!
   * assert fails
   */
  //double change = 0.000000000000000001;
  
  V.push_back(std::make_pair(Point(change+V[0].first[0],change+V[0].first[1],change+V[0].first[2]), 3));
  
  Delaunay T( V.begin(),V.end() );
  assert(T.is_valid());
  CGAL_assertion( T.number_of_vertices() == NUM_INSERTED_POINTS);
  // check that the info was correctly set.
  Delaunay::Finite_vertices_iterator vit;
  std::vector<Vb::Vertex_handle> vn;
  

  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
  {
    //T.adjacent_vertices(vit, std::back_inserter(vn));
    if( V[ vit->info() ].first != vit->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  std::cout << "iterate over neighbours:" << std::endl;
  
  std::vector<Delaunay::Vertex_handle> neighbour_handles;
  for ( vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
  {
    unsigned k_mt = vit->info();
    std::cout << "consider vertex: " << k_mt << std::endl;
    T.finite_adjacent_vertices(vit, std::back_inserter(neighbour_handles));
    std::cout << "neighbours: ";
    for( auto element: neighbour_handles)
    {
      std::cout << element->info() << ", ";
    }
    std::cout << std::endl;
    neighbour_handles.clear();
  }
  std::cout << "OK" << std::endl;
  /** result:
   * indices do not appear twice!
   */
  return 0;
}
