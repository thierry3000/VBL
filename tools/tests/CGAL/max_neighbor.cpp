// from https://doc.cgal.org/latest/Triangulation_3/Triangulation_3_2info_insert_with_pair_iterator_8cpp-example.html
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/point_generators_3.h>
#include <vector>
#include <iterator>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point                                             Point;
int main()
{
  const int NUM_INSERTED_POINTS = 1e6;
  CGAL::Random_points_in_cube_3<Point> rnd(1.);
  // Construction from a vector of 1,000,000 points
  std::ofstream myfile;
  myfile.open ("neighbours.txt");
  const int howManyTimesToCheck = 200;
  
  for(int i =0; i<howManyTimesToCheck; ++i)
  {
    std::vector<std::pair<Point, unsigned>> V;
    V.reserve(NUM_INSERTED_POINTS);
    for (unsigned i = 0; i != NUM_INSERTED_POINTS; ++i)
    {
      V.push_back(std::make_pair(*rnd++, i));
      //std::cout << "added point with index: " << i << std::endl;
    }
    
    Delaunay T( V.begin(),V.end() );
    assert(T.is_valid());
    CGAL_assertion( T.number_of_vertices() == NUM_INSERTED_POINTS);
    // check that the info was correctly set.
    Delaunay::Finite_vertices_iterator vit;
    std::vector<Vb::Vertex_handle> vn;
    
    std::vector<Delaunay::Vertex_handle> neighbours;
    int max_neighbours = 0;
    
    for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    {
      T.finite_adjacent_vertices(vit, std::back_inserter(neighbours));
      if(neighbours.size() > max_neighbours)
      {
	max_neighbours = neighbours.size();
      }
      //T.adjacent_vertices(vit, std::back_inserter(vn));
      if( V[ vit->info() ].first != vit->point() ){
	std::cerr << "Error different info" << std::endl;
	exit(EXIT_FAILURE);
      }
      neighbours.clear();
    }
    std::cout << max_neighbours << std::endl;
    myfile << max_neighbours << "\n";
  }
  
  //std::cout << "iterated over neighbours:" << std::endl;
  
  myfile.close();
//   std::vector<Delaunay::Vertex_handle> neighbour_handles;
//   for ( vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
//   {
//     unsigned k_mt = vit->info();
//     std::cout << "consider vertex: " << k_mt << std::endl;
//     T.finite_adjacent_vertices(vit, std::back_inserter(neighbour_handles));
//     std::cout << "neighbours: ";
//     for( auto element: neighbour_handles)
//     {
//       std::cout << element->info() << ", ";
//     }
//     std::cout << std::endl;
//     neighbour_handles.clear();
//   }
  
  
  std::cout << "OK" << std::endl;
  /** result:
   * maximum is 82
   */
  return 0;
}
