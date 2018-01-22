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
  const int NUM_INSERTED_POINTS = 2;
  CGAL::Random_points_in_cube_3<Point> rnd(1.);
  // Construction from a vector of 1,000,000 points
  //std::vector<Point> V;
  std::vector<std::pair<Point, unsigned>> V;
  V.reserve(NUM_INSERTED_POINTS);
  for (int i = 0; i != NUM_INSERTED_POINTS; ++i)
    V.push_back(std::make_pair(*rnd++, (unsigned)i));
  
  Delaunay T( V.begin(),V.end() );
  assert(T.is_valid());
  CGAL_assertion( T.number_of_vertices() == NUM_INSERTED_POINTS );
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
  std::cout << "OK" << std::endl;
  return 0;
}
