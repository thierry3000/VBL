#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/point_generators_3.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  // Delaunay T3
  typedef CGAL::Triangulation_data_structure_3< 
    CGAL::Triangulation_vertex_base_with_info_3<unsigned,K>, 
    CGAL::Triangulation_cell_base_3<K>, 
    CGAL::Parallel_tag>                          Tds;
  typedef CGAL::Delaunay_triangulation_3<K, Tds> Triangulation;
  typedef Triangulation::Point          Point;

#endif
int binomialCoeff(int n, int k)
{
    int res = 1;
    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;
    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}

void complicated_calculation(Triangulation::Vertex_handle vn)
{
  std::cout << " doing a complicated calculation on: " << vn->info() << std::endl;
  int b=0;
  //sleep(20);
  const int zahl = 200;
  for (int p=0;p<zahl;p++)
  {
    for(int t=0; t<zahl; t++)
    {
      for(int i=0;i< zahl; i++)
      {
        b=b+binomialCoeff(t,i);
        //std::cout << b <<std::endl;
      }
    }
  }
}
class ApplyFoo {
    Triangulation::Vertex_handle *const my_a;
public:
    void operator()( const tbb::blocked_range<size_t>& r ) const 
    {
        Triangulation::Vertex_handle *a = my_a;
        for( size_t i=r.begin(); i!=r.end(); ++i ) 
           complicated_calculation(a[i]);
    }
    ApplyFoo( Triangulation::Vertex_handle a[] ) :
        my_a(a)
    {}
};
void ParallelApplyFoo( Triangulation::Vertex_handle a[], size_t n)
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,n), ApplyFoo(a));
}

int main()
{
#ifdef CGAL_LINKED_WITH_TBB
  
  const int NUM_INSERTED_POINTS = 15;
  CGAL::Random_points_in_cube_3<Point> rnd(1.);
  // Construction from a vector of 1,000,000 points
  std::vector<std::pair<Point, unsigned>> V;
  V.reserve(NUM_INSERTED_POINTS);
  for (unsigned int i = 0; i != NUM_INSERTED_POINTS; ++i)
    V.push_back(std::make_pair(*rnd++, i));
  
  // Construct the locking data-structure, using the bounding-box of the points
  Triangulation::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);
  // Construct the triangulation in parallel
  Triangulation T(V.begin(), V.end(), &locking_ds);
  assert(T.is_valid());
  
  // Triangulation structure is filled, now traverse it
  // this works and is from the examples!
  Triangulation::Finite_vertices_iterator vit;
  CGAL::Compact_container<Triangulation::Vertex_handle> myContainer;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
  {
    if( V[ vit->info() ].first != vit->point() )
    {
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
    myContainer.insert(vit);
    //complicated_calculation(vit);
  }
  tbb::parallel_for(tbb::blocked_range<size_t>(0,myContainer.size()), ApplyFoo(&myContainer[0]));
  std::cout << "OK" << std::endl;
#endif //CGAL_LINKED_WITH_TBB
  return 0;
}

