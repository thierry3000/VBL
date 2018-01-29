#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Compact_container.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>

#include "tbb/concurrent_vector.h"

const int NUM_INSERTED_POINTS = 1e6;

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

void complicated_calculation(Triangulation::Vertex_handle *vn)
{
  unsigned long k = (*vn)->info();
  unsigned long c = 2*k;
  std::cout << " doing a complicated calculation on: " << (*vn)->info() << std::endl;
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
void complicated_calculation2(Triangulation::Vertex_handle &vn)
{
  unsigned long k = vn->info();
  unsigned long c = 2*k;
  if( k < (unsigned long) NUM_INSERTED_POINTS)
  {
//     std::cout << " doing a complicated calculation on: " << vn->info() << std::endl;
//     std::cout.flush();
    printf(" doing a complicated calculation on: %i\n", vn->info());
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
}

class ApplyFoo {
    //Triangulation::Vertex_handle *const my_a;
    CGAL::Compact_container<Triangulation::Vertex_handle> const *my_a;
public:
    void operator()( const tbb::blocked_range<size_t>& r ) const 
    {
        for( size_t i=r.begin(); i!=r.end(); ++i )
        {
          Triangulation::Vertex_handle a = (*my_a)[i];
           complicated_calculation2(a);
        }
    }
    void apply_complicated(int i)
    {
      Triangulation::Vertex_handle a = (*my_a)[i];
      complicated_calculation2(a);
    }
    ApplyFoo( CGAL::Compact_container<Triangulation::Vertex_handle> *a ):my_a(a){}
};



void ParallelApplyFoo( CGAL::Compact_container<Triangulation::Vertex_handle> *a)
{
  int n = a->size();
  tbb::parallel_for(tbb::blocked_range<size_t>(0,n), ApplyFoo(a));
}
// void ParallelApplyFoo2( CGAL::Compact_container<Triangulation::Vertex_handle> &a)
// {
//   int n = a.size();
//   tbb::parallel_for(tbb::blocked_range<size_t>(0,n), ApplyFoo(a));
// }

// void SerialApplyFoo( CGAL::Compact_container<Triangulation::Vertex_handle*> *a)
// {
//   int n = a->size();
//   for(int p=0; p<n; p++)
//   {
//     ApplyFoo(a).apply_complicated(p);
//   }
// }

void SerialApplyFoo2( CGAL::Compact_container<Triangulation::Vertex_handle> *a)
{
  for(int p=0; p<a->size(); ++p)
  {
    std::cout << "p : " << p << std::endl;
    complicated_calculation2( (*a)[p]);
  }
}

int main()
{
#ifdef CGAL_LINKED_WITH_TBB
  
  /* NUM_INSERTED_POINTS = 15 is still finishing but strange index output
   * NUM_INSERTED_POINTS = 16 -> Segmentation fault
   * gdb output:
   * 0x0000555555644b6b in complicated_calculation (vn=...)
   * at /home/usersHR/thierry/git_codes/VBL/tools/tests/parallel_insertion_in_delaunay_3.cpp:38
   * 38	  std::cout << " doing a complicated calculation on: " << vn->info() << std::endl;
   * ***
   * vn->info() is not accessible
   */
  
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
  std::array<Triangulation::Vertex_handle, NUM_INSERTED_POINTS> myContainer2;
  std::array<Triangulation::Vertex_handle*, NUM_INSERTED_POINTS> myContainer3;
  tbb::concurrent_vector<Triangulation::Vertex_handle> *myContainer4 = new tbb::concurrent_vector<Triangulation::Vertex_handle>;
  CGAL::Compact_container<Triangulation::Vertex_handle> *myContainer5 = new CGAL::Compact_container<Triangulation::Vertex_handle>();
  
  myContainer4->resize(NUM_INSERTED_POINTS);
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
  {
    if( V[ vit->info() ].first != vit->point() )
    {
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
    myContainer.insert(vit);
    Triangulation::Vertex_handle vertex=vit;
    
    myContainer3[vit->info()] = &vertex;
    myContainer2[vit->info()] = vit;
    //myContainer4.push_back(&vertex);
    //std::cout << vit->info() << std::endl;
    myContainer4->at(vit->info()) = vertex;
    myContainer5->insert(vit);
    
    //complicated_calculation(&vertex);
  }
  Triangulation::Vertex_handle bla = *myContainer3[42];
  unsigned long q = bla->info();
  //std::cout << bla->info() << std::endl;
  ParallelApplyFoo(myContainer5);
  //SerialApplyFoo2(myContainer5);
  
  std::cout << "OK" << std::endl;
#endif //CGAL_LINKED_WITH_TBB
  return 0;
}

