#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Compact_container.h>
#include <CGAL/Concurrent_compact_container.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>

#include "tbb/concurrent_vector.h"
#include "tbb/tbb.h"
#include "tbb/parallel_do.h"

#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
boost::random::random_device gen;
boost::random::uniform_01<> uni_f;

double myRandom()
{
  return uni_f(gen);
}

const int NUM_INSERTED_POINTS = 10;

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

Triangulation *T;


void Foo2(const Triangulation::Vertex &iter)
{
  //Triangulation::Vertex_handle h = &iter;
  std::vector<Triangulation::Vertex_handle> vn;
  Triangulation::Vertex_handle h = T->insert(iter.point());
  T->finite_adjacent_vertices(h, std::back_inserter(vn));
  unsigned long k = iter.info();
  unsigned long c = 2*k;
//     std::cout << " doing a complicated calculation on: " << vn->info() << std::endl;
//     std::cout.flush();
  printf(" doing a complicated new calculation on: %i\n", iter.info());
  for( auto bla: vn)
  {
    printf("%i has neighbor %i\n", iter.info(), bla->info());
  }
  const int zahl = 2e6 ;
  std::vector<double> vec_of_rand;
  int i=0;
  for (int p=0;p<zahl;p++)
  {
    double ran = myRandom();
    vec_of_rand.push_back(ran);
  }
  std::sort(vec_of_rand.begin(), vec_of_rand.end());
}
void Foo(const Triangulation::Vertex_handle a)
{
  //Triangulation::Vertex_handle a = (*my_a)[i];
  unsigned long k = a->info();
  unsigned long c = 2*k;
//     std::cout << " doing a complicated calculation on: " << vn->info() << std::endl;
//     std::cout.flush();
  printf(" doing a complicated new calculation on: %i\n", a->info());
  const int zahl = 2e6 ;
  std::vector<double> vec_of_rand;
  int i=0;
  for (int p=0;p<zahl;p++)
  {
    double ran = myRandom();
    vec_of_rand.push_back(ran);
  }
  std::sort(vec_of_rand.begin(), vec_of_rand.end());
}
// void Foo(const Triangulation::Vertex a)
// {
//   //Triangulation::Vertex_handle a = (*my_a)[i];
//   unsigned long k = a.info();
//   unsigned long c = 2*k;
// //     std::cout << " doing a complicated calculation on: " << vn->info() << std::endl;
// //     std::cout.flush();
//   printf(" doing a complicated new calculation on: %i\n", a.info());
//   const int zahl = 2e5 ;
//   std::vector<double> vec_of_rand;
//   int i=0;
//   for (int p=0;p<zahl;p++)
//   {
//     double ran = myRandom();
//     vec_of_rand.push_back(ran);
//   }
//   std::sort(vec_of_rand.begin(), vec_of_rand.end());
// }

class ApplyFoo {
public:
    void operator()( const Triangulation::Vertex_handle &vit ) const;
    
    ApplyFoo( ){};
};
void ApplyFoo::operator()( const Triangulation::Vertex_handle &vit) const 
{
  Foo(vit);
}

void ParallelApplyToConcurrent(const CGAL::Concurrent_compact_container<Triangulation::Vertex_handle> &a)
{
  printf("starting ParallelApplyToConcurrent\n");
  tbb::parallel_do(a.begin(), a.end(), ApplyFoo());
}

class ApplyFoo2 {
public:
    void operator()( const Triangulation::Vertex &iter ) const;
    
    ApplyFoo2( ){};
};
void ApplyFoo2::operator()( const Triangulation::Vertex &iter) const 
{
  Foo2(iter);
}
void ParallelApply()
{
  printf("starting ParallelApply\n");
  tbb::parallel_do(T->finite_vertices_begin(), T->finite_vertices_end(), ApplyFoo2());
}

#if 0
class ApplyFooT {
public:
    void operator()( const Triangulation::Vertex_handle &vit ) const;
    
    ApplyFooT( ){};
};
void ApplyFooT::operator()( const Triangulation::Vertex_handle &vit) const 
{
  //Foo(vit);
}

// void ParallelApplyToTri(const Triangulation::Finite_vertices_iterator &list)
// {
//   printf("starting ParallelApplyToConcurrent\n");
//   tbb::parallel_do(list.begin(), list.end(), ApplyFooT());
// }
#endif

void SerialApplyFoo( const CGAL::Concurrent_compact_container<Triangulation::Vertex_handle> &a)
{
  for(auto b:a)
  {
    Foo(b);
  }
}
void SerialApplyFoo2( )
{
  for(auto vit=T->finite_vertices_begin(); vit!=T->finite_vertices_end(); ++vit)
  {
    Foo(vit);
  }
}
void cgal_concurrent_container_test()
{
  tbb::concurrent_vector<int> bla;
  for(int i =0; i<100; ++i)
  {
    bla.push_back(i);
  }
  for( auto a:bla)
  {
    int k=a*a;
  }
  printf("finished cgal_concurrent_container_test\n");
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
  T = new Triangulation(V.begin(), V.end(), &locking_ds);
  assert(T->is_valid());
  
  // Triangulation structure is filled, now traverse it
  // this works and is from the examples!
  //Triangulation::Finite_vertices_iterator vit;
  CGAL::Compact_container<Triangulation::Vertex_handle> myContainer;
  std::array<Triangulation::Vertex_handle, NUM_INSERTED_POINTS> myContainer2;
  std::array<Triangulation::Vertex_handle*, NUM_INSERTED_POINTS> myContainer3;
  tbb::concurrent_vector<Triangulation::Vertex_handle> *myContainer4 = new tbb::concurrent_vector<Triangulation::Vertex_handle>;
  CGAL::Compact_container<Triangulation::Vertex_handle> *myContainer5 = new CGAL::Compact_container<Triangulation::Vertex_handle>();
  CGAL::Concurrent_compact_container<Triangulation::Vertex_handle> *myContainer6 = new CGAL::Concurrent_compact_container<Triangulation::Vertex_handle>();
  CGAL::Concurrent_compact_container<Triangulation::Vertex_handle> myContainer7;
  
  myContainer4->resize(NUM_INSERTED_POINTS);
  for (auto vit = T->finite_vertices_begin(); vit != T->finite_vertices_end(); ++vit)
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
    myContainer6->insert(vit);
    myContainer7.insert(vit);
    
    //complicated_calculation(&vertex);
  }
  Triangulation::Vertex_handle bla = *myContainer3[0];
  unsigned long q = bla->info();
  
  auto abc = T->finite_vertices_begin();
  std::cout << abc->info() << std::endl;
#if 0
  tbb::parallel_do(T.finite_vertices_begin(), T.finite_vertices_end(), ApplyFoo());
#endif
#if 0
  SerialApplyFoo(myContainer7);
#endif
#if 0
  //result of sizeof is 48 which means 48 byte
  std::cout << "size of vertex handle: " << sizeof(*((*myContainer5)[0])) << std::endl;
  ParallelApplyToConcurrent(myContainer7);
#endif
  
#if 1
  ParallelApply();
#endif
#if 0
  SerialApplyFoo2();
#endif
  
#if 0
  ParallelApplyFoo(myContainer5);
#endif
#if 0
  cgal_concurrent_container_test();
#endif
  std::cout << "OK" << std::endl;
#endif //CGAL_LINKED_WITH_TBB
  return 0;
}

