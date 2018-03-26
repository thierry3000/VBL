
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <tbb/concurrent_vector.h>
#include <omp.h>
#include <stdlib.h>     //for using the function sleep
#include <random>

# include <boost/random.hpp>
boost::random::mt19937 rng;         // produces randomness out of thin air
                                    // see pseudo-random number generators

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K>      Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;

int intRand(const int & min, const int & max) {
    static thread_local std::mt19937 generator;
    //std::uniform_int_distribution<int> distribution(min,max);
    std::uniform_real_distribution<float> distribution(min,max);
    return distribution(generator);
}

int main()
{
  int iam = 0, np = 1;
  #pragma omp parallel default(shared) private(iam, np)
  {
    #if defined (_OPENMP)
      np = omp_get_num_threads();
      iam = omp_get_thread_num();
    #endif
    printf("Hello from thread %d out of %d\n", iam, np);
  }
  // construction from a list of points :
  //std::list<Point> L;
#if 0
  tbb::concurrent_vector<double> bla;
#pragma omp parallel for
  for( int k = 0;k<10000000;++k)
  {
    int x = k;
    int y = k;
    int z = k;
    #if defined (_OPENMP)
      np = omp_get_num_threads();
      iam = omp_get_thread_num();
    #endif
    x = np*k;
    z = iam*k;
    //printf("Hello from thread %d out of %d\n", iam, np);
    //sleep(1000);
    bla.push_back(x+z);
  }
#endif
  //std::list<Point> L;
  tbb::concurrent_vector<Point> L;
  //insert n points
  uint number_of_points = 10000000;
  double scale = 1000;
  boost::random::uniform_01<> six_f;
  //printf("here");
#pragma omp parallel for
  for( int k = 0;k<number_of_points;++k)
  {
//     int x = six_f(rng)*scale;
//     int y = six_f(rng)*scale;
//     int z = six_f(rng)*scale;
    int x = scale*intRand(0,1);
    int y = scale*intRand(0,1);
    int z = scale*intRand(0,1);
    //printf("here2\n");
    L.push_back(Point(x,y,z));
  }
  printf("finished creating points\n");
  std::cout.flush();
  Triangulation T(L.begin(), L.end());
  printf("finished creating Triangulation\n");
  Triangulation::size_type n = T.number_of_vertices();
//   // insertion from a vector :
//   std::vector<Point> V(3);
//   V[0] = Point(0,0,1);
//   V[1] = Point(1,1,1);
//   V[2] = Point(2,2,2);
//   n = n + T.insert(V.begin(), V.end());
  assert( n == number_of_points );       // 6 points have been inserted
  assert( T.is_valid() ); // checking validity of T
  Locate_type lt;
  int li, lj;
  Point p(0,0,0);
  Cell_handle c = T.locate(p, lt, li, lj);
  // p is the vertex of c of index li :
  assert( lt == Triangulation::VERTEX );
  assert( c->vertex(li)->point() == p );
  Vertex_handle v = c->vertex( (li+1)&3 );
  // v is another vertex of c
  Cell_handle nc = c->neighbor(li);
  // nc = neighbor of c opposite to the vertex associated with p
  // nc must have vertex v :
  int nli;
  assert( nc->has_vertex( v, nli ) );
  // nli is the index of v in nc
  std::ofstream oFileT("output",std::ios::out);
  // writing file output;
  oFileT << T;
  Triangulation T1;
  std::ifstream iFileT("output",std::ios::in);
  // reading file output;
  iFileT >> T1;
  assert( T1.is_valid() );
  assert( T1.number_of_vertices() == T.number_of_vertices() );
  assert( T1.number_of_cells() == T.number_of_cells() );
  return 0;
}
