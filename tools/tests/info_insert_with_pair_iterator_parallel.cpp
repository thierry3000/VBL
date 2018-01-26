#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>


#include <CGAL/Concurrent_compact_container.h>

#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>

#include <iostream>
#include <fstream>


#include "info_insert_with_iterator_parallel.h"

#include "tbb/tbb.h"

void Foo(int &i)
{
  std::cout << "in Foo i: " << i << std::endl;
  sleep(5);
}
void SerialApplyFoo( int a[], size_t n)
{
  for( size_t i=0; i!=n; ++i)
  {
    Foo(a[i]);
  }
}
class ApplyFoo{
  int *const my_a;
public:
  void operator()(const tbb::blocked_range<size_t>&r) const{
    int *a = my_a;
    for(size_t i=r.begin(); i!=r.end(); ++i)
    {
      std::cout << "apply Foo in parallel for i= " << i << std::endl;
      Foo(a[i]);
    }
  }
  ApplyFoo( int a[]):
    my_a(a)
    {}
};
void ParallelApplyFoo(int a[], size_t n)
{
  parallel_for(tbb::blocked_range<size_t>(0,n), ApplyFoo(a));
};

void ApplyGeometricCalculation::operator()(const tbb::blocked_range<size_t> &r) const
{
  Vertex_handle *a = my_current_vertex_handle;
  
  for( size_t i=r.begin(); i!=r.end(); ++i)
    apply_geometry(a[i], p_DelTri);
};
void ApplyGeometricCalculation2::operator()(const tbb::blocked_range<size_t> &r) const
{
  Vertex_handle *a = my_current_vertex_handle;
  
  for( size_t i=r.begin(); i!=r.end(); ++i)
    apply_geometry(a[i], p_DelTri);
};
void apply_geometry(Vertex_handle &i, Delaunay *p_DelTri)// Foo
{
  if(i->is_valid())
  {
    std::cout << 6 << std::endl;
  }
//   if(!p_DelTri->is_infinite(i))
//   {
//     printf("I apply geometry calculations on vertex %i\n", i->info());
//   }
//   else
//   {
//     printf("INFINTE at: %i\n", i->info());
//   }
  
};
void ParallelApplyCalculation( Vertex_handle a[], Delaunay *_p, size_t n)
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,n), ApplyGeometricCalculation(a,_p));
}
int main()
{
  //std::array<int,10> aVector;
  const int n = 42;
  std::vector<int> aVector;
  for( int i=0;i<42;i++)
  {
    //aVector[i] = 2*i;
    aVector.push_back(2*i);
  }
  //ParallelApplyFoo(&aVector[0], aVector.size());
  //tbb::parallel_for(tbb::blocked_range<size_t>(0,aVector.size()), ApplyFoo(&aVector[0]));
  //SerialApplyFoo(&aVector[0], aVector.size());
  
  
#ifdef CGAL_LINKED_WITH_TBB
  
  //note for 5e5 cells we have an memory increase of few hundred MB
  
  CGAL::Random_points_in_cube_3<Point> rnd(1.);
  // Construction from a vector of 1,000,000 points
  //std::vector<Point> V;
  std::vector<std::pair<Point, unsigned>> V;
  V.reserve(NUM_INSERTED_POINTS);
  for (int i = 0; i < NUM_INSERTED_POINTS; ++i)
    V.push_back(std::make_pair(*rnd++, i));
  
  // Construct the locking data-structure, using the bounding-box of the points
  Delaunay *T;
  Delaunay::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1, -1, -1, 1., 1., 1.), 50);
  // Construct the triangulation in parallel
  if(NUM_INSERTED_POINTS>10)
  {
    T = new Delaunay(V.begin(), V.end(), &locking_ds);
  }
  else
  {
    T =new Delaunay(V.begin(), V.end());
  }
  //Triangulation T(V.begin(), V.end());
  assert(T->is_valid());
  CGAL_assertion( T->number_of_vertices() == NUM_INSERTED_POINTS );
  printf("Triangulation created!\n");
#if 1
  //this vector stores all the neigbourhood information
  
  CGAL::Compact_container<Vertex_handle> myContainer;
  Delaunay::Finite_vertices_iterator vit;
  Delaunay::All_vertices_iterator vit2;
  CGAL::Compact_container<Vertex_handle> *myContainer2 = new CGAL::Compact_container<Vertex_handle>();
  
  //for (vit = T->finite_vertices_begin(); vit != T->finite_vertices_end(); ++vit)
  for (vit2 = T->all_vertices_begin(); vit2 != T->all_vertices_end(); vit2++)
  {
    int k = vit2->info();
    std::cout << "k: " << k << std::endl;
    
    //arr_vn[k] = std::shared_ptr<std::vector<Vertex_handle>>(new std::vector<Vertex_handle>);
    //std::vector<Vertex_handle> *vn = new std::vector<Vertex_handle>;
    //T->incident_vertices(vit, std::back_inserter(*arr_vn[k]));	// list of neighbors
#if 0
    if( V[ vit2->info() ].first != vit2->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    myContainer.insert(vit2);
    myContainer2->insert(vit2);
  }
  std::cout << "myContainer.size: " << myContainer.size() << std::endl;
  std::cout << "myContainer2->size: " << myContainer2->size() << std::endl;
  
//   for (int v=0; v< myContainer2->size(); v++)
//   {
//     printf("index: %i, value: %i \n",v, myContainer[v]->info());
//   }

  //now we try to access in parallel
// #pragma omp parallel for
//   for ( int i=0; i<NUM_INSERTED_POINTS; i++)
//   {
//     int size_of_neighbourhood = (*arr_vn[i]).size();
//     //std::cout << size_of_neighbourhood << std::endl;
//   }
  #pragma omp parallel for
    for ( int i=0; i<myContainer2->size(); i++)
    {
      int k = (*myContainer2)[i]->info();
//       std::cout << k << std::endl;
      printf("k: %i\n", k);
    }
  #endif
  //tbb::parallel_for(tbb::blocked_range<size_t>(0,myContainer2->size()), ApplyGeometricCalculation(&((*myContainer2)[0]),T));
  
  tbb::parallel_for(tbb::blocked_range<size_t>(0,myContainer.size()), ApplyGeometricCalculation2(&myContainer[0],T));
  
  std::cout << "OK" << std::endl;
#endif //CGAL_LINKED_WITH_TBB
  return 0;
}
