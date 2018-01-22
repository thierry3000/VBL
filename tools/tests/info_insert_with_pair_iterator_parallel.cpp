#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>

#include <iostream>
#include <fstream>


#include "info_insert_with_iterator_parallel.h"
int main()
{
// #pragma omp parallel
//   {
//     #pragma omp for
//     for(int i=0;i<8;i++)
//     {
//       std::cout<<"check if omp works at all?" << std::endl;
//     }
//   }
  
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
  Delaunay::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1, -1, -1, 1., 1., 1.), 20);
  // Construct the triangulation in parallel
  Delaunay T(V.begin(), V.end());
  //Triangulation T(V.begin(), V.end());
  assert(T.is_valid());
  CGAL_assertion( T.number_of_vertices() == NUM_INSERTED_POINTS );
#if 1
  //this vector stores all the neigbourhood information
  
  Delaunay::Finite_vertices_iterator vit;

  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
  {
    int k = vit->info();
    arr_vn[k] = std::shared_ptr<std::vector<Vertex_handle>>(new std::vector<Vertex_handle>);
    //std::vector<Vertex_handle> *vn = new std::vector<Vertex_handle>;
    T.incident_vertices(vit, std::back_inserter(*arr_vn[k]));	// list of neighbors
    if( V[ vit->info() ].first != vit->point() ){
      std::cerr << "Error different info" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  //now we try to access in parallel
#pragma omp parallel for
  for ( int i=0; i<NUM_INSERTED_POINTS; i++)
  {
    int size_of_neighbourhood = (*arr_vn[i]).size();
    //std::cout << size_of_neighbourhood << std::endl;
  }
  #endif
  std::cout << "OK" << std::endl;
#endif //CGAL_LINKED_WITH_TBB
  return 0;
}
