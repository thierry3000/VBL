#ifndef INFO_INSERT_WITH_ITERATOR_PARALLEL
#define INFO_INSERT_WITH_ITERATOR_PARALLEL

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Exact_predicates_inexact_constructions_kernel     	K;

  typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> 	Vb;
  typedef CGAL::Fixed_alpha_shape_vertex_base_3<K, Vb>      			AsVb;
  typedef CGAL::Fixed_alpha_shape_cell_base_3<K>        				Fb;
  typedef CGAL::Triangulation_data_structure_3<AsVb,Fb,CGAL::Parallel_tag> Tds;
  typedef CGAL::Delaunay_triangulation_3<K,Tds> Delaunay;
  typedef CGAL::Fixed_alpha_shape_3<Delaunay>Fixed_alpha_shape_3;

  // ****************

  typedef K::Point_3                                          Point;
  typedef Fixed_alpha_shape_3::Vertex_handle					Vertex_handle;
  typedef Delaunay::Finite_vertices_iterator			Finite_vertices_iterator;

const int NUM_INSERTED_POINTS = 2e5;

std::array<std::shared_ptr<std::vector<Vertex_handle>>, NUM_INSERTED_POINTS> arr_vn; // vector of Vertex_handle's of neighbors
#endif
#endif
