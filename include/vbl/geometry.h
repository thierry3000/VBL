#ifndef GEOMETRY_H
#define GEOMETRY_H  // header guard
// header for the geometry part
#include <CGAL/Default.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>

#include <CGAL/Compact_container.h>
#include <CGAL/Concurrent_compact_container.h>
#include <boost/shared_ptr.hpp>
#include <tbb/mutex.h>

//#define myDebugComments
#define useSerialApproach

typedef CGAL::Exact_predicates_inexact_constructions_kernel     	K;

typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> 	Vb;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<K, Vb>      			AsVb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<K>        				Fb;

#ifdef useSerialApproach
  typedef CGAL::Triangulation_data_structure_3<AsVb,Fb>               Tds;
#else
  typedef CGAL::Triangulation_data_structure_3<AsVb,Fb,CGAL::Parallel_tag>               Tds;
#endif
typedef CGAL::Delaunay_triangulation_3<K,Tds>                      	Triangulation_3;
typedef CGAL::Fixed_alpha_shape_3<Triangulation_3>					Fixed_alpha_shape_3;

// ****************

typedef K::Point_3                                          Point;
typedef Fixed_alpha_shape_3::Vertex_handle					Vertex_handle;
typedef Triangulation_3::Finite_vertices_iterator			Finite_vertices_iterator;
typedef Triangulation_3::All_vertices_iterator			All_vertices_iterator;

std::vector< std::pair<Point,unsigned> > v;	// vector of points with info passed to CGAL



Triangulation_3 *DelTri;

#ifndef useSerialApproach
class ApplyGeometricCalculation{
public:
  friend class CellsSystem;
  void operator()(const Triangulation_3::Vertex_handle &item) const;
  ApplyGeometricCalculation( ){};
};

void apply_geometry(const Triangulation_3::Vertex_handle item);
void ParallelApplyGeometricCalculation( const CGAL::Concurrent_compact_container<Triangulation_3::Vertex_handle> &list);
#endif


#endif

#endif //#ifndef GEOM2_H
