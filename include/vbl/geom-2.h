#ifndef GEOM2_H
#define GEOM2_H  // header guard
// nuovo header per la geometria, che contiene i metodi per l'inserimento efficiente dell'info block
// si noti che rispetto ai metodi ottenuti dal gruppo di discussione di CGAL, qui e' stato disabilitato il metodo 
// di random shuffling (e' inutile, visto che le cellule sono gia' in posizioni generiche)

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_with_info_3<int, K>		VbI;
typedef CGAL::Alpha_shape_vertex_base_3<K,VbI>					Vb;
typedef CGAL::Alpha_shape_cell_base_3<K>						Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>				Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds>					Dt;
typedef CGAL::Alpha_shape_3<Dt>									Alpha_shape_3;

typedef Dt::Finite_vertices_iterator					Finite_vertices_iterator;
typedef Alpha_shape_3::Vertex_handle					Vertex_handle;
typedef Dt::Point										Point;

// ****************

typedef Dt::Vertex_iterator								Vertex_iterator;
typedef Dt::Point_iterator								Point_iterator;
typedef Dt::Cell_handle									Cell_handle;

typedef Dt::Vertex										Vertex;
// typedef Dh::Cell										Cell;
typedef Dt::Facet										Facet;
typedef Dt::Edge										Edge;
typedef Dt::Finite_edges_iterator						Finite_edges_iterator;
typedef Dt::Finite_facets_iterator						Finite_facets_iterator;
typedef Dt::Finite_cells_iterator						Finite_cells_iterator;
typedef Dt::Cell_circulator								Cell_circulator;
typedef Dt::Facet_circulator							Facet_circulator;

typedef K::Segment_3									Segment;
typedef K::Ray_3										Ray;
typedef K::Object_3										Object;

typedef CGAL::Vector_3<K>								Vector;

//spatial sort traits to use with a pair of point pointers and integer.
template<class Triangulation>
struct Traits_for_spatial_sort:public Triangulation::Geom_traits{
  typedef typename Triangulation::Geom_traits Gt;
  typedef std::pair<const typename Triangulation::Point*,int> Point_3;
  
  struct Less_x_3{
    bool operator()(const Point_3& p,const Point_3& q) const {
      return typename Gt::Less_x_3()(*(p.first),*(q.first));
    }
  };

  struct Less_y_3{
    bool operator()(const Point_3& p,const Point_3& q) const {
      return typename Gt::Less_y_3()(*(p.first),*(q.first));
    }
  };

  struct Less_z_3{
    bool operator()(const Point_3& p,const Point_3& q) const {
      return typename Gt::Less_z_3()(*(p.first),*(q.first));
    }
  };  
  
  Less_x_3  less_x_3_object () const {return Less_x_3();}
  Less_y_3 	less_y_3_object () const {return Less_y_3();}
  Less_z_3 	less_z_3_object () const {return Less_z_3();}
};


//function inserting points into a triangulation
//and setting the info field to the order in the input list.
template <class Triangulation,class Point_iterator>
void
build_triangulation_with_indices(
	Point_iterator begin,
	Point_iterator end,
	Triangulation& T)
{
  std::vector<std::pair<const typename Triangulation::Point*,int> > points;
  int index=-1;
  for (Point_iterator it=begin;it!=end;++it){
    points.push_back(std::make_pair(&(*it),++index));
  }
 // std::random_shuffle (points.begin(), points.end());
  spatial_sort (points.begin(),points.end(), Traits_for_spatial_sort<Triangulation>());

  typename Triangulation::Cell_handle hint;
  for (typename std::vector< std::pair<const typename Triangulation::Point*,int> >::const_iterator
        p = points.begin();p != points.end(); ++p)
  {
    typename Triangulation::Locate_type lt;
    typename Triangulation::Cell_handle c;
    int li, lj;
    c = T.locate (*(p->first), lt, li, lj, hint);

    typename Triangulation::Vertex_handle v = T.insert (*(p->first), lt, c, li, lj);
    if ( v==typename Triangulation::Vertex_handle() )
      hint=c;
    else{
      v->info() = p->second ;
      hint=v->cell();
    }
  }
}  

#endif //#ifndef GEOM2_H
