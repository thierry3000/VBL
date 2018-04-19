//
//file: examples/Interpolation/nn_coordinates_2.C 
//
#include <CGAL/basic.h>
#include <utility>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/natural_neighbor_coordinates_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<float, K>    Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
typedef std::vector< std::pair< Delaunay::Vertex_handle, K::FT  > >                                                      Vertex_handle_vector;

int main()
{
  
  std::vector<std::pair<K::Point_3, float>> V;
  /** this creates a cube of length dist 
   * with associtated value of 1 at zero z level 
   * and associtated value of 2 at z=1
   */
  float dist = 10;
  V.push_back(std::make_pair(K::Point_3(0,0,0), 1));
  V.push_back(std::make_pair(K::Point_3(1*dist,0,0), 1));
  V.push_back(std::make_pair(K::Point_3(0,1*dist,0), 1));
  V.push_back(std::make_pair(K::Point_3(1*dist,1*dist,0), 1));
  
  V.push_back(std::make_pair(K::Point_3(0,0,1*dist), 2));
  V.push_back(std::make_pair(K::Point_3(1*dist,0,1*dist), 2));
  V.push_back(std::make_pair(K::Point_3(0,1*dist,1*dist), 2));
  V.push_back(std::make_pair(K::Point_3(1*dist,1*dist,1*dist), 2));

  Delaunay dt(V.begin(),V.end());

  // this is the querry point
  K::Point_3 p(2, 7, 3);
  // stores the coordinates of the querry point in the new base
  Vertex_handle_vector coords;
  double norm_coeff;
  /**
   * execute sibson_natural neighbors coordinates
   */
  CGAL::Triple<std::back_insert_iterator<Vertex_handle_vector>,K::FT, bool> result = 
    CGAL::sibson_natural_neighbor_coordinates_3(dt, p,std::back_inserter(coords),norm_coeff);
  if(!result.third){
    std::cout << "The coordinate computation was not successful." 
	      << std::endl;
    std::cout << "The point (" <<p << ") lies outside the convex hull."
	      << std::endl;
  }
  
  std::cout << "Coordinate computation successful." << std::endl;
  std::cout << "Normalization factor: " <<norm_coeff << std::endl;  
  
  std::cout << "try Interpolation" << std::endl;
  
  for( auto aCoordinate : coords)
  {
    std::cout << "coordinate in xyz: " << aCoordinate.first->point()  << "coordinate in natural neighbors: " << aCoordinate.second << std::endl;
  }
  std::cout << "Linear combination of natural neighbors with Laplace natural coordinates";
  std::cout << " + correctness (ensured only with an exact number type supporting sqrt)" << std::endl;
  std::cout << is_correct_natural_neighborhood(dt,p,
                                                 coords.begin(),
                                                 coords.end(),
                                                 norm_coeff)
              << std::endl;
  if(!is_correct_natural_neighborhood(dt,p,
                                                 coords.begin(),
                                                 coords.end(),
                                                 norm_coeff))
  {
    std::cout << "bad!!!" << std::endl;
  }
  
  K::FT value_at_p(0);
  for(auto aPoint : coords)
  {
    value_at_p += aPoint.second * aPoint.first->info();
  }
  std::cout << "value at p: " << value_at_p/norm_coeff << std::endl;
  
  return 0; 
}
