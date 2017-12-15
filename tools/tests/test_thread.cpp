// compile with 
//
// g++-mp-6 -fopenmp -frounding-math -O2 -DNDEBUG -I/opt/local/include -I/usr/local/include -L/opt/local/lib -lCGAL -lgmpxx -lmpfr -lgmp -lmpfi -lboost_thread-mt -lm test.cpp -o test
// g++-mp-7 -fopenmp -frounding-math -O2 -DNDEBUG -I/opt/local/include -I/usr/local/include -L/opt/local/lib -lCGAL -lgmpxx -lmpfr -lgmp -lmpfi -lboost_thread-mt -lm test.cpp -o test


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/point_generators_3.h>
#include <vector>
#include <iostream>

#if 1
#ifdef CGAL_LINKED_WITH_TBB
  #define _parallel
#endif
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
#ifdef _parallel
  typedef CGAL::Triangulation_cell_base_3<K> CB;
  typedef CGAL::Triangulation_data_structure_3<Vb,CB,CGAL::Parallel_tag>                    Tds;
#else
  typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
#endif
  
//Use the Fast_location tag. Default or Compact_location works too.
	
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point                                             Point;
	
int main()
{
	  int np=0;
	  
	  std::cout << "Enter number of points: ";
	  std::cin >> np;
	  
	  int nr=0;
	  
	  std::cout << "Enter number of repeats: ";
	  std::cin >> nr;
	  
	  for(int l=0; l < nr; l++)
	  	{
	  
		  std::vector< std::pair<Point,unsigned> > points;
      points.reserve(np);
      CGAL::Random_points_in_cube_3<Point> rnd(1.);
		  for(int k=0; k < np; k++)
			{
// 				double xc = distribution(generator);
// 				double yc = distribution(generator);
// 				double zc = distribution(generator);
				points.push_back( std::make_pair(*rnd++,k) );
				//if( np < 10 ) std::cout << "(x,y,z) = (" << xc << "," << yc << "," << zc << ")" << std::endl;
      }
#ifdef _parallel
      // Construct the locking data-structure, using the bounding-box of the points
      Delaunay::Lock_data_structure locking_ds(CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);
      Delaunay T( points.begin(),points.end(), &locking_ds );
#else
      Delaunay T( points.begin(),points.end() );
#endif
	  
			  // CGAL_assertion( T.number_of_vertices() == 6 );
	  
			  // check that the info was correctly set.
			  Delaunay::Finite_vertices_iterator vit;
			  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
        {
          if( points[ vit->info() ].first != vit->point() )
          {
            std::cerr << "Error different info" << std::endl;
            exit(EXIT_FAILURE);
          }
          //some debug test
          std::cout << vit->info() << std::endl;
        }
			  std::cout << l << "\t OK \r";
			  std::cout.flush();
	  
	  	}
	  
	  return 0;
}
